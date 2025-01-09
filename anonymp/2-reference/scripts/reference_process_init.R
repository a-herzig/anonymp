# this script is executed by scripts/reference_process_init.sh

start_time <- Sys.time()

library(data.table)
library(dqrng)
library(R.utils)

POSITIONS_PATH <- "tmp/positions.txt"
REFERENCE_PATH <- "res/chr15_5popSim_4B11_Ref.vcf.gz"

## Load Input

# process arguments
args <- commandArgs(trailingOnly = TRUE)
chunk_name <- args[1]
if (is.na(chunk_name)) {
  chunk_name <- readLines("inbox/1-user-2-reference-chunks.txt", n = 1)
}
remove(args)
cat("process chunk", chunk_name, "\n")

# read keys from user
user_path <- paste("inbox/1-user-2-reference-chunk", chunk_name, ".Rdata", sep = "")
fromuser <- readRDS(user_path)

# read chunk positions
positions <- fread(POSITIONS_PATH, header = FALSE)
positions_minindex <- min(which(fromuser$minbound <= positions))
positions_maxindex <- max(which(fromuser$maxbound >= positions))
remove(positions)
nsnp <- positions_maxindex - positions_minindex + 1L

# read reference table
reference_haplotypes <- fread(
  text = gsub("|", "\t", readLines(REFERENCE_PATH, 10L + positions_maxindex)[(10L + positions_minindex):(10L + positions_maxindex)], fixed = TRUE, useBytes = TRUE),
  header = FALSE, sep = "\t",
  colClasses = c(rep("NULL", 9L), rep("integer", 90000L))
)
reference_haplotypes <- unname(as.matrix(reference_haplotypes))
nhaplotype <- ncol(reference_haplotypes)

## Encode reference

genotype_length <- 0x1p10L
haplotype_length <- 0x1p17L

# not reverted shuffle halotypes
haplotype_shuffle_key <- dqsample.int(nhaplotype, replace = FALSE)
reference_haplotypes_e4 <- reference_haplotypes[, haplotype_shuffle_key]
remove(haplotype_shuffle_key)

# filter visible snps
reference_haplotypes_e3 <- reference_haplotypes_e4[fromuser$visible_snps, ]

# add fake snps
# write only 0 in fake snps to don't favor some haplotypes at the cost of others because this could lead to non-optimal preselect choices
reference_haplotypes_e9 <- matrix(0L, genotype_length, ncol(reference_haplotypes_e3))
reference_haplotypes_e9[seq_along(fromuser$visible_snps), ] <- reference_haplotypes_e3

# add fake haplotypes
nfake_haplotype <- haplotype_length - ncol(reference_haplotypes)
fake_haplotypes <- matrix(dqsample.int(2L, nfake_haplotype * genotype_length, replace = TRUE), genotype_length)
reference_haplotypes_e6 <- cbind(reference_haplotypes_e9, fake_haplotypes)

# reverted shuffle haplotypes
haplotype_shuffle_key <- dqsample.int(haplotype_length, replace = FALSE)
reference_haplotypes_e8 <- reference_haplotypes_e6[, haplotype_shuffle_key]

# encryption 1
stopifnot(length(fromuser$snp_shuffle_key) == genotype_length)
stopifnot(nrow(reference_haplotypes_e8) == genotype_length)
reference_haplotypes_e1 <- reference_haplotypes_e8[fromuser$snp_shuffle_key, ]

# encryption 2
reference_haplotypes_e2 <- sweep(reference_haplotypes_e1, 1L, fromuser$snp_noise_key, "+") %% 2L
stopifnot(nrow(reference_haplotypes_e2) == genotype_length)
stopifnot(ncol(reference_haplotypes_e2) == haplotype_length)

## Share data

user_path <- paste("outbox/2-reference-1-user-chunk", chunk_name, ".Rdata", sep = "")
reference_path <- paste("outbox/2-reference-2-reference-chunk", chunk_name, ".bin.gz", sep = "")
compare_path <- paste("outbox/2-reference-3-compare-chunk", chunk_name, ".bin.gz", sep = "")
ppm_path <- paste("outbox/2-reference-4-ppm-chunk", chunk_name, ".Rdata", sep = "")

# share nhaplotype to user for final sum
saveRDS(list(nhaplotype = nhaplotype), user_path)

# keep reference_haplotypes_e4 for final step
reference_file <- gzfile(reference_path, "wb", compression = 3L)
writeBin(length(reference_haplotypes_e4), reference_file, size = 8L)
writeBin(as.vector(reference_haplotypes_e4), reference_file, size = 1L)
close(reference_file)

# share reference_haplotypes_e2 to compare server
compare_file <- gzfile(compare_path, "wb", compression = 3L)
writeBin(as.vector(reference_haplotypes_e2), compare_file, size = 1L)
close(compare_file)

#
saveRDS(list(haplotype_shuffle_key = unname(haplotype_shuffle_key), nhaplotype = nhaplotype), ppm_path)

duration <- difftime(Sys.time(), start_time, units = "secs")
duration_line <- paste("reference_init", chunk_name, duration, sep = ",")
duration_path <- "outbox/2-reference-duration_record.csv"
write(duration_line, file = duration_path, append = TRUE)
