# this script is executed by scripts/user_process_init.sh
# share an encoded genotype to compare server
# share encode keys to reference server

start_time <- Sys.time()

library(data.table)
library(dqrng)
library(R.utils)

CHUNKS_PATH <- "res/coordinates.I51.2_15.txt"
TARGET_PATH <- "res/chr15_5popSim_4B11_Cases_FULL.vcf.gz"
VISIBLE_POS_PATH <- "res/chip.txt"
share_dists_to_4_ppm <- TRUE

## Load Input

# process arguments
args <- commandArgs(trailingOnly = TRUE)
chunk_name <- args[1]
target_index <- as.integer(args[2])
chunk_code <- as.integer(args[3])
remove(args)
cat("process", chunk_name, "target =", target_index, "chunk =", chunk_code, "\n")

# read chunk positions
chunks_data <- fread(CHUNKS_PATH, sep = "\t", header = FALSE, colClasses = c("integer", "NULL", "character", "character", "NULL", "NULL", "NULL"), col.names = c("names", "bounds", "range"))
chunk_data <- chunks_data[chunks_data$names == chunk_code]
bounds <- as.integer(strsplit(strsplit(chunk_data$bounds, ":")[[1L]][2L], "-")[[1L]])
minbound <- bounds[1L]
maxbound <- bounds[2L]
range <- as.integer(strsplit(strsplit(chunk_data$range, ":")[[1L]][2L], "-")[[1L]])
minrange <- range[1L]
maxrange <- range[2L]
remove(chunks_data, chunk_data, bounds, range)

# read target table
# col 1 = pos, other cols = haplotypes
table <- fread(TARGET_PATH, sep = "\t", header = FALSE, skip = 10L, colClasses = c("NULL", "integer", rep("NULL", 6L + target_index), "character", rep("NULL", 2000L - target_index)))
positions <- table[, 1L]
snp_index <- which(minbound <= positions & positions <= maxbound)
positions <- unlist(unname(table[snp_index, 1L]))
target <- as.integer(sapply(strsplit(unlist(table[snp_index, 2L]), "|"), function(v) v[1L]))
remove(table)

nsnp <- length(snp_index)

# read visible snps in chip.txt
visible_snps <- fread(VISIBLE_POS_PATH, sep = "\t", colClasses = c("NULL", "integer"))
# compute snp specific data to share it
visible_snps <- match(unlist(visible_snps), positions)
visible_snps <- visible_snps[!is.na(visible_snps)]
nvsnp <- length(visible_snps)

if(share_dists_to_4_ppm) {
  # compute ponderated distance
  stopifnot(length(positions) == nsnp)
  visible_positions <- positions[visible_snps]
  visible_distances <- visible_positions[2L:nvsnp] - visible_positions[1L:(nvsnp - 1L)]
  distances <- positions[2L:nsnp] - positions[1L:(nsnp - 1L)]
  visible_pond_dists <- as.integer(round(0x1p24 * nvsnp / sum(visible_distances) * visible_distances))
  pond_dists <- as.integer(round(0x1p24 * nsnp / sum(distances) * distances))
  remove(visible_positions, visible_distances, distances)
} else {
  visible_pond_dists <- rep.int(0x1p24L, nvsnp - 1L)
  pond_dists <- rep.int(0x1p24L, nsnp - 1L)
}
stopifnot(length(visible_pond_dists) == nvsnp - 1L)
stopifnot(length(pond_dists) == nsnp - 1L)
stopifnot(round(mean(visible_pond_dists), -6) == round(0x1p24L, -6))
stopifnot(round(mean(pond_dists), -6) == round(0x1p24L, -6))

# compute genotype using visible_snp_index and target
genotype <- array(NA, dim = nsnp)
genotype[visible_snps] <- target[visible_snps]

## Encode genotype
genotype_length <- 1024L
stopifnot(length(visible_snps) <= genotype_length)

# encryption 3
genotype_e3 <- vector("integer", length = genotype_length)
genotype_e3[seq_along(visible_snps)] <- target[visible_snps]

# encryption 1
snp_shuffle_key <- dqsample.int(genotype_length, replace = FALSE)
genotype_e1 <- genotype_e3[snp_shuffle_key]

# encryption 2
snp_noise_key <- dqsample.int(2L, genotype_length, replace = TRUE)
genotype_e2 <- (genotype_e1 + snp_noise_key) %% 2L

## Share data

user_path <- paste("outbox/1-user-1-user-chunk", chunk_name, ".Rdata", sep = "")
reference_path <- paste("outbox/1-user-2-reference-chunk", chunk_name, ".Rdata", sep = "")
compare_path <- paste("outbox/1-user-3-compare-chunk", chunk_name, ".bin.gz", sep = "")
ppm_path <- paste("outbox/1-user-4-ppm-chunk", chunk_name, ".Rdata", sep = "")

datatouser <- list(
  nsnp = nsnp,
  target_index = target_index,
  chunk_index = chunk_code,
  total_start_time = start_time,
  minrange = minrange,
  maxrange = maxrange,
  visible_snps = visible_snps,
  target = target,
  genotype = genotype,
  positions = positions
)
saveRDS(datatouser, file = user_path)

# share visible_snps & snp_shuffle_key & snp_noise_key to GDI
datatoreference <- list(
  minbound = minbound,
  maxbound = maxbound,
  snp_shuffle_key = unname(snp_shuffle_key),
  snp_noise_key = unname(snp_noise_key),
  visible_snps = unname(visible_snps)
)
saveRDS(datatoreference, file = reference_path)

# share genotype_e2 to compare server
compare_file <- gzfile(compare_path, "wb", compression = 3L)
writeBin(genotype_e2, compare_file, size = 1L)
close(compare_file)

# share visible_snps & snp_shuffle_key to PPM server
datatoppm <- list(
  visible_snps = visible_snps,
  nsnp = nsnp,
  snp_shuffle_key = unname(snp_shuffle_key),
  pond_dists = unname(pond_dists),
  visible_pond_dists = unname(visible_pond_dists)
)
saveRDS(datatoppm, file = ppm_path)

duration <- difftime(Sys.time(), start_time, units = "secs")
duration_line <- paste("user_init", chunk_name, duration, sep = ",")
duration_path <- "outbox/1-user-duration_record.csv"
write(duration_line, file = duration_path, append = TRUE)
