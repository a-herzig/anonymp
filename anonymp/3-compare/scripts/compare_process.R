# this script is executed by scripts/compare_process.sh

start_time <- Sys.time()

## Load Input

# process arguments
args <- commandArgs(trailingOnly = TRUE)
chunk_name <- args[1]
if (is.na(chunk_name)) {
  chunk_name <- readLines("inbox/1-user-3-compare-chunks.txt", n = 1)
}
remove(args)
cat("process chunk", chunk_name, "\n")
preselect_nhaplotype <- 512L

user_path <- paste("inbox/1-user-3-compare-chunk", chunk_name, ".bin.gz", sep = "")
reference_path <- paste("inbox/2-reference-3-compare-chunk", chunk_name, ".bin.gz", sep = "")

genotype_length <- 0x1p10L

# receive genotype_e2 from user
user_file <- gzfile(user_path, "rb")
genotype_e2 <- readBin(user_file, "integer", n = 0x1p10L, size = 1L, signed = FALSE)
close(user_file)
stopifnot(length(genotype_e2) == genotype_length)
# receive reference_haplotypes_e2 from reference
reference_file <- gzfile(reference_path, "rb")
reference_haplotypes_e2 <- readBin(reference_file, "integer", n = 0x1p27L, size = 1L, signed = FALSE)
close(reference_file)
dim(reference_haplotypes_e2) <- c(genotype_length, 0x1p17L)

## Encode genotype

# compare genotype to reference_haplotypes
comparison_e1 <- 1 * sweep(reference_haplotypes_e2, 1L, genotype_e2, "==")

# preselect haplotypes with the better "Same" score
scores <- colSums(comparison_e1)
maximal_score <- sort(scores, decreasing = TRUE)[preselect_nhaplotype]
preselect <- which(scores >= maximal_score)

# hide unselected haplotypes
comparison_e9 <- comparison_e1[, preselect]

## Share data
# share comparison_e9 to PPM server
# comparison matrix size is not constant, writeBin would be difficult

ppm_path <- paste("outbox/3-compare-4-ppm-chunk", chunk_name, ".Rdata", sep = "")
saveRDS(list(comparison_e9 = comparison_e9, preselect_e9 = preselect), ppm_path)

duration <- difftime(Sys.time(), start_time, units = "secs")
duration_line <- paste("compare", chunk_name, duration, sep = ",")
duration_path <- "outbox/3-compare-duration_record.csv"
write(duration_line, file = duration_path, append = TRUE)
