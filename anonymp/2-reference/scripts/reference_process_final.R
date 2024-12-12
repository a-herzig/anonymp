# this script is executed by scripts/reference_process_final.sh
# Tranform reference haplotypes to share them to Imputation Server

start_time <- Sys.time()

library(dqrng)

chunk_size <- 2L^29L

## Load Input

# process arguments
args <- commandArgs(trailingOnly = TRUE)
chunk_name <- args[1]
if (is.na(chunk_name)) {
  chunk_name <- readLines("inbox/1-user-2-reference-chunks.txt", n = 1)
}
remove(args)
cat("process chunk", chunk_name, "\n")

reference_path <- paste("inbox/2-reference-2-reference-chunk", chunk_name, ".bin.gz", sep = "")
reference_file <- gzfile(reference_path, "rb")
reference_length <- readBin(reference_file, "integer", 1, size = 8)
reference_haplotypes_e4 <- readBin(reference_file, "integer", reference_length, size = 1)
close(reference_file)

## Encode genotype

# encryption 7
cat("prepopulation chunk size : 2 ^", log(length(reference_haplotypes_e4)) / log(2L), "\n")
stopifnot(length(reference_haplotypes_e4) <= chunk_size)
nfake <- chunk_size - length(reference_haplotypes_e4)
fake_ref <- rbinom(nfake, 1L, 0.4)
reference_haplotypes_e7 <- c(reference_haplotypes_e4, fake_ref)

# encryption 5
# read a 32 bit signed integer to use it as a seed
seed_path <- paste("inbox/4-ppm-2-reference-rand-chunk", chunk_name, ".txt", sep = "")
# set.seed only accepts 32 bit signed integers (from -2 ^ 31 + 1 to 2 ^ 31 - 1)
dqset.seed(as.integer(readLines(seed_path)))
remove(seed_path)
full_shuffle_key <- dqsample.int(chunk_size)
reference_haplotypes_e5 <- reference_haplotypes_e7[full_shuffle_key]

## Share data

# share reference_haplotypes_e5 to imputation product server
product_path <- paste("outbox/2-reference-5-product-chunk", chunk_name, ".bin", sep = "")

writeBin(reference_haplotypes_e5, product_path, size = 1L)

duration <- difftime(Sys.time(), start_time, units = "secs")
duration_line <- paste("reference_final", chunk_name, duration, sep = ",")
duration_path <- "outbox/2-reference-duration_record.csv"
write(duration_line, file = duration_path, append = TRUE)
