# this script is executed by scripts/reference_process_final.sh
# Tranform reference haplotypes to share them to Imputation Server

start_time <- Sys.time()

library(dqrng)

# chunk_size <- 2L^29L

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
reference_data <- readRDS(reference_path)
reference_haplotypes_e4 <- reference_data$reference_haplotypes
fake_dosage <- reference_data$fake_dosage # = dF
nsnp <- nrow(reference_haplotypes_e4)
nhaplotype_wfake <- ncol(reference_haplotypes_e4)
chunk_size <- length(reference_haplotypes_e4)

# read a 32 bit signed integer to use it as a seed
seed_path <- paste("inbox/4-ppm-2-reference-rand-chunk", chunk_name, ".txt", sep = "")
# set.seed only accepts 32 bit signed integers (from -2 ^ 31 + 1 to 2 ^ 31 - 1)
dqset.seed(as.integer(readLines(seed_path)))
remove(seed_path)
full_shuffle_key <- dqsample.int(chunk_size)
dqset.seed(NULL)
## Encode genotype

# 24/05/2025 : new noise idea from Anthony's script, = B1 in Anthony's script
# this is reverted in 1-User but the ppm interpolation is more difficult to re-order, = B1
ppm_noise <- dqrunif(chunk_size, min = 0.01, max = 1)
ppm_noise <- rep.int(0, chunk_size) # DEBUG : remove noise on imputation value
stopifnot(length(ppm_noise) == chunk_size)
# 24/05/2025 : in row shuffle
# = o1
full_shuffle_key_order <- order(full_shuffle_key)

# final fictive dosage, = dF2
fake_dosage_adjusted <- fake_dosage + rowSums(matrix(ppm_noise[full_shuffle_key_order], nrow = nsnp, ncol = nhaplotype_wfake))
stopifnot(length(fake_dosage_adjusted) == nsnp)
fake_dosage_adjusted <- rep.int(0, nsnp) # DEBUG : there is no dosage effect because no fake values and no noise

# encryption 5
reference_haplotypes_e5 <- c(reference_haplotypes_e4)[full_shuffle_key]
stopifnot(length(reference_haplotypes_e5) == chunk_size)

## Share data

# share reference_haplotypes_e5 to imputation product server
product_path <- paste("outbox/2-reference-5-product-chunk", chunk_name, ".Rdata", sep = "")
stopifnot(is.vector(reference_haplotypes_e5))
stopifnot(length(reference_haplotypes_e5) == chunk_size)
stopifnot(is.vector(ppm_noise))
stopifnot(length(ppm_noise) == chunk_size)
saveRDS(list(reference_haplotypes = reference_haplotypes_e5, ppm_noise = ppm_noise), product_path)
# writeBin(reference_haplotypes_e5, product_path, size = 1L)
user_path <- paste("outbox/2-reference-1-user-chunk", chunk_name, ".Rdata", sep = "")
stopifnot(is.vector(fake_dosage_adjusted))
stopifnot(length(fake_dosage_adjusted) == nsnp)
saveRDS(list(dosage = fake_dosage_adjusted), user_path)

duration <- difftime(Sys.time(), start_time, units = "secs")
duration_line <- paste("reference_final", chunk_name, duration, sep = ",")
duration_path <- "outbox/2-reference-duration_record.csv"
write(duration_line, file = duration_path, append = TRUE)
