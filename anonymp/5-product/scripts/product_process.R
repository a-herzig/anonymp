# this script is executed by scripts/product_process.sh
# Compute the expected value of each SNP to share it to User

start_time <- Sys.time()

chunk_size <- 2L^29L

## Load Input

# process arguments
args <- commandArgs(trailingOnly = TRUE)
chunk_name <- args[1]
if (is.na(chunk_name)) {
  chunk_name <- readLines("inbox/1-user-5-product-chunks.txt", n = 1)
}
remove(args)
cat("process chunk", chunk_name, "\n")

ppm_path <- paste("inbox/4-ppm-5-product-chunk", chunk_name, ".Rdata", sep = "")
# ppm_path <- paste("inbox/4-ppm-5-product-chunk", chunk_name, ".bin", sep = "")
reference_path <- paste("inbox/2-reference-5-product-chunk", chunk_name, ".Rdata", sep = "")

# receive ppm and positions from PPM server
fromppm <- readRDS(ppm_path)
ppm <- fromppm$ppm
# ppm <- readBin(ppm_path, "integer", n = chunk_size, size = 2, signed = FALSE)
# receive reference_haplotypes_e5 from GDI
fromreference <- readRDS(reference_path)
reference <- fromreference$reference_haplotypes
ppm_noise <- fromreference$ppm_noise
remove(fromppm, fromreference, ppm_path, reference_path)

## Impute

stopifnot(is.vector(reference))
stopifnot(is.vector(ppm))
stopifnot(is.vector(ppm_noise))
stopifnot(length(reference) == length(ppm_noise))
stopifnot(length(reference) == length(ppm))
imputation <- reference * ppm + ppm_noise

## Share data imputation_e5 to user
user_path <- paste("outbox/5-product-1-user-chunk", chunk_name, ".Rdata", sep = "")
# user_path <- paste("outbox/5-product-1-user-chunk", chunk_name, ".bin", sep = "")
saveRDS(list(imputation_e5 = imputation), user_path)
# writeBin(imputation_e5, user_path, size = 2)

duration <- difftime(Sys.time(), start_time, units = "secs")
duration_line <- paste("product", chunk_name, duration, sep = ",")
duration_path <- "outbox/5-product-duration_record.csv"
write(duration_line, file = duration_path, append = TRUE)
