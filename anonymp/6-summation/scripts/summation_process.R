# this script is executed by scripts/summation_process.sh
# Untransform expected value

start_time <- Sys.time()

library(data.table)

## Load Input

# process arguments
args <- commandArgs(trailingOnly = TRUE)
chunk_name <- args[1]
if (is.na(chunk_name)) {
  chunk_name <- readLines("inbox/1-user-6-summation-chunks.txt", n = 1)
}
remove(args)
cat("process chunk", chunk_name, "\n")

user_path <- paste("inbox/1-user-6-summation-chunk", chunk_name, ".Rdata", sep = "")
ppm_path <- paste("inbox/4-ppm-6-summation-chunk", chunk_name, ".Rdata", sep = "")
product_path <- paste("inbox/5-product-6-summation-chunk", chunk_name, ".Rdata", sep = "")

fromuser <- readRDS(user_path)
fromppm <- readRDS(ppm_path)
fromproduct <- readRDS(product_path)
remove(user_path, ppm_path, product_path)


# The bellow variables are all relative to the chunk
nsnp <- fromuser$nsnp


full_shuffle_key_order <- fromppm$full_shuffle_key_order

# receive imputation_e5 from imputation product server
imputation_e5 <- fromproduct$imputation_e5

remove(fromuser, fromppm, fromproduct)

## Decode imputation

# sum SNPs probabilities
stopifnot(length(imputation_e5) %% nsnp == 0)
nhaplotype_wfake <- length(imputation_e5) %/% nsnp

imputation <- rowSums(matrix(imputation_e5[full_shuffle_key_order], nrow = nsnp, ncol = nhaplotype_wfake))

summation_path <- paste("outbox/6-summation-1-user-chunk", chunk_name, ".Rdata", sep = "")
saveRDS(list(imputation), summation_path)

duration <- difftime(Sys.time(), start_time, units = "secs")
duration_line <- paste("summation", chunk_name, duration, sep = ",")
duration_path <- "outbox/6-summation_process_record.csv"
write(duration_line, file = duration_path, append = TRUE)
