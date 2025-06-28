# this script is executed by scripts/user_process_final.sh
# Untransform expected value

start_time <- Sys.time()

library(data.table)

## Load Input

# process arguments
args <- commandArgs(trailingOnly = TRUE)
chunk_name <- args[1]
if (is.na(chunk_name)) {
  chunk_name <- readLines("inbox/1-user-1-user-chunks.txt", n = 1)
}
remove(args)
cat("process chunk", chunk_name, "\n")

user_path <- paste("inbox/1-user-1-user-chunk", chunk_name, ".Rdata", sep = "")
reference_path <- paste("inbox/2-reference-1-user-chunk", chunk_name, ".Rdata", sep = "")
ppm_path <- paste("inbox/4-ppm-1-user-chunk", chunk_name, ".Rdata", sep = "")
product_path <- paste("inbox/5-product-1-user-chunk", chunk_name, ".Rdata", sep = "")

fromuser <- readRDS(user_path)
fromreference <- readRDS(reference_path)
fromppm <- readRDS(ppm_path)
fromproduct <- readRDS(product_path)
remove(user_path, reference_path, ppm_path, product_path)

# nhaplotype <- fromreference$nhaplotype
dosage <- fromreference$dosage

# The bellow variables are all relative to the chunk
nsnp <- fromuser$nsnp
minrange <- fromuser$minrange
maxrange <- fromuser$maxrange
positions <- fromuser$positions
target <- fromuser$target
visible_snps <- fromuser$visible_snps
genotype <- fromuser$genotype
target_index <- fromuser$target_index
chunk_index <- fromuser$chunk_index
total_start_time <- fromuser$total_start_time

# ppm_shuffle_key_order <- fromppm$ppm_shuffle_key
# ppm_select_pos <- fromppm$ppm_select_pos
full_shuffle_key_order <- fromppm$full_shuffle_key_order

# receive imputation_e5 from imputation product server
imputation_e5 <- fromproduct$imputation_e5

remove(fromuser, fromreference, fromppm, fromproduct)

## Decode imputation

if (FALSE) {
# remove ppm_shuffle_key shuffle on imputation
imputation_select <- imputation_e5[ppm_shuffle_key_order]
# build a large imputation matrix on all SNPs and all haplotypes
imputation_matrix <- vector("integer", nsnp * nhaplotype)
imputation_matrix[ppm_select_pos] <- imputation_select
dim(imputation_matrix) <- c(nsnp, nhaplotype)
}

# sum SNPs probabilities
stopifnot(length(imputation_e5) %% nsnp == 0)
nhaplotype_wfake <- length(imputation_e5) %/% nsnp
stopifnot(length(dosage) == nsnp)
imputation <- rowSums(matrix(imputation_e5[full_shuffle_key_order], nrow = nsnp, ncol = nhaplotype_wfake)) - dosage
# remove(imputation_select, imputation_matrix)

# which snps are included ? i.e. reported in the target
included_snps <- which(minrange <= positions & positions <= maxrange)

# filter included snps
imputation <- imputation[included_snps]

## Compute stats
# The instructions bellow only operate on included SNP

included_nsnp <- length(included_snps)
included_visible_snps <- match(visible_snps, included_snps)
included_visible_snps <- included_visible_snps[!is.na(included_visible_snps)]
result <- 1 * ((imputation >= 0.5) == target[included_snps])
result[included_visible_snps] <- 0L
is_unknown <- array(1, included_nsnp)
is_unknown[included_visible_snps] <- 0L
score <- 100L * cumsum(result) / cumsum(is_unknown)
summary <- data.frame(
  position = positions[included_snps],
  target = target[included_snps],
  genotype = genotype[included_snps],
  imputed = 1 * (imputation >= 2L^15L),
  compare = 1 * ((imputation >= 2L^15L) == target[included_snps]),
  succeed = cumsum(result),
  unknown = cumsum(is_unknown),
  score = paste(round(score, 2L), "%", sep = "")
)
print(summary[included_nsnp, ])
summary_path <- paste("outbox/1-user-summary-chunk", chunk_name, ".csv", sep = "")
fwrite(summary, summary_path)

stats <- data.frame(
  code = chunk_name,
  target = target_index,
  chunk = chunk_index,
  start = total_start_time,
  stop = Sys.time(),
  succeed = sum(result),
  unknown = sum(is_unknown)
)
fwrite(stats, "outbox/1-user-stats.csv", append = TRUE)

duration <- difftime(Sys.time(), start_time, units = "secs")
duration_line <- paste("user_final", chunk_name, duration, sep = ",")
duration_path <- "outbox/1-user-duration_record.csv"
write(duration_line, file = duration_path, append = TRUE)
