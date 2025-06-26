# Compute PPM to share it to Imputation Server

start_time <- Sys.time()

library(matrixStats)
library(dqrng)

RECOMBINATION_RATE <- 1 / 256
ERROR_RATE <- 1 / 256

## Load Input

# process arguments
args <- commandArgs(trailingOnly = TRUE)
chunk_name <- args[1]
if (is.na(chunk_name)) {
  chunk_name <- readLines("inbox/1-user-4-ppm-chunks.txt", n = 1)
}
remove(args)
cat("process chunk", chunk_name, "\n")

# receive visible_snps and snp_shuffle_key from user
user_path <- paste("inbox/1-user-4-ppm-chunk", chunk_name, ".Rdata", sep = "")
fromuser <- readRDS(user_path)
snp_shuffle_key <- fromuser$snp_shuffle_key
visible_snps <- fromuser$visible_snps
nsnp <- fromuser$nsnp
logvisible_pond_dists <- log(fromuser$visible_pond_dists) - 24 * log(2)
logpond_dists <- log(fromuser$pond_dists) - 24 * log(2)
stopifnot(round(mean(logvisible_pond_dists), -1) == 0)
stopifnot(round(mean(logpond_dists), -1) == 0)

# read keys from reference
reference_path <- paste("inbox/2-reference-4-ppm-chunk", chunk_name, ".Rdata", sep = "")
fromreference <- readRDS(reference_path)
haplotype_shuffle_key <- fromreference$haplotype_shuffle_key
nhaplotype <- fromreference$nhaplotype
fake_ppm_matrix <- fromreference$fake_ppm_matrix

# read comparison table
compare_path <- paste("inbox/3-compare-4-ppm-chunk", chunk_name, ".Rdata", sep = "")
fromcompare <- readRDS(compare_path)
comparison_e9 <- fromcompare$comparison_e9
preselect_e9 <- fromcompare$preselect_e9

remove(user_path, fromuser, reference_path, fromreference, compare_path, fromcompare)

## Compute PPM

#' Compute Posterior Probability Matrix using a dynamic Hidden Markov Model
hmm_ppm <- function(target_compare) {
  # genetic imputation using simple HMM package

  # target_compare is a matrix where col=haplotype, row=snp
  # TRUE -> same allele between target and haplotype on this snp
  # FALSE -> different allele between target and haplotype on this snp

  # recombination_rate is the probability of recombination between two adjacent snps
  # error_rate is the probability to have an error between T/F in target_compare for each value
  nhaplotype <- ncol(target_compare)
  nsnp <- nrow(target_compare)

  # implementation follows "A revealing Introduction to Hidden Markov Models"
  # HMM model
  logStartProb <- -log(nhaplotype)

  # recombination rate
  # expected number of recombination between t and t-1
  lambda <- 1 - exp(-RECOMBINATION_RATE)
  log_trans_probs <- log(diag(nhaplotype) * (1 - lambda) + lambda / nhaplotype)
  stopifnot(round(rowSums(exp(log_trans_probs)), 12L) == 1)

  logEmissionProbs <- log(c(1 - ERROR_RATE, ERROR_RATE))
  stopifnot(round(sum(exp(logEmissionProbs)), 12) == 1)

  # observations
  observations <- 2L - 1L * target_compare

  # convert comparisons to probas
  # compute forward probability
  alpha <- matrix(NA, nrow = nsnp, ncol = nhaplotype)
  alpha[1, ] <- logStartProb + logEmissionProbs[observations[1, ]]

  # alpha_{t, i} = sum_j(alpha_{t - 1, j} * trans_prob_{j, i}) * emission_probs_{t, i}
  # equivalent log equation : log_alpha_{t, i} = logsumexp_j(log_alpha_{t - 1, j} + log_trans_prob_{j, i}) + log_emission_probs_{t, i}
  # (log : sum -> logsumexp; prod -> sum)
  # R vectorisation : replace loops by operations on matrices
  for (t in 2:nsnp) {
    # compute alpha for SNP t, depend on t - 1
    temp <- sweep(log_trans_probs + logvisible_pond_dists[t - 1], 1L, alpha[t - 1L, ], "+")
    alpha[t, ] <- colLogSumExps(temp) + logEmissionProbs[observations[t, ]]
  }

  # compute backward probability
  beta <- matrix(NA, nrow = nsnp, ncol = nhaplotype)
  beta[nsnp, ] <- 0L

  # beta_{t, i} = sum_j(trans_probs_{i, j} * emission_probs_{t + 1, j} * beta_{t + 1, j})
  # <=> log_beta_{t, i} = logsumexp_j(log_trans_probs_{i, j} * log_emission_probs_{t + 1, j} + log_beta_{t + 1, j})
  for (t in (nsnp - 1):1) {
    # compute beta for SNP t, depend on t + 1
    temp <- logEmissionProbs[observations[t + 1L, ]] + beta[t + 1L, ]
    beta[t, ] <- colLogSumExps(sweep(log_trans_probs + logvisible_pond_dists[t], 1L, temp, "+"))
  }

  logppm <- alpha + beta
  logppm <- sweep(logppm, 1L, rowLogSumExps(logppm), "-")
  stopifnot(round(rowSums(exp(logppm)), digits = 12) == 1)
  return(logppm)
}

#' PPM interpolation
interpolate_ppm <- function(raw_logppm, visible_snps, nsnp) {
  nhaplotype <- ncol(raw_logppm)

  interpolated_logppm <- matrix(0, nrow = nsnp, ncol = nhaplotype)
  for (i in seq_len(nsnp)) {
    if (i %in% visible_snps) {
      interpolated_logppm[i, ] <- raw_logppm[which(visible_snps == i), ]
    } else if (i <= min(visible_snps)) {
      interpolated_logppm[i, ] <- raw_logppm[1, ]
    } else if (i >= max(visible_snps)) {
      interpolated_logppm[i, ] <- raw_logppm[length(visible_snps), ]
    } else {
      jinf <- which(visible_snps <= i)
      if (length(jinf) > 1L) {
        jinf <- max(jinf)
      }
      jsup <- which(visible_snps >= i)
      if (length(jsup) > 1L) {
        jsup <- min(jsup)
      }
      j <- c(jinf, jsup)
      logp <- raw_logppm[j, ]
      logw <- log(rev(abs(visible_snps[j] - i)))
      # add distance to weight
      snp_side <- cut(seq_len(nsnp), breaks = c(i, visible_snps[j]), labels = c("inf", "sup"), right = FALSE)
      logw[1] <- logw[1] - sum(logpond_dists[snp_side == "inf"], na.rm = TRUE)
      logw[2] <- logw[2] - sum(logpond_dists[snp_side == "sup"], na.rm = TRUE)
      interpolated_logppm[i, ] <- colLogSumExps(logp + logw) - logSumExp(logw)
    }
  }
  return(interpolated_logppm)
}

# decode comparison_e1
comparison_e8 <- comparison_e9[, order(haplotype_shuffle_key[preselect_e9])]
comparison <- comparison_e8[order(snp_shuffle_key)[seq_along(visible_snps)], ]

# which are the preselect haplotypes
preselect <- which(order(haplotype_shuffle_key) %in% preselect_e9)
# check no fake haplotype has been selected by 3-compare
# if it happens, restart the imputation on this chunk later
# if we permit that to happen, final result is inconstant; remove fake haplotypes, they are at the end
# comparison <- comparison[, preselect <= nhaplotype]
stopifnot(max(preselect) <= nhaplotype)

# compute raw ppm
# raw ppm on preselect only
raw_logppm <- hmm_ppm(comparison)

# interpolate logppm
interpolated_logppm <- na.fail(interpolate_ppm(raw_logppm, visible_snps, nsnp))
stopifnot(nrow(interpolated_logppm) == nsnp)
stopifnot(ncol(interpolated_logppm) == length(preselect))

# position vector has to point out of preselect referential
# because the reference_haplotype matrix of the product will
# not be filtered in preselect
standard_ppm <- matrix(0L, nrow = nsnp, ncol = nhaplotype)
standard_ppm[, preselect] <- exp(interpolated_logppm)

# 24/05/2025 : integrate fake PPM from 2-reference
stopifnot(nrow(fake_ppm_matrix) == nsnp)
stopifnot(ncol(fake_ppm_matrix) == 0x1p10L)
ppm_wfake <- cbind(standard_ppm, fake_ppm_matrix)
ppm_wfake <- standard_ppm # DEBUG do not add fake ppm matrix
stopifnot(nrow(ppm_wfake) == nsnp)
# stopifnot(ncol(ppm_wfake) == nhaplotype + 0x1p10L) DEBUG fake ppm matrix is not included
chunk_size <- length(ppm_wfake)
if (FALSE) { # multi-line comment on optimisations to upgrade the program easily

# filter the ppm to remove negligible probabilities
# ppm is now represented with a pair of vectors (value, position)
ppm_select <- integer_ppm[integer_ppm >= 1L]
ppm_select_pos <- which(integer_ppm >= 1L)
stopifnot(length(ppm_select) == length(ppm_select_pos))

# add fake probas to have all chunks of the same size
chunk_size <- 2L^29L
schunk_size <- 2L^21L
cat("prepopulation chunk size : 2 ^", log(length(integer_ppm)) / log(2L), "\n")
stopifnot(length(integer_ppm) <= chunk_size)
cat("prepopulation chunk size nonnull : 2 ^", log(length(ppm_select)) / log(2L), "\n")
stopifnot(length(ppm_select) <= schunk_size)

nfake <- schunk_size - length(ppm_select)
fake_ppm <- dqsample.int(2L^16L, nfake, replace = TRUE)
fake_ppm_pos <- dqsample(setdiff(seq_len(chunk_size), ppm_select_pos), nfake, replace = FALSE)
stopifnot(length(fake_ppm) == nfake)
stopifnot(length(fake_ppm_pos) == nfake)
ppm_wfake <- c(ppm_select, fake_ppm)
ppm_wfake_pos <- c(ppm_select_pos, fake_ppm_pos)
stopifnot(length(ppm_wfake) == length(ppm_wfake_pos))
stopifnot(length(ppm_wfake) == schunk_size)
stopifnot(max(ppm_wfake_pos) <= chunk_size)
}

# read a 32 bit signed integer to use it as a seed
seed_path <- paste("tmp/rand-chunk", chunk_name, ".txt", sep = "")
# set.seed only accepts 32 bit signed integers (from -2 ^ 31 + 1 to 2 ^ 31 - 1)
dqset.seed(as.integer(readLines(seed_path)))
remove(seed_path)

# encryption 5
# move positions of vectors
# create a second key reproducible using seed shared to 2-Reference, = ro in Anthony's script
full_shuffle_key <- dqsample.int(chunk_size, replace = FALSE)
dqset.seed(seed = NULL) # the following random operations are not seed-predicted
# update position vector
ppm_wfake_shuffled <- c(ppm_wfake)[full_shuffle_key]

nhaplotype_wfake <- nhaplotype + 0x1p10L
nhaplotype_wfake <- nhaplotype # DEBUG keep number of haplotypes without the fake ones
stopifnot(length(ppm_wfake_shuffled) == nsnp * nhaplotype_wfake)
ppm_shared <- ppm_wfake_shuffled

# shuffle PPM vectors
if (FALSE) {
ppm_shuffle_key <- dqsample.int(schunk_size, replace = FALSE)
ppm_shared <- ppm_wfake[ppm_shuffle_key]
ppm_shared_pos <- ppm_wfake_pos_shuffled[ppm_shuffle_key]
}

# 24/05/2025 : in row shuffle
# = o1
full_shuffle_key_order <- order(full_shuffle_key)
# this is not reverted in 1-User final script, but has no consequence on the final result because the sum in rowSums operation is commutative, = o2 in Anthony's script
full_shuffle_key_order_w_row_shuffle <- c(t(apply(matrix(full_shuffle_key_order, nrow = nsnp, ncol = nhaplotype_wfake), 1, sample, nhaplotype_wfake, FALSE)))
full_shuffle_key_order_w_row_shuffle <- full_shuffle_key_order # DEBUG do not use the shuffled order in rows
## Share data
user_path <- paste("outbox/4-ppm-1-user-chunk", chunk_name, ".Rdata", sep = "")
product_path <- paste("outbox/4-ppm-5-product-chunk", chunk_name, ".Rdata", sep = "")

# share preselect to user
stopifnot(is.vector(full_shuffle_key_order_w_row_shuffle))
stopifnot(length(full_shuffle_key_order_w_row_shuffle) == chunk_size)
saveRDS(list(full_shuffle_key_order = full_shuffle_key_order_w_row_shuffle), user_path)
# user and reference will receive seed to compute full_shuffle_key
# share interpolated_ppm_e5 to imputation product server
stopifnot(is.vector(ppm_shared))
stopifnot(length(ppm_shared) == chunk_size)
saveRDS(list(ppm = ppm_shared), product_path)
# write ppm in uint8 binary format
# reader needs to know the size before reading, it will be constant

duration <- difftime(Sys.time(), start_time, units = "secs")
duration_line <- paste("ppm", chunk_name, duration, sep = ",")
duration_path <- "outbox/4-ppm-duration_record.csv"
write(duration_line, file = duration_path, append = TRUE)
