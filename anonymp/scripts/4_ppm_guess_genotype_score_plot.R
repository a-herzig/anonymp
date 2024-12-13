library(tidyverse)

codes <- readr::read_lines("4-ppm/inbox/1-user-4-ppm-chunks.txt")
nchunk <- length(codes)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) >= 1) {
  output_file <- args[1]
} else {
  output_file <- "tmp/4_ppm_guess_genotype_score_plot.pdf"
}

score_tidy <- data.frame(
  code = rep(codes, 2L),
  chunk = 0L,
  target = 0L,
  tactic = "undefined",
  score = 0
)

for (code in codes) {
  chunk_suffix <- base::paste("-chunk", code, ".Rdata", sep = "")
  cmp <- readRDS(paste("4-ppm/inbox/3-compare-4-ppm", chunk_suffix, sep = ""))
  usr <- readRDS(paste("4-ppm/inbox/1-user-4-ppm", chunk_suffix, sep = ""))
  ref <- readRDS(paste("4-ppm/inbox/2-reference-4-ppm", chunk_suffix, sep = ""))
  # usual 4-ppm input decoding
  nvsnp <- length(usr$visible_snps)
  comparison <- cmp$comparison_e9[order(usr$snp_shuffle_key)[seq_len(nvsnp)], order(ref$haplotype_shuffle_key[cmp$preselect_e9])]
  # 4-ppm try to guess the genotype
  guessed_genotype <- 1 - 1 * (rowSums(comparison) >= nvsnp * 0.5)
  # compare guessed to original
  usr_usr <- readRDS(paste("1-user/inbox/1-user-1-user", chunk_suffix, sep = ""))
  stopifnot(all(usr_usr$visible_snps == usr$visible_snps))
  stopifnot(length(usr_usr$visible_snps) == nvsnp)
  # register scores in a tidy to plot them easily
  genotype <- usr_usr$genotype[usr_usr$visible_snps]
  mask <- score_tidy$code == code
  score_tidy$tactic[mask] <- c("comparison matrix", "all reference nucleotides")
  score_tidy$score[mask] <- c(
    sum(genotype == guessed_genotype) / nvsnp,
    sum(genotype == 0) / nvsnp
  )
  score_tidy$chunk[mask] = usr_usr$chunk_index
  score_tidy$target[mask] = usr_usr$target_index
}

ratio <- sum(score_tidy$score[score_tidy$tactic == "comparison matrix"]) / sum(score_tidy$score[score_tidy$tactic == "all reference nucleotides"])
cat("The ratio comparison-matrix-tactic to all-reference-nucleotides-tactic is ", ratio, ".\n", sep = "")

pdf(NULL)
ggplot(score_tidy, aes(factor(chunk), score, fill = tactic)) +
  geom_col(aes(group = target), position = "dodge") +
  ylim(0, 1) +
  theme(legend.position = "bottom")
ggsave(output_file, width = 256, height = 128, units = "mm")
