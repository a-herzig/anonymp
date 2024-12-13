#!/usr/bin/env Rscript
library(purrr)
library(tidyr)
library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(data.table, warn.conflicts = FALSE)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

compute_score <- function(subtable) {
  return(100 * sum(subtable$succeed) / sum(subtable$unknown))
}

stats_table <- data.table::fread(input_file)
score_table <- stats_table %>%
  dplyr::group_by(target) %>%
  nest() %>%
  dplyr::arrange(target) %>%
  dplyr::mutate(score = sapply(data, compute_score))

pdf(NULL)
ggplot(score_table, aes(factor(target), score)) +
  geom_point() +
  ylim(99.5, 100)
ggsave(output_file)
