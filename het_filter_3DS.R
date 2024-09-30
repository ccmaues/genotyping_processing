#!/usr/bin/env Rscript
library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
dat <- fread(input_file, data.table = FALSE)
valid <- dat[dat$F <= mean(dat$F) + 3 * sd(dat$F) & dat$F >= mean(dat$F) - 3 * sd(dat$F), ]
fwrite(valid[, c("FID", "IID")], "10_het.valid.sample", sep = "\t")
quit()
