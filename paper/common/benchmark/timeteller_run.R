# TimeTeller (Vlachou et al. 2024; github.com/VadimVasilyev1994/TimeTeller, commit 7c8e7ab) on one train/test split.
#   Rscript timeteller_run.R <train.csv> <test.csv> <out.csv> <cores>
# train.csv: sample, time (h, binned), genes...; test.csv: sample, genes... (log expression).
# Package defaults (intergene normalisation, 3 local PCs, epsilon 0.4, eta 0.35); log_thresh by the package's rule
# (`choose_logthresh` on the training data, then the largest threshold leaving < 2 % of samples with a flat
# likelihood, as `choose_logthresh_plot`). Writes sample, time_pred (1st likelihood peak, h), theta, max_loglik.
suppressPackageStartupMessages(library(TimeTeller))
args <- commandArgs(trailingOnly = TRUE)
train <- read.csv(args[1], check.names = FALSE)
test <- read.csv(args[2], check.names = FALSE)
cores <- as.integer(args[4])
genes <- setdiff(colnames(train), c("sample", "time"))
mat <- function(d) { m <- t(as.matrix(d[, genes])); colnames(m) <- d$sample; m }
quiet <- function(x) suppressMessages(suppressWarnings(capture.output(r <- x, file = nullfile())))

fit <- function(thresh) train_model(exp_matrix = mat(train), genes = genes, group_1 = train$sample, time = train$time,
                                    log_thresh = thresh, parallel_comp = cores > 1, cores = cores)
quiet(model <- fit(-5))
quiet(sel <- choose_logthresh(model, train_or_test = "train"))
ok <- sel$LogThresh[sel$Perc_Flat < 0.02]
thresh <- if (length(ok)) max(ok) else min(sel$LogThresh)
quiet(res <- test_model(model, exp_matrix = mat(test), test_group_1 = test$sample, log_thresh = thresh,
                        parallel_comp = cores > 1, cores = cores))
r <- res$Test_Data$Results_df  # get_tt_res(res, "test"), not exported
write.csv(data.frame(sample = test$sample, time_pred = r$time_1st_peak, theta = r$Theta, max_loglik = r$max_1st_peak,
                     log_thresh = thresh), args[3], row.names = FALSE)
