# A369 prep: RMA (affy 1.88, hgu133plus2cdf) of the STAGE CEL files (GSE40231_RAW.tar from GEO).
# The deposited series-matrix values carry wrong probe identities (tissue markers do not separate the tissues;
# checked 2026-11-24), so the arrays are normalised here from the raw files. Output: log2 RMA probeset values.
#   Rscript rma_stage.R <raw_dir> <out.tsv.gz>
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages({library(affy); library(hgu133plus2cdf)})
files <- sort(list.files(args[1], pattern = "\\.CEL\\.gz$", full.names = TRUE))
ab <- ReadAffy(filenames = files, compress = TRUE)
stopifnot(cdfName(ab) == "HG-U133_Plus_2")
es <- rma(ab)
x <- exprs(es)
colnames(x) <- sub("_.*$", "", basename(files))
con <- gzfile(args[2], "w")
write.table(data.frame(probe = rownames(x), x, check.names = FALSE), con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
cat("arrays", ncol(x), "probesets", nrow(x), "\n")
