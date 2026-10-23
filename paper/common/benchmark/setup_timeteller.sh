#!/usr/bin/env bash
# R 4.3 environment for TimeTeller (github.com/VadimVasilyev1994/TimeTeller, commit 7c8e7ab; needs R >= 4.2, the
# system R is 4.1). cosinor2 from CRAN, pkgndep 1.99.3 from the CRAN archive (later versions no longer export
# check_pkg, which TimeTeller imports). The user's ~/.Renviron is ignored so nothing is read from or written to ~/R.
set -e
T=${SCCC_DATA:?set SCCC_DATA}/tools
[ -d $T/TimeTeller ] || git clone https://github.com/VadimVasilyev1994/TimeTeller $T/TimeTeller
git -C $T/TimeTeller checkout -q 7c8e7ab
L=$T/timeteller_env/lib/R/library
export R_ENVIRON_USER=/dev/null R_LIBS_USER=$L
conda create -y -p $T/timeteller_env -c conda-forge -c bioconda --override-channels r-base=4.3 compilers make \
  r-circular r-compositions r-data.table r-doparallel r-dplyr r-foreach r-ggplot2 r-ggrepel r-matrix r-mvtnorm \
  r-plotly r-pracma r-psych r-purrr r-rfast r-rgl r-tibble r-tidyr r-remotes r-cosinor r-matrixstats r-gridextra \
  r-hmisc r-getoptlong r-htmltable r-hash r-biocmanager r-brew r-cowplot r-globaloptions r-formula \
  bioconductor-complexheatmap bioconductor-biocversion
$T/timeteller_env/bin/Rscript -e 'install.packages("cosinor2", lib = Sys.getenv("R_LIBS_USER"), repos = "https://cloud.r-project.org")'
curl -sO https://cran.r-project.org/src/contrib/Archive/pkgndep/pkgndep_1.99.3.tar.gz
$T/timeteller_env/bin/R CMD INSTALL -l $L pkgndep_1.99.3.tar.gz && rm pkgndep_1.99.3.tar.gz
$T/timeteller_env/bin/R CMD INSTALL -l $L $T/TimeTeller
$T/timeteller_env/bin/Rscript -e 'library(TimeTeller); cat("TimeTeller OK\n")'
