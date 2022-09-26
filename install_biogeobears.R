# then install BioGeoBEARS from source:
install.packages("Rcpp", dependencies=TRUE)
install.packages("RcppArmadillo", dependencies=TRUE)
install.packages("gdata", dependencies=TRUE)
install.packages("gtools", dependencies=TRUE)
install.packages("xtable", dependencies=TRUE)
install.packages("plotrix", dependencies=TRUE)
install.packages("vegan", dependencies=TRUE)
install.packages("FD", dependencies=TRUE)
install.packages("SparseM", dependencies=TRUE)
install.packages("ape", dependencies=TRUE)
install.packages("phylobase", dependencies=TRUE)
install.packages("rexpokit", dependencies=TRUE)
install.packages("cladoRcpp", dependencies=TRUE)


# From October 2018 onwards, install BioGeoBEARS from GitHub:
# https://github.com/nmatzke/BioGeoBEARS
install.packages("devtools")
library(devtools)
devtools::install_github(repo="nmatzke/BioGeoBEARS", INSTALL_opts="--byte-compile", dependencies=F)
