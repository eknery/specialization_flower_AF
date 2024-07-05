
# From October 2018 onwards, install BioGeoBEARS from GitHub:
# https://github.com/nmatzke/BioGeoBEARS
install.packages("devtools")
library("devtools")
devtools::install_github(repo="nmatzke/BioGeoBEARS", 
                         INSTALL_opts="--byte-compile", 
                         dependencies=T)
