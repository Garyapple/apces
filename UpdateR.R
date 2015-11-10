# installing/loading the package:
if(!require(installr)) { install.packages("installr"); require(installr)} #load / install+load installr

# step by step functions:
check.for.updates.R() # tells you if there is a new version of R or not.
install.R() # download and run the latest R installer
copy.packages.between.libraries() # copy your packages to the newest R installation from the one version before it (if ask=T, it will ask you between which two versions to perform the copying)
update.packages(checkBuilt=TRUE, ask=FALSE)

source("http://bioconductor.org/biocLite.R")
biocLite()   