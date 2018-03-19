# Author: Ioannis Moustakas, i.moustakas@lumc.nl (Based on a script by Szymon Kielbasa)
# Modified by: Davy Cats, d.cats@lumc.nl
# Title: Merge count files from featureCouns output
# Use: Rscript merge_counts.R columnIDToMergeOn columnIDBeingMerged listOfFilesToBeMerged... > outputFile

### Load Packages
library(dplyr)
library(reshape2)

### load arguments from the command line
args <- commandArgs(trailingOnly=TRUE)
idVars <- args[1]
measureVars <- args[2]
listOfFiles <- args[3:length(args)]

### Iterate over the list of files that are being merged and
### change the column name to the sample name
d <- do.call(rbind, lapply(listOfFiles, function(file){
    d <- read.table(file, header=TRUE, comment.char="#")
    colI <- grep(measureVars, colnames(d))
    colnames(d)[colI] <- strsplit(file, "/")[[1]][3]
    d <- d %>% melt(id.vars=idVars, measure.vars=colI,
                    variable.name="sample", value.name="count")
}))

### Reformat the data frame and output (in STDOUT) the merged table.
d <- d %>% dcast(paste0(idVars, " ~ sample"), value.var="count")
write.table(d, sep="\t", quote=FALSE, row.names=FALSE)
