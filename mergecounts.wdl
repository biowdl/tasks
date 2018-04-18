task MergeCounts {
    String? preCommand

    Array[File] inputFiles
    String outputFile
    String idVar
    String measurementVar

    # Based on a script by Szymon Kielbasa/Ioannis Moustakas
    command <<<
        set -e -o pipefail
        ${preCommand}
        R --no-save --slave <<CODE > ${outputFile}
            library(dplyr)
            library(reshape2)

            listOfFiles <- c("${sep='", "' inputFiles}")

            d <- do.call(rbind, lapply(listOfFiles, function(file){
                d <- read.table(file, header=TRUE, comment.char="#")
                colI <- grep(${measurementVar}, colnames(d))
                colnames(d)[colI] <- strsplit(file, "/")[[1]][3]
                d <- d %>% melt(id.vars=${idVar}, measure.vars=colI,
                    variable.name="sample", value.name="count")
            }))

            d <- d %>% dcast(paste0(${idVar}, " ~ sample"), value.var="count")
            write.table(d, sep="\t", quote=FALSE, row.names=FALSE)
        CODE
    >>>

    output {
        File mergedCounts = outputFile
    }

    runtime {
        memory: 4 + (2*length(inputFiles))
    }
}