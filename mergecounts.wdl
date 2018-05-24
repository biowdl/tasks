task MergeCounts {
    String? preCommand

    Array[File] inputFiles
    String outputFile
    Int featureColumn
    Int valueColumn
    Boolean inputHasHeader

    # Based on a script by Szymon Kielbasa/Ioannis Moustakas
    command <<<
        set -e -o pipefail
        mkdir -p ${sub(outputFile, basename(outputFile) + "$", "")}
        ${preCommand}
        R --no-save --slave <<CODE > ${outputFile}
            library(dplyr)
            library(reshape2)

            listOfFiles <- c("${sep='", "' inputFiles}")

            valueI <- ${valueColumn}
            featureI <- ${featureColumn}
            header <- ${true="TRUE" false="FALSE" inputHasHeader}

            d <- do.call(rbind, lapply(listOfFiles, function(file){
                d <- read.table(file, header=header, comment.char="#")

                splitPath <- strsplit(file, "/")[[1]]
                colnames(d)[valueI] <- sub("\\\.[^\\\.]*$", "",
                    splitPath[length(splitPath)])
                colnames(d)[featureI] <- "feature"

                d <- d %>% melt(id.vars=featureI, variable.name="sample", value.name="count")
            }))

            d <- d %>% dcast(feature ~ sample, value.var="count")
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