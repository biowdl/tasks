version 1.0

task MergeCounts {
    input {
        String? preCommand

        Array[File] inputFiles
        String outputFile
        Int featureColumn
        Int valueColumn
        Boolean inputHasHeader
        String featureAttribute = "gene_id"
        File referenceGtf
        Array[String]+? additionalAttributes
    }

    # Based on a script by Szymon Kielbasa/Ioannis Moustakas
    command <<<
        set -e -o pipefail
        mkdir -p ~{sub(outputFile, basename(outputFile) + "$", "")}
        ~{preCommand}
        R --no-save <<CODE
            library(dplyr)
            library(reshape2)
            library(refGenome)

            list.of.files <- c("~{sep='", "' inputFiles}")

            value.i <- ~{valueColumn}
            feature.i <- ~{featureColumn}
            header <- ~{true="TRUE" false="FALSE" inputHasHeader}
            feature.attribute <- "~{featureAttribute}"
            additional.attributes <- c(~{true='"' false="" defined(additionalAttributes)}~{sep='", "' additionalAttributes}~{true='"' false="" defined(additionalAttributes)})
            reference.gtf <- "~{referenceGtf}"
            output.path <- "~{outputFile}"

            d <- do.call(rbind, lapply(list.of.files, function(file){
                d <- read.table(file, sep="\t", header=header, comment.char="#")

                filename <- basename(file)
                colnames(d)[value.i] <- sub("\\\.[^\\\.]*$", "", filename)
                colnames(d)[feature.i] <- "feature"

                d <- d %>% melt(id.vars=feature.i, variable.name="sample",
                    value.name="count")
            }))

            d <- d %>% dcast(feature ~ sample, value.var="count")

            gtf <- ensemblGenome(dirname(reference.gtf))
            read.gtf(gtf, basename(reference.gtf))

            gtf.table <- gtf@ev$gtf
            gtf.table <- gtf.table[order(gtf.table[,feature.attribute]),]
            gtf.table <- gtf.table[!duplicated(gtf.table[,feature.attribute]),]
            id.table <- gtf.table[, c(feature.attribute, additional.attributes), drop=F]
            output.table <- merge(id.table, d, all.y = T, by.y="feature",
                by.x=feature.attribute)

            write.table(output.table, file=output.path, sep="\t", quote=FALSE,
                row.names=FALSE, na="")
        CODE
    >>>

    output {
        File mergedCounts = outputFile
    }

    runtime {
        memory: 4 + (2*length(inputFiles))
        docker: "biowdl/mergecounts:latest"
    }
}