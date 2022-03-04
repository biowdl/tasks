# Copyright (c) 2021 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

version 1.0

task DeconstructSigs {
    input {
        File signaturesMatrix
        File signaturesReference
        String outputPath = "./signatures.rds"

        Int timeMinutes = 15
        String memory = "4G"
        String dockerImage = "quay.io/biocontainers/r-deconstructsigs:1.9.0--r41hdfd78af_1"
    }

    command {
        R --no-echo << EOF
            library(deconstructSigs)
            tumor <- read.table("~{signaturesMatrix}", check.names=F)
            ref <- data.frame(t(read.table("~{signaturesReference}", check.names=F, header=T, row.names="Type")), check.names=F)
            tumor <- tumor[,colnames(ref)]

            sigs <- whichSignatures(tumor.ref=tumor, row.names(tumor), signatures.ref=ref, contexts.needed=T)
            saveRDS(sigs, "~{outputPath}")
        EOF
    }

    output {
        File signatureRDS = outputPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        signaturesMatrix: {description: "A table containing columns represtenting mutation types (matching the types in the signatures reference) and one row with the counts for each of these types for the sample of intrest.",
                           category: "required"}
        signaturesReference: {description: "A table describing the mutational signatures, formatted like those provided by COSMIC.",
                              category: "required"}
        outputPath: {description: "The location the output will be written to.", category: "common"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}