version 1.0

# Copyright (c) 2017 Leiden University Medical Center
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

task HTSeqCount {
    input {
        Array[File]+ inputBams
        File gtfFile
        String outputTable = "output.tsv"
        String order = "pos"
        String stranded = "no"
        Array[String] additionalAttributes = []

        String? featureType
        String? idattr

        Int nprocesses = 1
        String memory = "8GiB"
        Int timeMinutes = 1440 #10 + ceil(size(inputBams, "GiB") * 60) FIXME
        String dockerImage = "quay.io/biocontainers/htseq:0.12.4--py37hb3f55d8_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputTable})"
        htseq-count \
        --nprocesses ~{nprocesses} \
        -r ~{order} \
        -s ~{stranded} \
        ~{"--type " + featureType} \
        ~{"--idattr " + idattr} \
        ~{true="--additional-attr " false="" length(additionalAttributes) > 0 }~{sep=" --additional-attr " additionalAttributes} \
        ~{sep=" " inputBams} \
        ~{gtfFile} \
        -c ~{outputTable}
    }

    output {
        File counts = outputTable
    }

    runtime {
        cpu: nprocesses
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBams: {description: "The input BAM files.", category: "required"}
        gtfFile: {description: "A GTF/GFF file containing the features of interest.", category: "required"}
        outputTable: {description: "The path to which the output table should be written.", category: "common"}
        order: {description: "Equivalent to the -r option of htseq-count.", category: "advanced"}
        stranded: {description: "Equivalent to the -s option of htseq-count.", category: "common"}
        additionalAttributes: {description: "Equivalent to the --additional-attr option of htseq-count.", category: "advanced"}
        featureType: {description: "Equivalent to the --type option of htseq-count.", category: "advanced"}
        idattr: {description: "Equivalent to the --idattr option of htseq-count.", category: "advanced"}
        nprocesses: {description: "Number of processes to run htseq with.", category: "advanced"}
        memory: {description: "The amount of memory the job requires in GB.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        counts: {description: "Count table based on input BAM file."}
    }
}
