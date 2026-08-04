version 1.0

# Copyright (c) 2026 Leiden University Medical Center
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

task Index {
    input {
        File fasta
        String outputPath = basename(fasta, ".fasta") + ".gem_index"

        Int threads = 10
        String memory = "32GiB"
        String dockerImage = "quay.io/biocontainers/gem2:20200110--h9ee0642_1"
        Int timeMinutes = 1440
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPath})
        gem-mappability \
        -T ~{threads} \
        -i ~{fasta} \
        -o ~{outputPath}
    }

    output {
        File gemIndex = outputPath
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        fasta: {description: "The fasta file to index.", category: "required"}
        outputPath: {description: "The path to write the gem index file to.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Mappability {
    input {
        File gemIndex
        String outputPath = basename(gemIndex, ".gem_index") + ".gem_mappability_~{readLength}"
        Int readLength = 150

        Int threads = 10
        String memory = "32GiB"
        String dockerImage = "quay.io/biocontainers/gem2:20200110--h9ee0642_1"
        Int timeMinutes = 1440
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPath})
        gem-mappability \
        -t ~{threads} \
        -i ~{gemIndex} \
        -o ~{outputPath} \
        -l ~{readLength}
    }

    output {
        File mappability = outputPath
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        gemIndex: {description: "The gem index file to calculate mappability for.", category: "required"}
        outputPath: {description: "The path to write the mappability file to.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
