version 1.0

# Copyright (c) 2020 Leiden University Medical Center
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

task mergePacBio {
    input {
        Array[File]+ reports
        String mergedReport

        String memory = "4G"
        String dockerImage = "lumc/pacbio-merge:0.2"
    }

    command {
        set -e
        mkdir -p $(dirname ~{mergedReport})
        pacbio_merge \
        --reports ~{sep=" " reports} \
        --json-output ~{mergedReport}
    }

    runtime {
        memory: memory
        docker: dockerImage
    }

    output {
        File MergedReport = mergedReport
    }

    parameter_meta {
        # inputs
        reports: {description: "The PacBio report files to merge.", category: "required"}
        mergedReport: {description: "The location the merged PacBio report file should be written to.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}

task ccsChunks {
    input {
        Int chunkCount

        String memory = "4G"
        String dockerImage = "python:3.7-slim"
    }

    command {
        set -e
        python <<CODE
        for i in range(1, ~{chunkCount} + 1):
            print(i, ~{chunkCount}, sep="/")
        CODE
    }

    runtime {
        memory: memory
        docker: dockerImage
    }

    output {
        Array[String] chunks = read_lines(stdout())
    }

    parameter_meta {
        # inputs
        chunkCount: {description: "The number of chunks to create.", category: "required"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}
