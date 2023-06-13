version 1.0

# MIT License
#
# Copyright (c) 2023 Leiden University Medical Center
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

task FastqFilter {
    input {
        Array[File]+ fastq
        Array[String]+ outputPaths
        Int? minLength
        Int? maxLength

        String memory = "1GiB"
        Int timeMinutes = 1 + ceil(size(fastq, "G"))
        String dockerImage = "quay.io/biocontainers/fastq-filter:0.3.0--py39hf95cd2a_1"
    }

    command {
        set -e
        mkdir -p $(dirname ~{sep=" " outputPaths})
        fastq-filter \
        -o ~{sep=" -o " outputPaths} \
        ~{"-l " + minLength} \
        ~{"-L " + maxLength} \
        ~{sep=" " fastq}
    }

    output {
        Array[File] filtered = outputPaths
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        fastq: {description: "A list of fastq files to filter.", category: "required"}
        outputPaths: {description: "A list containing the output paths for each input fastq file.", category: "required"}
        minLength: {description: "Equivalent to fastq-filter's `--min-length` option.", category: "common"}
        maxLength: {description: "Equivalent to fastq-filter's `--max-length` option.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}