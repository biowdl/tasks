version 1.0

# Copyright (c) 2025 Leiden University Medical Center
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


task ValidateFastq {
    input {
        File inputRead1
        File? inputRead2

        String memory = "1GiB"
        Int timeMinutes = 5 + ceil(size(inputRead1, "GiB"))
        String dockerImage = "quay.io/biocontainers/biopet-validatefastq:0.1.1--hdfd78af_3"
    }

    command {
        set -e
        java -jar /usr/local/share/biopet-validatefastq-0.1.1-3/validatefastq-assembly-0.1.1.jar \
        --fastq1 ~{inputRead1} \
        ~{"--fastq2 " + inputRead2}
    }

    output {
    }

    runtime {
        cpu: 1
        memory: memory
        docker: dockerImage
        time_minutes: timeMinutes
    }

    parameter_meta {
        # inputs
        inputRead1: {description: "The location of the first FASTQ file (first reads for pairs, in case of paired-end sequencing).", category: "required"}
        inputRead2: {description: "The location of the paired end reads.", category: "common"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}
