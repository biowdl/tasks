version 1.0

# Copyright (c) 2019 Leiden University Medical Center
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

task Fastqsplitter {
    input {
        File inputFastq
        Array[String]+ outputPaths

        Int? compressionLevel
        Int? threadsPerFile

        # fastqplitter utilizes one thread per input file and one or
        # more threads per output file + one thread for the application.
        # Since a compression level of 1 is used, each output file
        # uses approx 0.5 cores.
        Int cores = 1 + ceil(0.5 * length(outputPaths))
        String dockerImage = "quay.io/biocontainers/fastqsplitter:1.1.0--py37h516909a_1"
    }

    # Busybox mkdir does not accept multiple paths.
    command <<<
        set -e
        for FILE in ~{sep=' ' outputPaths}
        do
            mkdir -p "$(dirname ${FILE})"
        done
        fastqsplitter \
        ~{"-c " + compressionLevel} \
        ~{"-t " + threadsPerFile} \
        -i ~{inputFastq} \
        -o ~{sep=' -o ' outputPaths}
    >>>

    output {
        Array[File] chunks = outputPaths
    }

    # Using very safe margins here. 10MB/300MB per outputfile is used for
    # single-threaded/multi-threaded compression.
    Float memoryPerFile = if select_first([threadsPerFile, 1]) > 1 then 0.40 else 0.02
    Int fastqsplitterMemory = ceil(0.100 + memoryPerFile * length(outputPaths))
    # Make sure a minimum of 2 GB is present to pull the singularity image.
    Int memory = if fastqsplitterMemory <= 2 then 2 else fastqsplitterMemory

    runtime {
        cpu: cores
        memory: "~{memory}GiB"
        docker: dockerImage
    }
}
