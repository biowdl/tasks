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

task ChunkedScatter {
    input {
        File inputFile
        String prefix = "./scatter"
        Int? chunkSize
        Int? overlap
        Int? minimumBasesPerFile

        Int timeMinutes = 2
        String dockerImage = "quay.io/biocontainers/chunked-scatter:0.1.0--py_0"
    }

    command {
        set -e
        mkdir -p ~{prefix}
        chunked-scatter \
        -p ~{prefix} \
        -i ~{inputFile} \
        ~{"-c " + chunkSize} \
        ~{"-o " + overlap} \
        ~{"-m " + minimumBasesPerFile}
    }

    output {
        Array[File] scatters = glob(prefix + "*.bed")
    }

    runtime {
        memory: "4G"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        inputFile: {description: "Either a bed file describing regiosn of intrest or a sequence dictionary.", category: "required"}
        prefix: {description: "The prefix for the output files.", category: "advanced"}
        chunkSize: {description: "Equivalent to chunked-scatter's `-c` option.", category: "advanced"}
        overlap: {description: "Equivalent to chunked-scatter's `-o` option.", category: "advanced"}
        minimumBasesPerFile: {description: "Equivalent to chunked-scatter's `-m` option.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}