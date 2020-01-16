version 1.0

# Copyright (c) 2018 Leiden University Medical Center
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

task InputConverter {
    input {
        File samplesheet
        String outputFile = "samplesheet.json"
        # File checking only works when:
        # 1. Paths are absolute
        # 2. When containers have the directory with the files mounted.
        # Therefore this functionality does not work well with cromwell.
        Boolean skipFileCheck=true
        Boolean checkFileMd5sums=false
        Boolean old=false
        String dockerImage = "quay.io/biocontainers/biowdl-input-converter:0.2.1--py_0"
    }

    command <<<
        set -e
        mkdir -p "$(dirname ~{outputFile})"
        biowdl-input-converter \
        -o ~{outputFile} \
        ~{true="--skip-file-check" false="" skipFileCheck} \
        ~{true="--check-file-md5sums" false="" checkFileMd5sums} \
        ~{true="--old" false="" old} \
        ~{samplesheet}
    >>>

    output {
        File json = outputFile
    }

    runtime {
        docker: dockerImage
    }

    parameter_meta {
        samplesheet: {description: "The samplesheet to be processed.", category: "required"}
        outputFile: {description: "The location the JSON representation of the samplesheet should be written to.",
                     category: "advanced"}
        skipFileCheck: {description: "Whether or not the existance of the files mentioned in the samplesheet should be checked.",
                        category: "advanced"}
        checkFileMd5sums: {description: "Whether or not the MD5 sums of the files mentioned in the samplesheet should be checked.",
                           category: "advanced"}
        old: {description: "Whether or not the old samplesheet format should be used.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
