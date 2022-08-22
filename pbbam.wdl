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

task Index {
    input {
        File bamFile

        String? outputBamPath

        String memory = "2GiB"
        Int timeMinutes = 1 + ceil(size(bamFile, "G") * 4)
        String dockerImage = "quay.io/biocontainers/pbbam:1.6.0--h5b7e6e0_0"
    }

    # Select_first is needed, otherwise womtool validate fails.
    String outputPath = select_first([outputBamPath, basename(bamFile)])
    String bamIndexPath = outputPath + ".pbi"

    command {
        set -e
        # Make sure outputBamPath does not exist.
        if [ ! -f ~{outputPath} ]
        then
            mkdir -p "$(dirname ~{outputPath})"
            ln ~{bamFile} ~{outputPath}
        fi
        pbindex ~{outputPath} ~{bamIndexPath}
    }

    output {
        File indexedBam = outputPath
        File index =  bamIndexPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bamFile: {description: "The BAM file for which an index should be made.", category: "required"}
        outputBamPath: {description: "The location where the BAM file should be written to. The index will appear alongside this link to the BAM file.", category: "common"}
        memory: {description: "The amount of memory needed for the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        indexedBam: {description: "The original input BAM file."}
        index: {description: "The index of the input BAM file."}
    }
}
