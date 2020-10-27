version 1.0

# Copyright (c) 2020 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

task Mapping {
    input {
        String presetOption
        Boolean sort=true
        String sample
        File referenceMMI
        File queryFile

        Int cores = 4
        String memory = "30G"
        Int timeMinutes = 1 + ceil(size(queryFile, "G") * 2000 / cores)
        String dockerImage = "quay.io/biocontainers/pbmm2:1.3.0--h56fc30b_1"
    }

    command {
        pbmm2 align \
        --preset ~{presetOption} \
        ~{true="--sort" false="" sort} \
        -j ~{cores} \
        ~{referenceMMI} \
        ~{queryFile} \
        --sample ~{sample} \
        ~{sample}.align.bam
    }

    output {
        File outputAlignmentFile = sample + ".align.bam"
        File outputIndexFile = sample + ".align.bam.bai"
    }

    runtime {
        cpu: cores
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        presetOption: {description: "This option applies multiple options at the same time.", category: "required"}
        sort: {description: "Sort the output bam file.", category: "advanced"}
        sample: {description: "Name of the sample"}
        referenceMMI: {description: "MMI file for the reference.", category: "required"}
        queryFile: {description: "BAM file with reads to align against the reference.", category: "required"}
        cores: {description: "The number of cores to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # output
        outputAlignmentFile: {description: "Mapped bam file."}
        outputIndexFile: {description: "Bam index file."}
    }
}
