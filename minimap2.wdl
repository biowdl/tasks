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

task Indexing {
    input {
        Boolean useHomopolymerCompressedKmer = false
        Int kmerSize = 15
        Int minimizerWindowSize = 10
        String outputPrefix
        File referenceFile

        Int? splitIndex

        Int cores = 1
        String memory = "4GiB"
        Int timeMinutes = 10
        String dockerImage = "quay.io/biocontainers/minimap2:2.20--h5bf99c6_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        minimap2 \
        ~{true="-H" false="" useHomopolymerCompressedKmer} \
        -k ~{kmerSize} \
        -w ~{minimizerWindowSize} \
        ~{"-d " + outputPrefix + ".mmi"} \
        -t ~{cores} \
        ~{"-I " + splitIndex} \
        ~{referenceFile}
    }

    output {
        File indexFile = outputPrefix + ".mmi"
    }

    runtime {
        cpu: cores
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        useHomopolymerCompressedKmer: {description: "Use homopolymer-compressed k-mer (preferrable for pacbio).", category: "advanced"}
        kmerSize: {description: "K-mer size (no larger than 28).", category: "advanced"}
        minimizerWindowSize: {description: "Minimizer window size.", category: "advanced"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        referenceFile: {description: "Reference fasta file.", category: "required"}
        splitIndex: {description: "Split index for every ~NUM input bases.", category: "advanced"}
        cores: {description: "The number of cores to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        indexFile: {description: "Indexed reference file."}
    }
}

task Mapping {
    input {
        String presetOption
        Int kmerSize = 15
        Boolean skipSelfAndDualMappings = false
        Boolean outputSam = false
        String outputPrefix
        Boolean addMDTagToSam = false
        Boolean secondaryAlignment = false
        File referenceFile
        File queryFile

        Int? maxIntronLength
        Int? maxFragmentLength
        Int? retainMaxSecondaryAlignments
        Int? matchingScore
        Int? mismatchPenalty
        String? howToFindGTAG

        Int cores = 4
        String memory = "30GiB"
        Int timeMinutes = 1 + ceil(size(queryFile, "G") * 200 / cores)
        String dockerImage = "quay.io/biocontainers/minimap2:2.20--h5bf99c6_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        minimap2 \
        -x ~{presetOption} \
        -k ~{kmerSize} \
        ~{true="-X" false="" skipSelfAndDualMappings} \
        ~{true="-a" false="" outputSam} \
        -o ~{outputPrefix} \
        ~{true="--MD" false="" addMDTagToSam} \
        --secondary=~{true="yes" false="no" secondaryAlignment} \
        -t ~{cores} \
        ~{"-G " + maxIntronLength} \
        ~{"-F " + maxFragmentLength} \
        ~{"-N " + retainMaxSecondaryAlignments} \
        ~{"-A " + matchingScore} \
        ~{"-B " + mismatchPenalty} \
        ~{"-u " + howToFindGTAG} \
        ~{referenceFile} \
        ~{queryFile}
    }

    output {
        File alignmentFile = outputPrefix
    }

    runtime {
        cpu: cores
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        presetOption: {description: "This option applies multiple options at the same time.", category: "common"}
        kmerSize: {description: "K-mer size (no larger than 28).", category: "advanced"}
        skipSelfAndDualMappings: {description: "Skip self and dual mappings (for the all-vs-all mode).", category: "advanced"}
        outputSam: {description: "Output in the sam format.", category: "common"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        addMDTagToSam: {description: "Adds a MD tag to the sam output file.", category: "common"}
        secondaryAlignment: {description: "Whether to output secondary alignments.", category: "advanced"}
        referenceFile: {description: "Reference fasta file.", category: "required"}
        queryFile: {description: "Input fasta file.", category: "required"}
        maxIntronLength: {description: "Max intron length (effective with -xsplice; changing -r).", category: "advanced"}
        maxFragmentLength: {description: "Max fragment length (effective with -xsr or in the fragment mode).", category: "advanced"}
        retainMaxSecondaryAlignments: {description: "Retain at most N secondary alignments.", category: "advanced"}
        matchingScore: {description: "Matching score.", category: "advanced"}
        mismatchPenalty: {description: "Mismatch penalty.", category: "advanced"}
        howToFindGTAG: {description: "How to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG.", category: "common"}
        cores: {description: "The number of cores to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        alignmentFile: {description: "Mapping and alignment between collections of dna sequences file."}
    }
}
