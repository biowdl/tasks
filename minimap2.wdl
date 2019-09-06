version 1.0

# Copyright (c) 2019 Sequencing Analysis Support Core - Leiden University Medical Center
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

task Indexing {
    input {
        File referenceFile
        String outputPrefix
        Boolean useHomopolymerCompressedKmer = false
        Int kmerSize = 15
        Int minimizerWindowSize = 10

        Int? splitIndex

        Int cores = 1
        Int memory = 4
        String dockerImage = "quay.io/biocontainers/minimap2:2.17--h84994c4_0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        minimap2 \
        ~{true="-H" false="" useHomopolymerCompressedKmer} \
        ~{"-k " + kmerSize} \
        ~{"-w " + minimizerWindowSize} \
        ~{"-I " + splitIndex} \
        ~{"-d " + outputPrefix + ".mmi"} \
        ~{"-t " + cores} \
        ~{referenceFile}
    }

    output {
        File outputIndexFile = outputPrefix + ".mmi"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        referenceFile: "Reference fasta file."
        outputPrefix: "Output directory path + output file prefix."
        useHomopolymerCompressedKmer: "Use homopolymer-compressed k-mer (preferrable for PacBio)."
        kmerSize: "K-mer size (no larger than 28)."
        minimizerWindowSize: "Minimizer window size."
        splitIndex: "Split index for every ~NUM input bases."

        outputIndexFile: "Indexed reference file."
    }
}

task Mapping {
    input {
        File queryFile
        File referenceFile
        String outputPrefix
        String presetOption
        Boolean outputSAM = false

        Int? maxFragmentLength
        Int? maxIntronLength
        Boolean? skipSelfAndDualMappings
        Int? retainMaxSecondaryAlignments
        Int? matchingScore
        Int? mismatchPenalty
        String? howToFindGTAG
        Boolean? secondaryAlignment

        Int cores = 4
        Int memory = 7
        String dockerImage = "quay.io/biocontainers/minimap2:2.17--h84994c4_0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        minimap2 \
        ~{"-x " + presetOption} \
        ~{true="-a" false="" outputSAM} \
        ~{"-G " + maxIntronLength} \
        ~{"-F " + maxFragmentLength} \
        ~{true="-X" false="" skipSelfAndDualMappings} \
        ~{"-N " + retainMaxSecondaryAlignments} \
        ~{"-A " + matchingScore} \
        ~{"-B " + mismatchPenalty} \
        ~{"-u " + howToFindGTAG} \
        --secondary=~{true="yes" false="no" secondaryAlignment} \
        ~{"-o " + outputPrefix} \
        ~{"-t " + cores} \
        ~{referenceFile} \
        ~{queryFile}
    }

    output {
        File outputAlignmentFile = outputPrefix
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        queryFile: "Input fasta file."
        referenceFile: "Reference fasta file."
        outputPrefix: "Output directory path + output file prefix."
        presetOption: "This option applies multiple options at the same time."
        outputSAM: "Output in the SAM format."
        maxFragmentLength: "Max fragment length (effective with -xsr or in the fragment mode)."
        maxIntronLength: "Max intron length (effective with -xsplice; changing -r)."
        skipSelfAndDualMappings: "Skip self and dual mappings (for the all-vs-all mode)."
        retainMaxSecondaryAlignments: "Retain at most INT secondary alignments."
        matchingScore: "Matching score."
        mismatchPenalty: "Mismatch penalty."
        howToFindGTAG: "How to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG."
        secondaryAlignment: "Whether to output secondary alignments."

        outputAlignmentFile: "Mapping and alignment between collections of DNA sequences file."
    }
}
