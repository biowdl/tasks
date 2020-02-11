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

task GffRead {
    input {
        File inputGff
        File genomicSequence
        File? genomicIndex  # Optional. GFFRead can create this by itself.
        String? exonsFastaPath
        String? CDSFastaPath
        String? proteinFastaPath
        String? filteredGffPath
        Boolean outputGtfFormat = false
        String dockerImage = "quay.io/biocontainers/gffread:0.9.12--0"
    }

    # The mkdirs below are hackish. It should be
    # ~{"mkir -p $(dirname " + somePath + ")"}
    # but this goes wrong. Cromwell will always use ')' even if somepath is not defined.
    # Which leads to crashing.
    command {
        set -e
        ~{"mkdir -p $(dirname " + CDSFastaPath}~{true=")" false="" defined(CDSFastaPath)}
        ~{"mkdir -p $(dirname " + exonsFastaPath}~{true=")" false="" defined(exonsFastaPath)}
        ~{"mkdir -p $(dirname " + proteinFastaPath}~{true=")" false="" defined(proteinFastaPath)}
        ~{"mkdir -p $(dirname " + filteredGffPath}~{true=")" false="" defined(filteredGffPath)}
        gffread \
        ~{inputGff} \
        -g ~{genomicSequence} \
        ~{"-w " + exonsFastaPath} \
        ~{"-x " + CDSFastaPath} \
        ~{"-y " + proteinFastaPath} \
        ~{"-o " + filteredGffPath} \
        ~{true="-T " false="" outputGtfFormat}
    }

    output {
        File? exonsFasta = exonsFastaPath
        File? CDSFasta = CDSFastaPath
        File? proteinFasta = proteinFastaPath
        File? filteredGff = filteredGffPath
    }

    runtime {
        docker: dockerImage
    }

    parameter_meta {
        inputGff: {description: "The input GFF file.", category: "required"}
        genomicSequence: {description: "The genome.", category: "required"}
        genomicIndex: {description: "The genome's index.", category: "advanced"}
        exonsFastaPath: {description: "The location the exons fasta should be written to.", category: "advanced"}
        CDSFastaPath: {description: "The location the CDS fasta should be written to.", category: "advanced"}
        proteinFastaPath: {description: "The location the protein fasta should be written to.", category: "advanced"}
        filteredGffPath: {description: "The location the filtered GFF should be written to.", category: "advanced"}
        outputGtfFormat: {description: "Equivalent to gffread's `-T` flag.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}