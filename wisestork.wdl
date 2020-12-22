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

task Count {
    input {
        File inputBam
        File inputBamIndex
        File reference
        File referenceIndex
        String outputBed = "output.bed"

        Int? binSize
        File? binFile

        String dockerImage = "quay.io/biocontainers/wisestork:0.1.2--pyh24bf2e0_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputBed})"
        wisestork count \
        ~{"--binsize " + binSize} \
       --reference ~{reference} \
       ~{"--bin-file " + binFile} \
       --output ~{outputBed} \
       --input ~{inputBam}
    }

    output {
        File bedFile = outputBed
    }

    runtime {
        docker: dockerImage
    }
}

task GcCorrect {
    input {
        File reference
        File referenceIndex
        File inputBed
        String outputBed = "output.bed"

        Int? binSize
        File? binFile
        Float? fracN
        Int? iter
        Float? fracLowess

        String dockerImage = "quay.io/biocontainers/wisestork:0.1.2--pyh24bf2e0_0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputBed})
        wisestork gc-correct \
        ~{"--binsize " + binSize} \
        --reference ~{reference} \
        ~{"--bin-file " + binFile} \
        --output ~{outputBed} \
        --input ~{inputBed} \
        ~{"--frac-n " + fracN} \
        ~{"--iter " + iter} \
        ~{"--frac-lowess " + fracLowess}
    }

    output {
        File bedFile = outputBed
    }

    runtime {
        docker: dockerImage
    }
}

task Newref {
    input {
        File reference
        File referenceIndex
        Array[File]+ inputBeds
        String outputBed = "output.bed"

        Int? binSize
        File? binFile
        Int? nBins

        Int memory = 2 + ceil(length(inputBeds) * 0.15)
        String dockerImage = "quay.io/biocontainers/wisestork:0.1.2--pyh24bf2e0_0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputBed})
        wisestork newref \
        ~{"--binsize " + binSize} \
        --reference ~{reference} \
        ~{"--bin-file " + binFile} \
        --output ~{outputBed} \
        -I ~{sep=" -I " inputBeds} \
        ~{"--n-bins " + nBins}
    }

    output {
        File bedFile = outputBed
    }

    runtime {
        memory: "~{memory}G"
        docker: dockerImage
    }
}

task Zscore {
    input {
        File reference
        File referenceIndex
        File inputBed
        File inputBedIndex
        File dictionaryFile
        File dictionaryFileIndex
        String outputBed = "output.bed"

        Int? binSize
        File? binFile

        String dockerImage = "quay.io/biocontainers/wisestork:0.1.2--pyh24bf2e0_0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputBed})
        wisestork zscore \
        ~{"--binsize " + binSize} \
        --reference ~{reference} \
        ~{"--bin-file " + binFile} \
        --output ~{outputBed} \
        --input ~{inputBed} \
        ~{"--dictionary-file " + dictionaryFile}
    }

    output {
        File bedFile = outputBed
    }

    runtime {
        docker: dockerImage
    }
}
