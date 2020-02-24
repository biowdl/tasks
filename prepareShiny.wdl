version 1.0

# Copyright (c) 2017 Sequencing Analysis Support Core - Leiden University Medical Center
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

task CreateSampleGetCount {
    input {
        File countTable
        String shinyDir = "/shiny"
        String dockerImage = "tomkuipers1402/shiny-py"

        String memory = "6G"
    }

    command {
        set -e
        mkdir -p ${shinyDir}
        sampleSheet.py \
            -i ${countTable} \
            -o ${shinyDir}
        cp ${countTable} ${shinyDir}/allCounts.tsv
    }

    output {
        String shinyOutput = shinyDir
    }

    runtime {
        memory: memory
        docker: dockerImage
    }
}

task CreateAnnotation {
    input {
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        File? referenceGtfFile
        String shinyDir = "/shiny"
        String dockerImage = "tomkuipers1402/shiny-py"

        String memory = "6G"
    }

    command {
        set -e
        mkdir -p ${shinyDir}
        annoGenerator.py \
            -r ${referenceFasta} \
            -g ${referenceGtfFile} \
            -o ${shinyDir}
    }

    output {
        String shinyOutput = shinyDir
    }

    runtime {
        memory: memory
        docker: dockerImage
    }
}
