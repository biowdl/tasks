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

task CreateSamplesheet {
    input {
        File countTable
        String shinyDir = "."

        Int threads = 1
        String memory = "6G"
        Int timeMinutes = 45
        String dockerImage = "tomkuipers1402/shiny-py:v1.0.1"
    }

    command {
        set -e
        mkdir -p ${shinyDir}
        sampleSheet.py \
            -i ${countTable} \
            -o ${shinyDir}
    }

    output {
        File dgeSamples = shinyDir + "/sampleSheet.tsv"
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        countTable: {description: "The created count table from HTseq.", category: "required"}
        shinyDir: {description: "The directory to write the output to.", category: "required"}
        
        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task CreateAnnotation {
    input {
        File referenceFasta
        File referenceGtfFile
        String shinyDir = "."

        Int threads = 1
        String memory = "6G"
        Int timeMinutes = 45
        String dockerImage = "tomkuipers1402/shiny-py:v1.0.1"
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
        File dgeAnnotation = shinyDir + "/annotation.tsv"
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        referenceFasta: {description: "The reference Fasta file.", category: "required"}
        referenceGtfFile: {description: "The reference GTF file.", category: "required"}
        shinyDir: {description: "The directory to write the output to.", category: "required"}
        
        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
