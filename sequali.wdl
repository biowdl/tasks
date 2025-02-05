version 1.0

# Copyright (c) 2024 Leiden University Medical Center
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

task Sequali {
    input {
        File reads
        File? mate_reads 
        String outDir = "."

        Int threads = 2
        String memory = "4GiB"
        String dockerImage = "quay.io/biocontainers/sequali:0.12.0--py312hf67a6ed_0"
        Int timeMinutes = 59
    } 

    command <<<
        set -e 
        mkdir -p $(dirname outputDir)
        sequali \
        --outdir ~{outDir} \
        --threads ~{threads} \
        ~{reads} \
        ~{mate_reads} 
    >>>
    
    output {
        File html = outDir + "/" + basename(reads) + ".html"
        File json = outDir + "/" + basename(reads) + ".json"
    }

    runtime {
        cpu: threads
        memory: memory
        docker: dockerImage
        time_minutes: timeMinutes
    }

    parameter_meta {
        # inputs
        reads: {description: "A FASTQ or BAM file.", category: "required"}
        mate_reads: {description: "FASTQ mate file"}
        threads: {description: "The number of cores to use.", category: "advanced"}

        outDir: {description: "The path to write the output to.", catgory: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        html: {description: "HTML report file."}
        json: {description: "JSON report file for use with MultiQC."}
    }
}