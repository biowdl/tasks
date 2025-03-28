version 1.0

# Copyright (c) 2025 Leiden University Medical Center
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

task Pileup {
    input {
        File bam
        File bamIndex
        String outputBed = "output.bedMethyl"
        String outputBedGraph = "m_CG0_combined.bedgraph"
        File referenceFasta
        File referenceFastaFai

        Int? intervalSize
        File? includeBed

        Boolean cpg = false
        Boolean combineMods = false
        Boolean combineStrands = false
        String? ignore
        String logFilePath = "modkit.log"

        Int threads = 8
        String memory = "4GiB"
        Int timeMinutes = 2880 / threads  # 2 Days / threads
        String dockerImage = "quay.io/biocontainers/ont-modkit:0.4.2--hcdda2d0_0"
    }

    command <<<
        set -e
        mkdir -p $(dirname ~{outputBed})
        mkdir -p $(dirname ~{logFilePath})
        modkit pileup \
        --threads ~{threads} \
        ~{"--interval-size " + intervalSize} \
        ~{"--include-bed " + includeBed} \
        ~{"--ignore " + ignore} \
        --ref ~{referenceFasta} \
        ~{true="--cpg" false="" cpg} \
        ~{true="--combine-mods" false="" combineMods} \
        ~{true="--combine-strands" false="" combineStrands} \
        --log-filepath ~{logFilePath} \
        ~{bam} \
         - | tee ~{outputBed} | awk -v OFS="\t" '{print $1, $2, $3, $11, $10}' > ~{outputBedGraph}
    >>>

    # You can use modkit pileup ${bam_path} - | tee out.bedmethyl | awk -v OFS="\t" '{print $1, $2, $3, $11, $10}' > out.bg to get both outputs at once without running anything twice.
    # https://github.com/nanoporetech/modkit/issues/210#issuecomment-2181706374

    output {
        File out = outputBed  # Normal mode
        File outFiles = outputBedGraph # Bedgraph mode
        File logFile = logFilePath
    }

    runtime {
        docker: dockerImage
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
    }

    parameter_meta {
        # input
        bam: {description: "The input alignment file", category: "required"}
        bamIndex: {description: "The index for the input alignment file", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        outputBed: {description: "The output name where the bedMethyl file should be placed.", category: "common"}
        outputBedgraph: {description: "The output name where the bedgraph file should be placed", category: "common"}

        intervalSize: {description: "Sets the interval size", category: "advanced"}
        includeBed: {description: "Bed file with regions to include", category: "advanced"}
        cpg: {description: "Whether to call only at cpg sites", category: "advanced"}
        combineMods: {description: "Whether to combine modifications in the output", category: "advanced"}
        combineStrands: {description: "Whether to combine strands in the output", category: "advanced"}
        ignore: {description: "Modification type to ignore. For example 'h'.", category: "advanced"}
        logFilePath: {description: "Path where the log file should be written.", category: "advanced"}

        threads: {description: "The number of threads to use for variant calling.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"} 
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    
        # output
        out: {description: "The output bed files. Not available when bedgraph = true."}
        outFiles: {description: "Output files when bedgraph = true."}
        logFile: {description: "The generated log file."}
    }
}
