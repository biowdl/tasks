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
        String outputBedGraph = "combined.bedgraph"
        File referenceFasta
        File referenceFastaFai

        Int? intervalSize
        File? includeBed
        String? filterThreshold
        String? filterPercentile

        Boolean cpg = false
        Boolean combineMods = false
        Boolean combineStrands = false
        String? ignore
        String logFilePath = "modkit.log"

        Int threads = 8
        String memory = "4GiB"
        Int timeMinutes = 2880 / threads  # 2 Days / threads
        String dockerImage = "quay.io/biocontainers/ont-modkit:0.4.3--hcdda2d0_0"
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
        ~{"--filter-percentile " + filterPercentile} \
        ~{"--filter-threshold " + filterThreshold} \
        --log-filepath ~{logFilePath} \
        ~{bam} \
         - | tee ~{outputBed} | awk -v OFS="\t" '{print $1, $2, $3, $11, $10 >> "~{outputBedGraph}_"$4"_"$6".bedGraph"}'
        # Separately generate the combined file as well, so users can have a choice.
        cat ~{outputBed} | awk -v OFS="\t" '{print $1, $2, $3, $11, $10}' > ~{outputBedGraph}
    >>>

    # You can use modkit pileup ${bam_path} - | tee out.bedmethyl | awk -v OFS="\t" '{print $1, $2, $3, $11, $10}' > out.bg to get both outputs at once without running anything twice.
    # https://github.com/nanoporetech/modkit/issues/210#issuecomment-2181706374

    output {
        File out = outputBed  # Normal mode
        File outGraph = outputBedGraph  # Normal mode
        Array[File] outFiles = glob(outputBedGraph + "*.bedGraph")  # Bedgraph mode
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
        outputBedGraph: {description: "The output name where the bedgraph file should be placed", category: "common"}

        intervalSize: {description: "Sets the interval size", category: "advanced"}
        includeBed: {description: "Bed file with regions to include", category: "advanced"}
        cpg: {description: "Whether to call only at cpg sites", category: "advanced"}
        combineMods: {description: "Whether to combine modifications in the output", category: "advanced"}
        combineStrands: {description: "Whether to combine strands in the output", category: "advanced"}
        ignore: {description: "Modification type to ignore. For example 'h'.", category: "advanced"}
        logFilePath: {description: "Path where the log file should be written.", category: "advanced"}
        filterThreshold: {description: "Global filter threshold can be specified with by a decimal number (e.g. 0.75). Otherwise the automatic filter percentile will be used.", category: "advanced"}
        filterPercentile: {description: "This defaults to 0.1, to remove the lowest 10% confidence modification calls, but can be manually adjusted", category: "advanced"}

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

task Summary {
    input {
        File bam
        File bamIndex

        String summary = "modkit.summary.txt"

        Boolean sample = true
        Int? numReads # = 10042
        Float? samplingFrac # = 0.1
        Int? seed

        Int threads = 4
        String memory = ceil(size(bam, "GiB") * 0.1) + 5 # Based on a linear model with some fudge (memory = 0.07540 * file_size - 0.6).
        Int timeMinutes = 60 # originally this was set at "2 Days / threads" but with 4 threads and that much ram, it's pretty fast.
        String dockerImage = "quay.io/biocontainers/ont-modkit:0.4.3--hcdda2d0_0"
    }

    command <<<
        set -e
        mkdir -p $(dirname ~{summary})

        modkit summary \
        --threads ~{threads} \
        ~{true="" false="--no-sampling" sample} \
        ~{"--num-reads " + numReads} \
        ~{"--sampling-frac " + samplingFrac} \
        ~{"--seed " + seed} \
        ~{bam} > ~{summary}
    >>>

    output {
        File summaryReport = summary # Normal mode
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

        sample: {description: "Allows you to disable sampling and report stats for the whole file.", category: "advanced"}
        numReads: {description: "By default a fixed amount of reads are read, you can set this to change the number of reads to sample.", category: "advanced"}
        samplingFrac: {description: "Use a fixed percentage of reads, rather than a fixed number of reads, for sampling.", category: "advanced"}
        seed: {description: "A seed can be provided for reproducibility in the sampling fraction case.", category: "advanced"}

        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # output
        summaryReport: {description: "The output modkit summary."}
    }
}

task SampleProbs {
    input {
        File bam
        File bamIndex

        String summary = "modkit-sample-probs"

        Boolean sample = true
        Int? numReads # = 10042
        Float? samplingFrac # = 0.1
        Int? seed

        Int threads = 4
        String memory = "32G"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/ont-modkit:0.4.3--hcdda2d0_0"
    }

    command <<<
        set -e
        mkdir -p ~{summary}

        modkit sample-probs \
            --threads ~{threads} \
            --out-dir ~{summary} \
            ~{true="" false="--no-sampling" sample} \
            ~{"--num-reads " + numReads} \
            ~{"--sampling-frac " + samplingFrac} \
            ~{"--seed " + seed} \
            --hist \
            ~{bam}
    >>>

    output {
        File reportCounts = "~{summary}/counts.html"
        File reportProportion = "~{summary}/proportion.html"
        File reportProbabilitiesTsv = "~{summary}/probabilities.tsv"
        File reportThresholdsTsv = "~{summary}/thresholds.tsv"
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
        summary: {description: "A folder for the outputs", category: "required"}

        sample: {description: "Allows you to disable sampling and report stats for the whole file.", category: "advanced"}
        numReads: {description: "By default a fixed amount of reads are read, you can set this to change the number of reads to sample.", category: "advanced"}
        samplingFrac: {description: "Use a fixed percentage of reads, rather than a fixed number of reads, for sampling.", category: "advanced"}
        seed: {description: "A seed can be provided for reproducibility in the sampling fraction case.", category: "advanced"}

        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # output
        reportCounts: {description: "The output html report of counts"}
        reportProportion: {description: "The output html report of proportions"}
        reportProbabilitiesTsv: {description: "The output TSV of Probabilities"}
        reportThresholdsTsv: {description: "The output TSV of thresholds"}
    }
}
