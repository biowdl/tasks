version 1.0

# Copyright (c) 2020 Leiden University Medical Center
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

task Format {
    input {
        Array[File]+ inputFiles
        String format = "fasta"
        String outputPath = "seq_data.sdf"

        String rtgMem = "8G"
        String memory = "9G"
        Int timeMinutes = 1 + ceil(size(inputFiles) * 2)
        String dockerImage = "quay.io/biocontainers/rtg-tools:3.10.1--0"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPath})
        rtg RTG_MEM=~{rtgMem} format -f ~{format} \
        -o ~{outputPath} \
        ~{sep=' ' inputFiles}
    }

    output {
        File sdf = outputPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputFiles: {description: "Input sequence files. May be specified 1 or more times.", category: "required"}
        format: {description: "Format of input. Allowed values are [fasta, fastq, fastq-interleaved, sam-se, sam-pe].", category: "advanced"}
        outputPath: {description: "Where the output should be placed.", category: "advanced"}
        rtgMem: {description: "The amount of memory rtg will allocate to the JVM.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        sdf: {description: "RTGSequence Data File (SDF) format version of the input file(s)."}
    }
}

task VcfEval {
    input {
        File baseline
        File baselineIndex
        File calls
        File callsIndex
        Boolean squashPloidy = false
        String outputMode = "split"
        String outputDir = "output/"
        File template
        Boolean allRecords = false
        Boolean decompose = false
        Boolean refOverlap = false

        File? evaluationRegions
        File? bedRegions
        String? sample

        String rtgMem = "8G"
        Int threads = 1  # Tool default is number of cores in the system 😱.
        String memory = "9G"
        Int timeMinutes = 1 + ceil(size([baseline, calls], "G") * 5)
        String dockerImage = "quay.io/biocontainers/rtg-tools:3.10.1--0"
    }

    command <<<
        set -e
        mkdir -p "$(dirname ~{outputDir})"
        rtg RTG_MEM=~{rtgMem} vcfeval \
        --baseline ~{baseline} \
        --calls ~{calls} \
        ~{"--evaluation-regions " + evaluationRegions} \
        ~{"--bed-regions " + bedRegions} \
        --output ~{outputDir} \
        --template ~{template} \
        ~{true="--all-records" false="" allRecords} \
        ~{true="--decompose" false="" decompose} \
        ~{true="--ref-overlap" false="" refOverlap} \
        ~{"--sample " + sample } \
        ~{true="--squash-ploidy" false="" squashPloidy} \
        ~{"--output-mode " + outputMode} \
        --threads ~{threads}
    >>>

    output {
        File falseNegativesVcf = outputDir + "/fn.vcf.gz"
        File falseNegativesVcfIndex = outputDir + "/fn.vcf.gz.tbi"
        File falsePositivesVcf = outputDir + "/fp.vcf.gz"
        File falsePositivesVcfIndex = outputDir + "/fp.vcf.gz.tbi"
        File summary = outputDir + "/summary.txt"
        File truePositivesBaselineVcf = outputDir + "/tp-baseline.vcf.gz"
        File truePositivesBaselineVcfIndex = outputDir + "/tp-baseline.vcf.gz.tbi"
        File truePositivesVcf = outputDir + "/tp.vcf.gz"
        File truePositivesVcfIndex = outputDir + "/tp.vcf.gz.tbi"
        File nonSnpRoc = outputDir + "/non_snp_roc.tsv.gz"
        File phasing = outputDir + "/phasing.txt"
        File weightedRoc = outputDir + "/weighted_roc.tsv.gz"
        Array[File] allStats = [falseNegativesVcf,
                                falseNegativesVcfIndex,
                                falsePositivesVcf,
                                falsePositivesVcfIndex,
                                truePositivesBaselineVcf,
                                truePositivesBaselineVcfIndex,
                                truePositivesVcf,
                                truePositivesVcfIndex,
                                summary,
                                nonSnpRoc,
                                phasing,
                                weightedRoc]
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        baseline: {description: "VCF file containing baseline variants.", category: "required"}
        baselineIndex: {description: "The baseline's VCF index.", category: "required"}
        calls: {description: "VCF file containing called variants.", category: "required"}
        callsIndex: {description: "The call's VCF index.", category: "required"}
        squashPloidy: {description: "treat heterozygous genotypes as homozygous ALT in both baseline and calls, to allow matches that ignore zygosity differences.", category: "common"}
        outputMode: {description: "output reporting mode. Allowed values are [split, annotate, combine, ga4gh, roc-only] (Default is split).", category: "advanced"}
        outputDir: {description: "Directory for output.", category: "advanced"}
        template: {description: "SDF of the reference genome the variants are called against.", category: "required"}
        allRecords: {description: "use all records regardless of FILTER status (Default is to only process records where FILTER is \".\" or \"PASS\").", category: "common"}
        decompose: {description: "decompose complex variants into smaller constituents to allow partial credit.", category: "common"}
        refOverlap: {description: "allow alleles to overlap where bases of either allele are same-as-ref (Default is to only allow VCF anchor base overlap).", category: "common"}
        sample: {description: "the name of the sample to select. Use <baseline_sample>,<calls_sample> to select different sample names for baseline and calls. (Required when using multi-sample VCF files).", category: "common"}
        bedRegions: {description: "if set, only read VCF records that overlap the ranges contained in the specified BED file.", category: "advanced"}
        evaluationRegions: {description: "if set, evaluate within regions contained in the supplied BED file, allowing transborder matches. To be used for truth-set high-confidence regions or other regions of interest where region boundary effects should be minimized.", category: "advanced"}
        rtgMem: {description: "The amount of memory rtg will allocate to the JVM.", category: "advanced"}
        threads: {description: "Number of threads. Default is 1.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        falseNegativesVcf: {description: "Variants from thebaselineVCF which were not correctly called."}
        falseNegativesVcfIndex: {description: "Index of the output VCF file `falseNegativesVcf`."}
        falsePositivesVcf: {description: "Variants from thecallsVCF which do not agree with baseline variants."}
        falsePositivesVcfIndex: {description: "Index of the output VCF file `falsePositivesVcf`."}
        summary: {description: "Summary statistic file."}
        truePositivesBaselineVcf: {description: "Variants from thebaselineVCF which agree with variants in thecalls VCF."}
        truePositivesBaselineVcfIndex: {description: "Index of the output VCF file `truePositivesBaselineVcf`."}
        truePositivesVcf: {description: "Variants from thecallsVCF which agree with variants in the baseline VCF."}
        truePositivesVcfIndex: {description: "Index of the output VCF file `truePositivesVcf`."}
        nonSnpRoc: {description: "ROC data derived from those variants which were not represented asSNPs."}
        phasing: {description: "Phasing file."}
        weightedRoc: {description: "ROC data derived from all analyzed call variants, regardless of their representation."}
        allStats: {description: "All output files combined in a array."}
    }
}
