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

task ParallelPaired {
    input {
        String outputDir
        File referenceFasta
        File referenceFastaFai
        File tumorBam
        File tumorBamIndex
        File normalBam
        File normalBamIndex

        File? classifierSNV
        File? classifierIndel
        File? inclusionRegion
        File? exclusionRegion
        File? mutect2VCF
        File? varscanSNV
        File? varscanIndel
        File? jsmVCF
        File? somaticsniperVCF
        File? vardictVCF
        File? museVCF
        File? lofreqSNV
        File? lofreqIndel
        File? scalpelVCF
        File? strelkaSNV
        File? strelkaIndel

        Int threads = 1
        Int timeMinutes = 60
        String dockerImage = "lethalfang/somaticseq:3.1.0"
    }

    command {
        /opt/somaticseq/somaticseq_parallel.py \
        ~{"--classifier-snv " + classifierSNV} \
        ~{"--classifier-indel " + classifierIndel} \
        --output-directory ~{outputDir} \
        --genome-reference ~{referenceFasta} \
        ~{"--inclusion-region " + inclusionRegion} \
        ~{"--exclusion-region " + exclusionRegion} \
        --threads ~{threads} \
        paired \
        --tumor-bam-file ~{tumorBam} \
        --normal-bam-file ~{normalBam} \
        ~{"--mutect2-vcf " + mutect2VCF} \
        ~{"--varscan-snv " + varscanSNV} \
        ~{"--varscan-indel " + varscanIndel} \
        ~{"--jsm-vcf " + jsmVCF} \
        ~{"--somaticsniper-vcf " + somaticsniperVCF} \
        ~{"--vardict-vcf " + vardictVCF} \
        ~{"--muse-vcf " + museVCF} \
        ~{"--lofreq-snv " + lofreqSNV} \
        ~{"--lofreq-indel " + lofreqIndel} \
        ~{"--scalpel-vcf " + scalpelVCF} \
        ~{"--strelka-snv " + strelkaSNV} \
        ~{"--strelka-indel " + strelkaIndel}
    }

    output {
        File indels = outputDir + if defined(classifierIndel)
            then "/SSeq.Classified.sINDEL.vcf"
            else "/Consensus.sINDEL.vcf"
        File snvs = outputDir + if defined(classifierSNV)
            then "/SSeq.Classified.sSNV.vcf"
            else "/Consensus.sSNV.vcf"
        File ensembleIndels = outputDir + "/Ensemble.sINDEL.tsv"
        File ensembleSNV = outputDir + "/Ensemble.sSNV.tsv"
    }

    runtime {
        cpu: threads
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        outputDir: {description: "The directory to write the output to.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        tumorBam: {description: "The tumor/case sample's BAM file.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor/case sample's BAM file.", category: "required"}
        normalBam: {description: "The normal/control sample's BAM file.", category: "required"}
        normalBamIndex: {description: "The index for the normal/control sample's BAM file.", category: "required"}
        classifierSNV: {description: "A somaticseq SNV classifier.", category: "common"}
        classifierIndel: {description: "A somaticseq Indel classifier.", category: "common"}
        inclusionRegion: {description: "A bed file describing regions to include.", category: "common"}
        exclusionRegion: {description: "A bed file describing regions to exclude.", category: "common"}
        mutect2VCF: {description: "A VCF as produced by mutect2.", category: "advanced"}
        varscanSNV: {description: "An SNV VCF as produced by varscan.", category: "advanced"}
        varscanIndel: {description: "An indel VCF as produced by varscan.", category: "advanced"}
        jsmVCF: {description: "A VCF as produced by jsm.", category: "advanced"}
        somaticsniperVCF: {description: "A VCF as produced by somaticsniper.", category: "advanced"}
        vardictVCF: {description: "A VCF as produced by vardict.", category: "advanced"}
        museVCF: {description: "A VCF as produced by muse.", category: "advanced"}
        lofreqSNV: {description: "An SNV VCF as produced by lofreq.", category: "advanced"}
        lofreqIndel: {description: "An indel VCF as produced by lofreq.", category: "advanced"}
        scalpelVCF: {description: "A VCF as produced by scalpel.", category: "advanced"}
        strelkaSNV: {description: "An SNV VCF as produced by strelka.", category: "advanced"}
        strelkaIndel: {description: "An indel VCF as produced by somaticsniper.", category: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        indels: {description: ""}
        snvs: {description: ""}
        ensembleIndels: {description: ""}
        ensembleSNV: {description: ""}
    }
}

task ParallelPairedTrain {
    input {
        File truthSNV
        File truthIndel
        String outputDir
        File referenceFasta
        File referenceFastaFai
        File tumorBam
        File tumorBamIndex
        File normalBam
        File normalBamIndex

        File? inclusionRegion
        File? exclusionRegion
        File? mutect2VCF
        File? varscanSNV
        File? varscanIndel
        File? jsmVCF
        File? somaticsniperVCF
        File? vardictVCF
        File? museVCF
        File? lofreqSNV
        File? lofreqIndel
        File? scalpelVCF
        File? strelkaSNV
        File? strelkaIndel

        Int threads = 1
        Int timeMinutes = 240
        String dockerImage = "lethalfang/somaticseq:3.1.0"
    }

    command {
        /opt/somaticseq/somaticseq_parallel.py \
        --somaticseq-train \
        --truth-snv ~{truthSNV} \
        --truth-indel ~{truthIndel} \
        --output-directory ~{outputDir} \
        --genome-reference ~{referenceFasta} \
        ~{"--inclusion-region " + inclusionRegion} \
        ~{"--exclusion-region " + exclusionRegion} \
        --threads ~{threads} \
        paired \
        --tumor-bam-file ~{tumorBam} \
        --normal-bam-file ~{normalBam} \
        ~{"--mutect2-vcf " + mutect2VCF} \
        ~{"--varscan-snv " + varscanSNV} \
        ~{"--varscan-indel " + varscanIndel} \
        ~{"--jsm-vcf " + jsmVCF} \
        ~{"--somaticsniper-vcf " + somaticsniperVCF} \
        ~{"--vardict-vcf " + vardictVCF} \
        ~{"--muse-vcf " + museVCF} \
        ~{"--lofreq-snv " + lofreqSNV} \
        ~{"--lofreq-indel " + lofreqIndel} \
        ~{"--scalpel-vcf " + scalpelVCF} \
        ~{"--strelka-snv " + strelkaSNV} \
        ~{"--strelka-indel " + strelkaIndel}
    }

    output {
        File consensusIndels = outputDir + "/Consensus.sINDEL.vcf"
        File consensusSNV = outputDir + "/Consensus.sSNV.vcf"
        File ensembleIndels = outputDir + "/Ensemble.sINDEL.tsv"
        File ensembleSNV = outputDir + "/Ensemble.sSNV.tsv"
        File ensembleIndelsClassifier = outputDir + "/Ensemble.sINDEL.tsv.ntChange.Classifier.RData"
        File ensembleSNVClassifier = outputDir + "/Ensemble.sSNV.tsv.ntChange.Classifier.RData"
    }

    runtime {
        cpu: threads
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        truthSNV: {description: "A VCF of true SNVs.", category: "required"}
        truthIndel: {description: "A VCF of true indels.", category: "required"}
        outputDir: {description: "The directory to write the output to.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        tumorBam: {description: "The tumor/case sample's BAM file.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor/case sample's BAM file.", category: "required"}
        normalBam: {description: "The normal/control sample's BAM file.", category: "required"}
        normalBamIndex: {description: "The index for the normal/control sample's BAM file.", category: "required"}
        inclusionRegion: {description: "A bed file describing regions to include.", category: "common"}
        exclusionRegion: {description: "A bed file describing regions to exclude.", category: "common"}
        mutect2VCF: {description: "A VCF as produced by mutect2.", category: "advanced"}
        varscanSNV: {description: "An SNV VCF as produced by varscan.", category: "advanced"}
        varscanIndel: {description: "An indel VCF as produced by varscan.", category: "advanced"}
        jsmVCF: {description: "A VCF as produced by jsm.", category: "advanced"}
        somaticsniperVCF: {description: "A VCF as produced by somaticsniper.", category: "advanced"}
        vardictVCF: {description: "A VCF as produced by vardict.", category: "advanced"}
        museVCF: {description: "A VCF as produced by muse.", category: "advanced"}
        lofreqSNV: {description: "An SNV VCF as produced by lofreq.", category: "advanced"}
        lofreqIndel: {description: "An indel VCF as produced by lofreq.", category: "advanced"}
        scalpelVCF: {description: "A VCF as produced by scalpel.", category: "advanced"}
        strelkaSNV: {description: "An SNV VCF as produced by strelka.", category: "advanced"}
        strelkaIndel: {description: "An indel VCF as produced by somaticsniper.", category: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        consensusIndels: {description: ""}
        consensusSNV: {description: ""}
        ensembleIndels: {description: ""}
        ensembleSNV: {description: ""}
        ensembleIndelsClassifier: {description: ""}
        ensembleSNVClassifier: {description: ""}
    }
}

task ParallelSingle {
    input {
        File bam
        File bamIndex
        String outputDir
        File referenceFasta
        File referenceFastaFai

        File? classifierSNV
        File? classifierIndel
        File? inclusionRegion
        File? exclusionRegion
        File? mutect2VCF
        File? varscanVCF
        File? vardictVCF
        File? lofreqVCF
        File? scalpelVCF
        File? strelkaVCF

        Int threads = 1
        Int timeMinutes = 60
        String dockerImage = "lethalfang/somaticseq:3.1.0"
    }

    command {
        /opt/somaticseq/somaticseq_parallel.py \
        ~{"--classifier-snv " + classifierSNV} \
        ~{"--classifier-indel " + classifierIndel} \
        --output-directory ~{outputDir} \
        --genome-reference ~{referenceFasta} \
        ~{"--inclusion-region " + inclusionRegion} \
        ~{"--exclusion-region " + exclusionRegion} \
        --threads ~{threads} \
        single \
        --bam-file ~{bam} \
        ~{"--mutect2-vcf " + mutect2VCF} \
        ~{"--varscan-vcf " + varscanVCF} \
        ~{"--vardict-vcf " + vardictVCF} \
        ~{"--lofreq-vcf " + lofreqVCF} \
        ~{"--scalpel-vcf " + scalpelVCF} \
        ~{"--strelka-vcf " + strelkaVCF}
    }

    output {
        File indels = outputDir + if defined(classifierIndel)
            then "/SSeq.Classified.sINDEL.vcf"
            else "/Consensus.sINDEL.vcf"
        File snvs = outputDir + if defined(classifierSNV)
            then "/SSeq.Classified.sSNV.vcf"
            else "/Consensus.sSNV.vcf"
        File ensembleIndels = outputDir + "/Ensemble.sINDEL.tsv"
        File ensembleSNV = outputDir + "/Ensemble.sSNV.tsv"
    }

    runtime {
        cpu: threads
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bam: {description: "The input BAM file.", category: "required"}
        bamIndex: {description: "The index for the input BAM file.", category: "required"}
        outputDir: {description: "The directory to write the output to.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        classifierSNV: {description: "A somaticseq SNV classifier.", category: "common"}
        classifierIndel: {description: "A somaticseq Indel classifier.", category: "common"}
        inclusionRegion: {description: "A bed file describing regions to include.", category: "common"}
        exclusionRegion: {description: "A bed file describing regions to exclude.", category: "common"}
        mutect2VCF: {description: "A VCF as produced by mutect2.", category: "advanced"}
        varscanVCF: {description: "A VCF as produced by varscan.", category: "advanced"}
        vardictVCF: {description: "A VCF as produced by vardict.", category: "advanced"}
        lofreqVCF: {description: "A VCF as produced by lofreq.", category: "advanced"}
        scalpelVCF: {description: "A VCF as produced by scalpel.", category: "advanced"}
        strelkaVCF: {description: "A VCF as produced by strelka.", category: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        indels: {description: ""}
        snvs: {description: ""}
        ensembleIndels: {description: ""}
        ensembleSNV: {description: ""}
    }
}

task ParallelSingleTrain {
    input {
        File bam
        File bamIndex
        File truthSNV
        File truthIndel
        String outputDir
        File referenceFasta
        File referenceFastaFai

        File? inclusionRegion
        File? exclusionRegion
        File? mutect2VCF
        File? varscanVCF
        File? vardictVCF
        File? lofreqVCF
        File? scalpelVCF
        File? strelkaVCF

        Int threads = 1
        Int timeMinutes = 240
        String dockerImage = "lethalfang/somaticseq:3.1.0"
    }

    command {
        /opt/somaticseq/somaticseq_parallel.py \
        --somaticseq-train \
        --truth-snv ~{truthSNV} \
        --truth-indel ~{truthIndel} \
        --output-directory ~{outputDir} \
        --genome-reference ~{referenceFasta} \
        ~{"--inclusion-region " + inclusionRegion} \
        ~{"--exclusion-region " + exclusionRegion} \
        --threads ~{threads} \
        single \
        --bam-file ~{bam} \
        ~{"--mutect2-vcf " + mutect2VCF} \
        ~{"--varscan-vcf " + varscanVCF} \
        ~{"--vardict-vcf " + vardictVCF} \
        ~{"--lofreq-vcf " + lofreqVCF} \
        ~{"--scalpel-vcf " + scalpelVCF} \
        ~{"--strelka-vcf " + strelkaVCF}
    }

    output {
        File consensusIndels = outputDir + "/Consensus.sINDEL.vcf"
        File consensusSNV = outputDir + "/Consensus.sSNV.vcf"
        File ensembleIndels = outputDir + "/Ensemble.sINDEL.tsv"
        File ensembleSNV = outputDir + "/Ensemble.sSNV.tsv"
        File ensembleIndelsClassifier = outputDir + "/Ensemble.sINDEL.tsv.ntChange.Classifier.RData"
        File ensembleSNVClassifier = outputDir + "/Ensemble.sSNV.tsv.ntChange.Classifier.RData"
    }

    runtime {
        cpu: threads
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bam: {description: "The input BAM file.", category: "required"}
        bamIndex: {description: "The index for the input BAM file.", category: "required"}
        truthSNV: {description: "A VCF of true SNVs.", category: "required"}
        truthIndel: {description: "A VCF of true indels.", category: "required"}
        outputDir: {description: "The directory to write the output to.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        inclusionRegion: {description: "A bed file describing regions to include.", category: "common"}
        exclusionRegion: {description: "A bed file describing regions to exclude.", category: "common"}
        mutect2VCF: {description: "A VCF as produced by mutect2.", category: "advanced"}
        varscanVCF: {description: "A VCF as produced by varscan.", category: "advanced"}
        vardictVCF: {description: "A VCF as produced by vardict.", category: "advanced"}
        lofreqVCF: {description: "A VCF as produced by lofreq.", category: "advanced"}
        scalpelVCF: {description: "A VCF as produced by scalpel.", category: "advanced"}
        strelkaVCF: {description: "A VCF as produced by strelka.", category: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        consensusIndels: {description: ""}
        consensusSNV: {description: ""}
        ensembleIndels: {description: ""}
        ensembleSNV: {description: ""}
        ensembleIndelsClassifier: {description: ""}
        ensembleSNVClassifier: {description: ""}
    }
}

task ModifyStrelka {
    input {
        File strelkaVCF
        String outputVCFName = basename(strelkaVCF, ".gz")

        Int timeMinutes = 20
        String dockerImage = "lethalfang/somaticseq:3.1.0"
    }

    command {
        set -e
        /opt/somaticseq/vcfModifier/modify_Strelka.py \
        -infile ~{strelkaVCF} \
        -outfile "modified_strelka.vcf"
        first_FORMAT_line_num=$(grep -n -m 1 '##FORMAT' "modified_strelka.vcf" | cut -d : -f 1)
        sed "$first_FORMAT_line_num"'i##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' "modified_strelka.vcf" > ~{outputVCFName}
    }

    output {
        File outputVcf = outputVCFName
    }

    runtime {
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        strelkaVCF: {description: "A vcf file as produced by strelka.", category: "required"}
        outputVCFName: {description: "The location the output VCF file should be written to.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputVcf: {description: ""}
    }
}
