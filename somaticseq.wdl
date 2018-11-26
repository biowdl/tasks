version 1.0

import "common.wdl" as common

task ParallelPaired {
    input {
        String installDir = "/opt/somaticseq" #the location in the docker image

        File? classifierSNV
        File? classifierIndel
        String outputDir
        Reference reference
        File? inclusionRegion
        File? exclusionRegion
        IndexedBamFile tumorBam
        IndexedBamFile normalBam
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
    }

    command {
        ~{installDir}/somaticseq_parallel.py \
        ~{"--classifier-snv " + classifierSNV} \
        ~{"--classifier-indel " + classifierIndel} \
        --output-directory ~{outputDir} \
        --genome-reference ~{reference.fasta} \
        ~{"--inclusion-region " + inclusionRegion} \
        ~{"--exclusion-region " + exclusionRegion} \
        --threads ~{threads} \
        paired \
        --tumor-bam-file ~{tumorBam.file} \
        --normal-bam-file ~{normalBam.file} \
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
        docker: "lethalfang/somaticseq:3.1.0"
    }
}

task ParallelPairedTrain {
    input {
        String installDir = "/opt/somaticseq" #the location in the docker image

        File truthSNV
        File truthIndel
        String outputDir
        Reference reference
        File? inclusionRegion
        File? exclusionRegion
        IndexedBamFile tumorBam
        IndexedBamFile normalBam
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
    }

    command {
        ~{installDir}/somaticseq_parallel.py \
        --somaticseq-train \
        --truth-snv ~{truthSNV} \
        --truth-indel ~{truthIndel} \
        --output-directory ~{outputDir} \
        --genome-reference ~{reference.fasta} \
        ~{"--inclusion-region " + inclusionRegion} \
        ~{"--exclusion-region " + exclusionRegion} \
        --threads ~{threads} \
        paired \
        --tumor-bam-file ~{tumorBam.file} \
        --normal-bam-file ~{normalBam.file} \
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
        docker: "lethalfang/somaticseq:3.1.0"
    }
}

task ParallelSingle {
    input {
        String installDir = "/opt/somaticseq" #the location in the docker image

        File? classifierSNV
        File? classifierIndel
        String outputDir
        Reference reference
        File? inclusionRegion
        File? exclusionRegion
        IndexedBamFile bam
        File? mutect2VCF
        File? varscanVCF
        File? vardictVCF
        File? lofreqVCF
        File? scalpelVCF
        File? strelkaVCF

        Int threads = 1
    }

    command {
        ~{installDir}/somaticseq_parallel.py \
        ~{"--classifier-snv " + classifierSNV} \
        ~{"--classifier-indel " + classifierIndel} \
        --output-directory ~{outputDir} \
        --genome-reference ~{reference.fasta} \
        ~{"--inclusion-region " + inclusionRegion} \
        ~{"--exclusion-region " + exclusionRegion} \
        --threads ~{threads} \
        single \
        --bam-file ~{bam.file} \
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
        docker: "lethalfang/somaticseq:3.1.0"
    }
}

task ParallelSingleTrain {
    input {
        String installDir = "/opt/somaticseq" #the location in the docker image

        File truthSNV
        File truthIndel
        String outputDir
        Reference reference
        File? inclusionRegion
        File? exclusionRegion
        IndexedBamFile bam
        File? mutect2VCF
        File? varscanVCF
        File? vardictVCF
        File? lofreqVCF
        File? scalpelVCF
        File? strelkaVCF

        Int threads = 1
    }

    command {
        ~{installDir}/somaticseq_parallel.py \
        --somaticseq-train \
        --truth-snv ~{truthSNV} \
        --truth-indel ~{truthIndel} \
        --output-directory ~{outputDir} \
        --genome-reference ~{reference.fasta} \
        ~{"--inclusion-region " + inclusionRegion} \
        ~{"--exclusion-region " + exclusionRegion} \
        --threads ~{threads} \
        single \
        --bam-file ~{bam.file} \
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
        docker: "lethalfang/somaticseq:3.1.0"
    }
}
