version 1.0

import "common.wdl" as common

task SomaticSeqWrapper {
    input {
        String? preCommand
        String? installDir

        String outputDir
        Reference ref
        File? inclusionRegion
        File? exclusionRegion
        File tumorBam
        File tumorIndex
        File normalBam
        File normalIndex
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
    }

    String toolCommand = if defined(installDir)
        then installDir + "/SomaticSeq.Wrapper.sh"
        else "SomaticSeq.Wrapper.sh"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        --output-dir ~{outputDir} \
        --genome-reference ~{ref.fasta} \
        ~{"--inclusion-region " +  inclusionRegion} \
        ~{"--exclusion-region " + exclusionRegion} \
        --tumor-bam ~{tumorBam} \
        --normal-bam ~{normalBam} \
        ~{"--mutect2 " + mutect2VCF} \
        ~{"--varscan-snv " + varscanSNV} \
        ~{"--varscan-indel " + varscanIndel} \
        ~{"--jsm " + jsmVCF} \
        ~{"--sniper " + somaticsniperVCF} \
        ~{"--vardict " + vardictVCF} \
        ~{"--muse " + museVCF} \
        ~{"--lofreq-snv " + lofreqSNV} \
        ~{"--lofreq-indel " + lofreqIndel} \
        ~{"--scalpel " + scalpelVCF} \
        ~{"--strelka-snv " + strelkaSNV} \
        ~{"--strelka-indel " + strelkaIndel}
    }

    output {
        File consensusIndels = outputDir + "/Consensus.sINDEL.vcf"
        File consensusSNV = outputDir + "/Consensus.sSNV.vcf"
        File ensembleIndels = outputDir + "/Ensemble.sINDEL.tsv"
        File ensembleSNV = outputDir + "/Ensemble.sSNV.tsv"
    }
}

task SsSomaticSeqWrapper {
    input {
        String? preCommand
        String? installDir

        String outputDir
        Reference ref
        File? inclusionRegion
        File? exclusionRegion
        File bam
        File bamIndex
        File? mutect2VCF
        File? varscanVCF
        File? vardictVCF
        File? lofreqVCF
        File? scalpelVCF
        File? strelkaVCF
    }

    String toolCommand = if defined(installDir)
        then installDir + "/ssSomaticSeq.Wrapper.sh"
        else "ssSomaticSeq.Wrapper.sh"

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{toolCommand} \
        --output-dir ~{outputDir} \
        --genome-reference ~{ref.fasta} \
        ~{"--inclusion-region " +  inclusionRegion} \
        ~{"--exclusion-region " + exclusionRegion} \
        --in-bam ~{bam} \
        ~{"--mutect2 " + mutect2VCF} \
        ~{"--varscan " + varscanVCF} \
        ~{"--vardict " + vardictVCF} \
        ~{"--lofreq " + lofreqVCF} \
        ~{"--scalpel " + scalpelVCF} \
        ~{"--strelka " + strelkaVCF}
    }

    output {
        File consensusIndels = outputDir + "/Consensus.ssINDEL.vcf"
        File consensusSNV = outputDir + "/Consensus.ssSNV.vcf"
        File ensembleIndels = outputDir + "/Ensemble.ssINDEL.tsv"
        File ensembleSNV = outputDir + "/Ensemble.ssSNV.tsv"
    }
}