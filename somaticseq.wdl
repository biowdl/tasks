version 1.0

task ParallelPaired {
    input {
        String installDir = "/opt/somaticseq" #the location in the docker image

        File? classifierSNV
        File? classifierIndel
        String outputDir
        File referenceFasta
        File referenceFastaFai
        File? inclusionRegion
        File? exclusionRegion
        File tumorBam
        File tumorBamIndex
        File normalBam
        File normalBamIndex
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
        String dockerImage = "lethalfang/somaticseq:3.1.0"
    }

    command {
        ~{installDir}/somaticseq_parallel.py \
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
        docker: dockerImage
    }
}

task ParallelPairedTrain {
    input {
        String installDir = "/opt/somaticseq" #the location in the docker image

        File truthSNV
        File truthIndel
        String outputDir
        File referenceFasta
        File referenceFastaFai
        File? inclusionRegion
        File? exclusionRegion
        File tumorBam
        File tumorBamIndex
        File normalBam
        File normalBamIndex
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
        String dockerImage = "lethalfang/somaticseq:3.1.0"
    }

    command {
        ~{installDir}/somaticseq_parallel.py \
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
        docker: dockerImage
    }
}

task ParallelSingle {
    input {
        String installDir = "/opt/somaticseq" #the location in the docker image

        File? classifierSNV
        File? classifierIndel
        String outputDir
        File referenceFasta
        File referenceFastaFai
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

        Int threads = 1
        String dockerImage = "lethalfang/somaticseq:3.1.0"
    }

    command {
        ~{installDir}/somaticseq_parallel.py \
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
        docker: dockerImage
    }
}

task ParallelSingleTrain {
    input {
        String installDir = "/opt/somaticseq" #the location in the docker image

        File truthSNV
        File truthIndel
        String outputDir
        File referenceFasta
        File referenceFastaFai
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

        Int threads = 1
        String dockerImage = "lethalfang/somaticseq:3.1.0"
    }

    command {
        ~{installDir}/somaticseq_parallel.py \
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
        docker: dockerImage
    }
}

task ModifyStrelka {
    input {
        String installDir = "/opt/somaticseq/vcfModifier" #the location in the docker image

        File strelkaVCF
        String outputVCFName = basename(strelkaVCF, ".gz")

        Int threads = 1
        String dockerImage = "lethalfang/somaticseq:3.1.0"
    }

    command {
        set -e

        ~{installDir}/modify_Strelka.py \
        -infile ~{strelkaVCF} \
        -outfile "modified_strelka.vcf"

        first_FORMAT_line_num=$(grep -n -m 1 '##FORMAT' "modified_strelka.vcf" | cut -d : -f 1)
        sed "$first_FORMAT_line_num"'i##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' "modified_strelka.vcf" > ~{outputVCFName}
    }

    output {
        File outputVcf = outputVCFName
    }

    runtime {
        cpu: threads
        docker: dockerImage
    }
}
