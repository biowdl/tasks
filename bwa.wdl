task mem {
    File read1
    File indexBase
    Array[File] indexFiles # Not used by the command, but needed so cromwell does provide the proper links.
    File? read2
    String? outputFile = "aligned.bam"
    String? preCommand
    Int? threads = 1
    String? memory = "4G"
    Int? minimumSeedLength
    Int? w
    Int? d
    Float? r
    Int? y
    Int? c
    Int? D
    Int? m
    Int? W
    Boolean? skipMateRescue
    Boolean? skipPairing
    Int? matchScore
    Int? mismatchPenalty
    String? gapOpenPenalty
    String? gapExtensionPenalty
    String? endClippingPenalty
    String? unpairedReadPairPenalty
    String? readType
    Boolean? smartPairing
    String? readGroupHeaderLine
    String? H
    Boolean? j
    Boolean? five
    Boolean? q
    Int? K
    Int? minimumOutputScore
    String? h
    Boolean? a
    Boolean? appendComment
    Boolean? V
    Boolean? Y
    Boolean? M
    String? I
    parameter_meta {
        referenceFiles: "Should contain all the index files from the index task"
    }
    command {
        set -e -o pipefail
        ${"mkdir -p $(dirname " + outputFile + ")"}
        ${preCommand}
        bwa mem \
        ${"-t " + threads } \
        ${"-k " + minimumSeedLength } \
        ${"-w " + w } \
        ${"-d " + d } \
        ${"-r " + r } \
        ${"-y " + y } \
        ${"-c " + c } \
        ${"-D " + D } \
        ${"-W " + W } \
        ${"-m " + m } \
        ${true="-s " false="" skipMateRescue } \
        ${true="-P " false="" skipPairing } \
        ${"-A " + matchScore } \
        ${"-B " + mismatchPenalty } \
        ${"-O " + gapOpenPenalty } \
        ${"-E " + gapExtensionPenalty } \
        ${"-L " + endClippingPenalty } \
        ${"-U " + unpairedReadPairPenalty } \
        ${"-x " + readType } \
        ${true="-p " false="" smartPairing} \
        ${"-r " + readGroupHeaderLine} \
        ${"-H " + H } \
        ${"-o " + outputFile } \
        ${true="-j" false="" j } \
        ${true="-5" false="" five } \
        ${true="-q" false="" q } \
        ${"-K " + K } \
        ${"-T " + minimumOutputScore } \
        ${"-h " + h } \
        ${true="-a" false="" a } \
        ${true="-C" false="" appendComment } \
        ${true="-V" false="" V } \
        ${true="-Y" false="" Y } \
        ${true="-M" false="" M } \
        ${"-I " + I } \
        ${indexBase} ${read1} ${read2} \
        | picard SortSam CREATE_INDEX=TRUE TMP_DIR=null \
        INPUT=/dev/stdin SORT_ORDER=coordinate OUTPUT=${outputFile}
    }
    output {
        File bamFile = select_first([outputFile])
        File bamIndex = sub(bamFile, ".bam$", ".bai")
    }
    runtime {
        cpu: select_first([threads])
        memory: select_first([memory])
        }
}

task index {
    File fasta
    String? preCommand
    String? constructionAlgorithm
    Int? blockSize
    String? outputDir
    String fastaFilename = basename(fasta)

    command {
        set -e -o pipefail
        ${"mkdir -p " + outputDir}
        ${preCommand}
        ln -sf ${fasta} ${outputDir + "/"}${fastaFilename}
        bwa index \
        ${"-a " + constructionAlgorithm} \
        ${"-b" + blockSize} \
        ${outputDir + "/"}${fastaFilename}
    }

    output {
        File indexBase = if (defined(outputDir)) then select_first([outputDir]) + "/" + fastaFilename else fastaFilename
        File indexedFasta = indexBase
        Array[File] indexFiles = [indexBase + ".bwt",indexBase + ".pac",indexBase + ".sa",indexBase + ".amb",indexBase + ".ann"]
    }
    parameter_meta {
        fasta: "Fasta file to be indexed"
        constructionAlgorithm: "-a STR    BWT construction algorithm: bwtsw, is or rb2 [auto]"
        blockSize: "-b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]"
        outputDir: "index will be created in this output directory"
    }
}