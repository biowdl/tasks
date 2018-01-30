task mem {
    File read1
    File referenceFile
    File? read2
    String? outputFile
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
    command {
        set -e -o pipefail
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
        ${referenceFile} ${read1} ${read2}
    }
    output {
        File samFile = if defined(outputFile) then outputFile else stdout()
    }
    runtime {
        cpu: select_first([threads])
        memory: select_first([memory])
        }
}