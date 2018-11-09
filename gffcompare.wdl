version 1.0

task GffCompare {
    input {
        String? preCommand
        File? inputGtfList
        Array[File]+? inputGtfFiles
        File referenceAnnotation
        String? outputDir
        String outPrefix = "gffcmp" # gffcmp is the default used by the program as well. This needs to be
        # defined in order for the output values to be consistent and correct.
        File? genomeSequences
        Int? maxDistanceFreeEndsTerminalExons
        Int? maxDistanceGroupingTranscriptStartSites
        String? namePrefix
        Boolean C = false
        Boolean A = false
        Boolean X = false
        Boolean K = false
        Boolean snCorrection = false
        Boolean precisionCorrection = false
        Boolean discardSingleExonTransfragsAndReferenceTranscripts = false
        Boolean discardSingleExonReferenceTranscripts = false
        Boolean noTmap = false
        Boolean verbose = false
        Boolean debugMode = false
    }
    # This allows for the creation of output directories"
    String dirPrefix= if defined(outputDir) then outputDir + "/" else ""
    String totalPrefix = dirPrefix + outPrefix

    parameter_meta {}

    command {
        set -e
        ~{preCommand}
        ~{"mkdir -p " + outputDir}
        gffcompare \
        -r ~{referenceAnnotation} \
        ~{"-o " + totalPrefix } \
        ~{"-s " + genomeSequences} \
        ~{"-e " + maxDistanceFreeEndsTerminalExons} \
        ~{"-d " + maxDistanceGroupingTranscriptStartSites} \
        ~{"-p " + namePrefix} \
        ~{true="-C" false="" C} \
        ~{true="-A" false="" A} \
        ~{true="-X" false="" X} \
        ~{true="-K" false="" K} \
        ~{true="-R" false="" snCorrection} \
        ~{true="-Q" false="" precisionCorrection} \
        ~{true="-M" false="" discardSingleExonTransfragsAndReferenceTranscripts} \
        ~{true="-N" false="" discardSingleExonReferenceTranscripts} \
        ~{true="-T" false="" noTmap} \
        ~{true="-V" false="" verbose} \
        ~{true="D" false="" debugMode} \
        ~{"-i " + inputGtfList} \
        ~{sep=" " + inputGtfFiles}
    }

    output {}
}