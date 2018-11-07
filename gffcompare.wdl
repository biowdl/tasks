version 1.0

task GffCompare {
    input {
        File? inputGtfList
        Array[File]+? inputGtfFiles
        File referenceAnnotation
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

    parameter_meta {}

    command {
        gffcompare \
        -r ~{referenceAnnotation} \
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