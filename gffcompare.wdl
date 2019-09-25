version 1.0

task GffCompare {
    input {
        File? inputGtfList
        Array[File] inputGtfFiles
        File referenceAnnotation
        String? outputDir
        String outPrefix = "gffcmp" # gffcmp is the default used by the program as well. This
        # needs to be defined in order for the output values to be consistent and correct.
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

        String dockerImage = "quay.io/biocontainers/gffcompare:0.10.6--h2d50403_0"

        # This workaround only works in the input section.
        # Issue addressed at https://github.com/openwdl/wdl/pull/263
        File? noneFile # This is a wdl workaround. Please do not assign!
    }
    # This allows for the creation of output directories
    String dirPrefix= if defined(outputDir)
        then outputDir + "/"
        else ""
    String totalPrefix = dirPrefix + outPrefix

    command {
        set -e
        ~{"mkdir -p " + outputDir}
        gffcompare \
        -r ~{referenceAnnotation} \
        ~{"-o '" + totalPrefix + "'"} \
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
        ~{sep=" " inputGtfFiles}
    }

    # Output of gffcompare is not stable. It depends on the number of files in the input.
    Int noFilesGtfList = if defined(inputGtfList)
        then length(read_lines(select_first([inputGtfList])))
        else 0
    Int noInputFiles = length(inputGtfFiles)
    Boolean oneFile = (noFilesGtfList + noInputFiles) == 1
    String annotatedName = if oneFile
        then "annotated"
        else "combined"

    # Check if a redundant .gtf will be created
    Boolean createRedundant = C || A || X

    output {
        File annotated = totalPrefix + "." + annotatedName + ".gtf"
        File loci = totalPrefix + ".loci"
        File stats = totalPrefix + ".stats"
        File tracking = totalPrefix + ".tracking"
        # noneFile is not stable. Please replace this as soon as wdl spec allows
        File? redundant = if createRedundant
            then totalPrefix + ".redundant.gtf"
            else noneFile
        File? missedIntrons = if debugMode
            then totalPrefix + ".missed_introns.gtf"
            else noneFile
    }

    runtime {
       docker: dockerImage
    }
}