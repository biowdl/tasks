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

task GffCompare {
    input {
        Array[File] inputGtfFiles
        # gffcmp is the default used by the program as well. This needs to be
        # defined in order for the output values to be consistent and correct.
        String outPrefix = "gffcmp"
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

        File? inputGtfList
        File? referenceAnnotation
        String? outputDir
        File? genomeSequences
        Int? maxDistanceFreeEndsTerminalExons
        Int? maxDistanceGroupingTranscriptStartSites
        String? namePrefix

        String memory = "4G"
        Int timeMinutes = 1 + ceil(size(inputGtfFiles, "G") * 30)
        String dockerImage = "quay.io/biocontainers/gffcompare:0.10.6--h2d50403_0"

        # This workaround only works in the input section.
        # Issue addressed at https://github.com/openwdl/wdl/pull/263.
        File? noneFile # This is a wdl workaround. Please do not assign!
    }

    # This allows for the creation of output directories.
    String dirPrefix = if defined(outputDir)
        then select_first([outputDir]) + "/"
        else ""
    String totalPrefix = dirPrefix + outPrefix

    command {
        set -e
        ~{"mkdir -p " + outputDir}
        gffcompare \
        ~{"-r " + referenceAnnotation} \
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
    String annotatedName = if oneFile && defined(referenceAnnotation)
        then "annotated"
        else "combined"

    # Check if a redundant .gtf will be created.
    Boolean createRedundant = C || A || X

    output {
        # noneFile is not stable. Please replace this as soon as wdl spec allows.
        File annotated = totalPrefix + "." + annotatedName + ".gtf"
        File loci = totalPrefix + ".loci"
        File stats = totalPrefix + ".stats"
        File tracking = totalPrefix + ".tracking"
        Array[File] allFiles = select_all([annotated, loci, stats, tracking, redundant, missedIntrons])
        File? redundant = if createRedundant
            then totalPrefix + ".redundant.gtf"
            else noneFile
        File? missedIntrons = if debugMode
            then totalPrefix + ".missed_introns.gtf"
            else noneFile
    }

    runtime {
        memory: memory
       time_minutes: timeMinutes
       docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputGtfFiles: {description: "The input GTF files.", category: "required"}
        referenceAnnotation: {description: "The GTF file to compare with.", category: "required"}
        outPrefix: {description: "The prefix for the output.", category: "advanced"}
        C: {description: "Equivalent to gffcompare's `-C` flag.", category: "advanced"}
        A: {description: "Equivalent to gffcompare's `-A` flag.", category: "advanced"}
        X: {description: "Equivalent to gffcompare's `-X` flag.", category: "advanced"}
        K: {description: "Equivalent to gffcompare's `-K` flag.", category: "advanced"}
        snCorrection: {description: "Equivalent to gffcompare's `-R` flag.", category: "advanced"}
        precisionCorrection: {description: "Equivalent to gffcompare's `-Q` flag.", category: "advanced"}
        discardSingleExonTransfragsAndReferenceTranscripts: {description: "Equivalent to gffcompare's `-M` flag.", category: "advanced"}
        discardSingleExonReferenceTranscripts: {description: "Equivalent to gffcompare's `-N` flag.", category: "advanced"}
        noTmap: {description: "Equivalent to gffcompare's `-T` flag.", category: "advanced"}
        verbose: {description: "Equivalent to gffcompare's `-V` flag.", category: "advanced"}
        debugMode: {description: "Equivalent to gffcompare's `-D` flag.", category: "advanced"}
        inputGtfList: {description: "Equivalent to gffcompare's `-i` option.", category: "advanced"}
        outputDir: {description: "The location the output should be written.", category: "common"}
        genomeSequences: {description: "Equivalent to gffcompare's `-s` option.", category: "advanced"}
        maxDistanceFreeEndsTerminalExons: {description: "Equivalent to gffcompare's `-e` option.", category: "advanced"}
        maxDistanceGroupingTranscriptStartSites: {description: "Equivalent to gffcompare's `-d` option.", category: "advanced"}
        namePrefix: {description: "Equivalent to gffcompare's `-p` option.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        annotated: {description: "Annotated GTF file."}
        loci: {description: "File describing the processed loci."}
        stats: {description: "Various statistics related to the “accuracy” (or a measure of agreement) of the input transcripts when compared to reference annotation data."}
        tracking: {description: "File matching up transcripts between samples."}
        allFiles: {description: "A collection of all output files."}
        redundant: {description: "File containing duplicate/redundant transcripts."}
        missedIntrons: {description: "File denoting missed introns."}
    }

    meta {
        WDL_AID: {
            exclude: ["noneFile"]
        }
    }
}
