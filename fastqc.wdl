task fastqc {
    File seqFile
    String outdirPath
    String? condaEnv
    Boolean? casava
    Boolean? nano
    Boolean? noFilter
    Boolean? extract = true
    Boolean? nogroup
    Int? minLength
    String? format
    Int? threads = 1
    File? contaminants
    File? adapters
    File? limits
    Int? kmers
    String? dir

    command {
    set -e -o pipefail
    ${"source activate " + condaEnv}
    mkdir -p ${outdirPath}
    fastqc \
    ${"--outdir " + outdirPath} \
    ${true="--casava" false="" casava} \
    ${true="--nano" false="" nano} \
    ${true="--nofilter" false="" noFilter} \
    ${true="--extract" false="" extract} \
    ${true="--nogroup" false="" nogroup} \
    ${"--min_length " + minLength } \
    ${"--format " + format} \
    ${"--threads " + threads} \
    ${"--contaminants " + contaminants} \
    ${"--adapters " + adapters} \
    ${"--limits " + limits} \
    ${"--kmers " + kmers} \
    ${"--dir " + dir} \
    ${seqFile}

    }

    output {
        # Apparently, the escape character needs to be escaped in regexes.
        # This regex chops of the extension and replaces it with _fastqc for the reportdir.
        # Just as fastqc does it.
        String reportDir = outdirPath + "/" + sub(basename(seqFile), "\\.[^\\.]*$", "_fastqc")
        File rawReport = reportDir + "/fastqc_data.txt"
        File htmlReport = reportDir + "/fastqc_report.html"
        File summary = reportDir + "/summary.txt"
        File adapterContent = reportDir + "/Images/adapter_content.png"
        File duplicationLevels = reportDir + "/Images/duplication_levels.png"
        File perBaseNContent = reportDir + "/Images/per_base_n_content.png"
        File perBaseQuality = reportDir + "/Images/per_base_quality.png"
        File perBaseSequenceContent = reportDir + "/Images/per_base_sequence_content.png"
        File perSequenceGCContent = reportDir + "/Images/per_sequence_gc_content.png"
        File perSequenceQuality = reportDir + "/Images/per_sequence_quality.png"
        File perTileQuality = reportDir + "/Images/per_tile_quality.png"
        File sequenceLengthDistribution = reportDir + "/Images/sequence_length_distribution.png"
    }

    runtime {
        cpu: select_first([threads])
    }
}

task extractAdapters {
    File extractAdaptersFastqcJar
    File inputFile
    String outputDir
    String? adapterOutputFilePath = outputDir + "/adapter.list"
    String? contamsOutputFilePath = outputDir + "/contaminations.list"
    Boolean? skipContams
    File? knownContamFile
    File? knownAdapterFile
    Float? adapterCutoff
    Boolean? outputAsFasta
    command {
    set -e
    mkdir -p ${outputDir}
    java -jar ${extractAdaptersFastqcJar} \
    --inputFile ${inputFile} \
    ${"--adapterOutputFile " + adapterOutputFilePath } \
    ${"--contamsOutputFile " + contamsOutputFilePath } \
    ${"--knownContamFile " + knownContamFile} \
    ${"--knownAdapterFile " + knownAdapterFile} \
    ${"--adapterCutoff " + adapterCutoff} \
    ${true="--skipContams" false="" skipContams} \
    ${true="--outputAsFasta" false="" outputAsFasta}
    }

    output {
        File adapterOutputFile = select_first([adapterOutputFilePath])
        File contamsOutputFile = select_first([contamsOutputFilePath])
        Array[String] adapterList = read_lines(select_first([adapterOutputFilePath]))
        Array[String] contamsList = read_lines(select_first([contamsOutputFilePath]))
    }
}