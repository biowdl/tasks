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
    ${"--outdir " + outdirPath}
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
        String reportDir = outdirPath + "/" + sub(basename(seqFile),"\..*$","_fastqc")
        File rawReport = reportDir + "/fastqc_data.txt"
        File htmlReport = reportDir + "fastqc_report.html"
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
        cpu: select_first(threads)
    }
}
