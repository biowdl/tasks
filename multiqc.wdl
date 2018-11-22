version 1.0

task MultiQC {
    input {
        String? preCommand
        File analysisDirectory
        Boolean force = false
        Boolean dirs = false
        Int? dirsDepth
        Boolean fullNames = false
        String? title
        String? comment
        String? fileName
        String? outDir
        String? template
        String? tag
        String? ignore
        String? ignoreSamples
        Boolean ignoreSymlinks = false
        File? sampleNames
        File? fileList
        Array[String]+? exclude
        Array[String]+? module
        Boolean dataDir = false
        Boolean noDataDir = false
        String? dataFormat
        Boolean zipDataDir = false
        Boolean export = false
        Boolean flat = false
        Boolean interactive = false
        Boolean lint = false
        Boolean pdf = false
        Boolean megaQCUpload = false # This must be actively enabled in my opinion. The tools default is to upload.
        File? config  # A directory
        String? clConfig
        Boolean verbose  = false
        Boolean quiet = false
    }

    String outputDir = if defined(outDir) then select_first([outDir]) else "."

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{if defined(outDir) then "mkdir -p " + outputDir else ""}
        multiqc \
        ~{true="--force" false="" force} \
        ~{true="--dirs" false="" dirs} \
        ~{"--dirs-depth " + dirsDepth} \
        ~{true="--fullnames" false="" fullNames} \
        ~{"--title " + title} \
        ~{"--comment " + comment} \
        ~{"--filename " + fileName} \
        ~{"--outdir " + outDir} \
        ~{"--template " + template} \
        ~{"--tag " + tag} \
        ~{"--ignore " + ignore} \
        ~{"--ignore-samples" + ignoreSamples} \
        ~{true="--ignore-symlinks" false="" ignoreSymlinks} \
        ~{"--sample-names " + sampleNames} \
        ~{"--file-list " + fileList} \
        ~{true="--exclude " false="" defined(exclude)}~{sep=" --exclude " exclude} \
        ~{true="--module " false="" defined(module)}~{sep=" --module " module} \
        ~{true="--data-dir" false="" dataDir} \
        ~{true="--no-data-dir" false="" noDataDir} \
        ~{"--data-format " + dataFormat} \
        ~{true="--zip-data-dir" false="" zipDataDir} \
        ~{true="--export" false="" export} \
        ~{true="--flat" false="" flat} \
        ~{true="--interactive" false="" interactive} \
        ~{true="--lint" false="" lint} \
        ~{true="--pdf" false="" pdf} \
        ~{false="--no-megaqc-upload" true="" megaQCUpload} \
        ~{"--config " + config} \
        ~{"--cl-config " + clConfig } \
        ~{analysisDirectory}
    }

    String reportFilename = if (defined(fileName)) then sub(select_first([fileName]), "\.html$", "") else "multi_qc"
    output {
        File multiqcReport = outputDir + "/" + reportFilename + ".html"
        File multiqcDataDir = outputDir + "/" +reportFilename + "_data"
        File multigcDataJson = multiqcDataDir + "/multiqc_data.json"
    }
}