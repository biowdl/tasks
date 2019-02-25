version 1.0

task MultiQC {
    input {
        String? preCommand
        File analysisDirectory
        Array[File] dependencies   # This must be used in order to run multiqc after these tasks.
        Boolean force = false
        Boolean dirs = false
        Int? dirsDepth
        Boolean fullNames = false
        String? title
        String? comment
        String? fileName
        String outDir = "."
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

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p ~{outDir}
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

    String reportFilename = if (defined(fileName))
        then sub(select_first([fileName]), "\.html$", "")
        else "multiqc"

    output {
        File multiqcReport = outDir + "/" + reportFilename + "_report.html"
        File multiqcDataDir = outDir + "/" +reportFilename + "_data"
    }
}
