version 1.0

task MultiQC {
    input {
        String? preCommand
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
        File? fileList
        Array[String] exclude
        Array[String] module
        Boolean dataDir = false
        Boolean noDataDir = false
        String? dataFormat
        Boolean zipDataDir = false
        Boolean export = false
        Boolean flat = false
        Boolean interactive = false
        Boolean lint = false
        Boolean pdf = false
        Boolean noMegaQCUpload = true
        File? config  # A directory
        String? clConfig
        Boolean verbose  = false
        Boolean quiet = false
    }
    command {
        set -e -o pipefail
        ~{preCommand}
        multiqc \
        ~{true="--force" false="" force} \
        ~{true="--dirs" false="" dirs} \
        ~{"--dirs-depth " + dirsDepth} \
        ~{true="--fullnames" false="" fullNames} \
        ~{"--title " + title } \
        ~{"--comment " + comment } \
        ~{"--filename " + fileName} \
        ~{
    }
    output {}
}