task fastq {
    Array[String]+ files # Array[String]+ instead of Array[File]+ to allow for globs
    String outputFile
    String? preCommand
    String? type
    String? start
    String? end
    Int? minLength
    Int? maxLength
    Boolean? highQuality
    Boolean? normalQuality
    String? group
    Boolean? gzip = true
    command {
    set -e -o pipefail
    mkdir -p $(dirname ${outputFile})
    ${preCommand}

    poretools fastq \
    ${"--type " + type} \
    ${"--start " + start } \
    ${"--end " + end } \
    ${"--min-length " + minLength } \
    ${"--max-length " + maxLength } \
    ${if highQuality then "--high-quality" else ""} \
    ${if normalQuality then "--normal-quality" else ""} \
    ${"--group " + group} \
    ${sep=" " files} ${if gzip then "| gzip " else ""}> ${outputFile}
    }

    output {
    File fastq = outputFile
    }
}