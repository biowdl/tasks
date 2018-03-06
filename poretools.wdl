task fastq {
    Array[File] files
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
    ${sep=" " files} > ${outputFile}
    }

    output {
    File fastq = outputFile
    }
}