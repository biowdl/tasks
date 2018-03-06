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
    ${true="--high-quality" false="" highQuality} \
    ${true="--normal-quality" false="" normalQuality} \
    ${"--group " + group} \
    ${sep=" " files} ${true="| gzip " false="" gzip}> ${outputFile}
    }

    output {
    File fastq = outputFile
    }
}