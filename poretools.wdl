task fastq {
    Array[File]+ files # Files should exist! Also accepts multiple directories (unlike poretools).
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

    (
    # Allow for multiple directory input by looping over files
    for file in ${sep=" " files}
    do
        poretools fastq \
        ${"--type " + type} \
        ${"--start " + start } \
        ${"--end " + end } \
        ${"--min-length " + minLength } \
        ${"--max-length " + maxLength } \
        ${true="--high-quality" false="" highQuality} \
        ${true="--normal-quality" false="" normalQuality} \
        ${"--group " + group} \
        $file
    done
    ) ${true="| gzip " false="" gzip}> ${outputFile}
    }

    output {
    File fastq = outputFile
    }
}