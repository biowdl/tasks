task Bowtie_se {
    File reference
    File read1
    String sam
    String? precommand
    Int? threads
    String? otherOptions

    command {
        ## https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
        set -eo pipefail
        mkdir -p
        ${precommand}
        bowtie \
        ${otherOptions}
        ${"--threads" + threads} \
        ${reference} \
        ${read1} > \
        ${sam}
    }
}

task Bowtie_pe {
    File reference
    File read1
    File read2
    String sam
    String? precommand
    Int? threads
    String? otherOptions

    command {
        set -eo pipefail
        mkdir -p
        ${precommand}
        bowtie \
        ${"--threads" + threads} \
        ${reference} \
        -1 ${read1} |
        -2 ${read2}> \
        ${sam}
    }
}