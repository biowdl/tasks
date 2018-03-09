task unicycler {
    String? preCommand
    File short1
    File short2
    File? unpaired
    File? long
    String out
    Int? minFastaLength
    Int? keep
    Boolean? vcf
    Int? threads = 1
    Int? memory = 4
    String? mode
    Int? linearSeqs
    Int? verbosity
    command {
        set -e -o pipefail
        mkdir -p ${out}
        ${preCommand}
        unicycler \
        --short1 ${short1}
        --short2 ${short2}

    }

    runtime {
        cpu: select_first([threads])
        memory: select_first([memory])
    }
}