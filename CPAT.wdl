version 1.0

task CPAT {
    input {
        File gene
        String outFilePath
        File hex
        File logitModel
        File? referenceGenome
        File? referenceGenomeIndex  # Should be added as input if
        # CPAT should not index the reference genome.
        Array[String]? startCodons
        Array[String]? stopCodons
        String dockerImage = "biocontainers/cpat:v1.2.4_cv1"
    }

    # Some WDL magic in the command section to properly output the start and stopcodons to the command.
    # select_first is needed in order to convert the optional arrays to non-optionals.
    command {
        set -e
        mkdir -p "$(dirname ~{outFilePath})"
        cpat.py \
        --gene ~{gene} \
        --outfile ~{outFilePath} \
        --hex ~{hex} \
        --logitModel ~{logitModel} \
        ~{"--ref " + referenceGenome} \
        ~{true="--start" false="" defined(startCodons)} ~{sep="," select_first([startCodons, [""]])} \
        ~{true="--stop" false="" defined(stopCodons)} ~{sep="," select_first([stopCodons, [""]])}
    }

    output {
        File outFile = outFilePath
    }

    runtime {
        docker: dockerImage
    }

    parameter_meta {
        gene: {description: "Equivalent to CPAT's `--gene` option.", category: "required"}
        outFilePath: {description: "Equivalent to CPAT's `--outfile` option.", category: "required"}
        hex: {description: "Equivalent to CPAT's `--hex` option.", category: "required"}
        logitModel: {description: "Equivalent to CPAT's `--logitModel` option.", category: "required"}
        referenceGenome: {description: "Equivalent to CPAT's `--ref` option.", category: "advanced"}
        referenceGenomeIndex: {description: "The index of the reference. Should be added as input if CPAT should not index the reference genome.",
                               category: "advanced"}
        startCodons: {description: "Equivalent to CPAT's `--start` option.", category: "advanced"}
        stopCodons: {description: "Equivalent to CPAT's `--stop` option.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

# There is also make_hexamer_tab.py and make_logitModel.py
# that can be added as tasks here.