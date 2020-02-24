version 1.0

task CreateSampleGetCount {
    input {
        File countTable
        String shinyDir = "/shiny"
        String dockerImage = "tomkuipers1402/shiny-py"

        String memory = "6G"
    }

    command {
        set -e
        mkdir -p ${shinyDir}
        sampleSheet.py \
            -i ${countTable} \
            -o ${shinyDir}
        cp ${countTable} ${shinyDir}/allCounts.tsv
    }

    output {
        String shinyOutput = shinyDir
    }

    runtime {
        memory: memory
        docker: dockerImage
    }
}

task CreateAnnotation {
    input {
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        File? referenceGtfFile
        String shinyDir = "/shiny"
        String dockerImage = "tomkuipers1402/shiny-py"

        String memory = "6G"
    }

    command {
        set -e
        mkdir -p ${shinyDir}
        annoGenerator.py \
            -r ${referenceFasta} \
            -g ${referenceGtfFile} \
            -o ${shinyDir}
    }

    output {
        String shinyOutput = shinyDir
    }

    runtime {
        memory: memory
        docker: dockerImage
    }
}
