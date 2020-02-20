version 1.0

task CollectAllFiles {
    input {
        File inputSheet
        File countTable
        File? annotation
        String outputSheet = "sampleSheet.tsv"
        String shinyDir = "/shiny"
        String dockerImage = "tomkuipers1402/pandas"
    }

    command {
        set -e
        mkdir -p ${shinyDir}
        sampleSheet.py \
            -i ${inputSheet} \
            -o ${shinyDir}/${outputSheet}
        ${shinyDir}/${countTable}
        ${shinyDir}/${annotation}
    }

    output {
        String shinyOutput = shinyDir
    }

    runtime {
        docker: dockerImage
    }
}
