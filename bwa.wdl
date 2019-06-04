version 1.0

task Mem {
    input {
        File read1
        File? read2
        BwaIndex bwaIndex
        String outputPath
        String? readgroup

        String? picardJar

        Int threads = 2
        Int memory = 8
        Int picardMemory = 4
        String dockerTag = "43ec6124f9f4f875515f9548733b8b4e5fed9aa6-0"
    }

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{outputPath})
        bwa mem \
        ~{"-t " + threads} \
        ~{"-R '" + readgroup}~{true="'" false="" defined(readgroup)} \
        ~{bwaIndex.fastaFile} \
        ~{read1} \
        ~{read2} \
        | picard -Xmx~{picardMemory}G SortSam \
        INPUT=/dev/stdin \
        OUTPUT=~{outputPath} \
        SORT_ORDER=coordinate \
        CREATE_INDEX=true
    }

    output {
        File outputBam = outputPath
        File outputBamIndex = sub(outputPath, ".bam$", ".bai")
    }

    runtime{
        cpu: threads
        memory: memory + picardMemory + picardMemory
        # A mulled container is needed to have both picard and bwa in one container.
        # This container contains: picard (2.18.7), bwa (0.7.17-r1188)
        docker: "quay.io/biocontainers/mulled-v2-002f51ea92721407ef440b921fb5940f424be842:" +
            dockerTag
    }
}

task Index {
    # Since this task uses `ln` this is not stable or usable with containers
    input {
        File fasta
        String? preCommand
        String? constructionAlgorithm
        Int? blockSize
        String? outputDir
    }

    String fastaFilename = basename(fasta)
    String outputFile = if (defined(outputDir)) then outputDir + "/" + fastaFilename else fasta

    command {
        set -e -o pipefail
        ~{"mkdir -p " + outputDir}
        ~{preCommand}
        if [[ ! '~{outputDir}' =  '' ]]
        then
            ln -sf ~{fasta} ~{outputDir + "/"}~{fastaFilename}
        fi
        bwa index \
        ~{"-a " + constructionAlgorithm} \
        ~{"-b" + blockSize} \
        ~{outputFile}
    }

    output {
        BwaIndex outputIndex = object {
            fastaFile: outputFile,
            indexFiles: [outputFile + ".bwt",
                outputFile + ".pac",
                outputFile + ".sa",
                outputFile + ".amb",
                outputFile + ".ann"]
        }
    }

    parameter_meta {
        fasta: "Fasta file to be indexed"
        constructionAlgorithm: "-a STR    BWT construction algorithm: bwtsw, is or rb2 [auto]"
        blockSize: "-b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]"
        outputDir: "index will be created in this output directory"
    }
}

struct BwaIndex {
    File fastaFile
    Array[File] indexFiles
}
