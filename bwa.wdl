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
        # A mulled container is needed to have both picard and bwa in one container.
        # This container contains: picard (2.18.7), bwa (0.7.17-r1188)
        String dockerImage = "quay.io/biocontainers/mulled-v2-002f51ea92721407ef440b921fb5940f424be842:43ec6124f9f4f875515f9548733b8b4e5fed9aa6-0"
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
        File outputBamIndex = sub(outputPath, "\.bam$", ".bai")
    }

    runtime{
        cpu: threads
        memory: memory + picardMemory + picardMemory
        docker: dockerImage
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

task Kit {
    input {
        File read1
        File? read2
        BwaIndex bwaIndex
        String outputPrefix
        String? readgroup
        Boolean sixtyFour = false

        Int threads = 2
        Int sortThreads = 2
        Int memory = 10
        String dockerImage = "biocontainers/bwakit:v0.7.15_cv1"
    }

    command {
        set -e
        bwa mem \
          -t ~{threads} \
          ~{"-R '" + readgroup}~{true="'" false="" defined(readgroup)} \
          ~{bwaIndex.fastaFile} \
          ~{read1} \
          ~{read2} \
          2> ~{outputPrefix}.log.bwamem | \
        k8 /opt/conda/bin/bwa-postalt.js \
          -p ~{outputPrefix}.hla \
          ~{bwaIndex.fastaFile}~{true=".64.alt" false=".alt" sixtyFour} | \
        samtools sort \
          -@ ~{sortThreads} \
          -m1G \
          - \
          -o ~{outputPrefix}.aln.bam
        samtools index ~{outputPrefix}.aln.bam ~{outputPrefix}.aln.bai
    }

    output {
        File outputBam = outputPrefix + ".aln.bam"
        File outputBamIndex = outputPrefix + ".aln.bai"
    }

    runtime{
        cpu: threads + sortThreads
        memory: memory
        docker: dockerImage
    }
}

struct BwaIndex {
    File fastaFile
    Array[File] indexFiles
}
