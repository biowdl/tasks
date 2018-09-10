version 1.0

import "common.wdl" as common

task Mem {
    input {
        String? preCommand
        FastqPair inputFastq
        BwaIndex bwaIndex
        String outputPath
        String? readgroup

        String? picardJar

        Int threads = 1
        Int memory = 8
        Int picardMemory = 4
    }

    String picardPrefix = if defined(picardJar)
        then "java -Xmx" + picardMemory + "G -jar " + picardJar
        else "picard -Xmx" + picardMemory + "G"

    # Post alt script from bwa
    String altCommand = if (defined(bwaIndex.altIndex)) then "| bwa-postalt " + bwaIndex.altIndex else ""

    # setNmMdAndUqTags is only required if alt sequences are added
    String setNmMdAndUqTagsCommand = picardPrefix + " SetNmMdAndUqTags " +
                                             " INPUT=/dev/stdin OUTPUT=" + outputPath +
                                             " CREATE_INDEX=true" +
                                             " R=" + bwaIndex.fastaFile

    String sortSamCommand = picardPrefix + " SortSam " +
                 " INPUT=/dev/stdin SORT_ORDER=coordinate " +
                 if(defined(bwaIndex.altIndex)) then " OUTPUT=/dev/stdout "
                 else " OUTPUT=" + outputPath + " CREATE_INDEX=true "

    String picardCommand = if (defined(bwaIndex.altIndex)) then sortSamCommand + " | " + setNmMdAndUqTagsCommand
    else sortSamCommand

    String readgroupArg = if (defined(readgroup)) then "-R '" + readgroup + "'" else ""

    command {
        set -e -o pipefail
        mkdir -p $(dirname ~{outputPath})
        ~{preCommand}
        bwa mem ~{"-t " + threads} \
        ~{readgroupArg} \
        ~{bwaIndex.fastaFile} \
        ~{inputFastq.R1} \
        ~{inputFastq.R2} \
        ~{altCommand} \
        | ~{picardCommand}
    }

    output {
        IndexedBamFile bamFile = object {
          file: outputPath,
          index: sub(outputPath, ".bam$", ".bai")
        }
    }

    runtime{
        cpu: threads
        memory: memory + picardMemory + picardMemory
    }
}

task Index {
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
            indexFiles: [outputFile + ".bwt",outputFile + ".pac",outputFile + ".sa",outputFile + ".amb",outputFile + ".ann"]
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
    File? altIndex
}
