version 1.0

import "common.wdl"

task BgzipAndIndex {
    input {
        File inputFile
        String outputDir
        String type = "vcf"

        String dockerTag = "0.2.6--ha92aebf_0"
    }

    String outputGz = outputDir + "/" + basename(inputFile) + ".gz"

    command {
        set -e
        bgzip -c ~{inputFile} > ~{outputGz}
        tabix ~{outputGz} -p ~{type}
    }

    output {
        File compressed = outputGz
        File index = outputGz + ".tbi"
    }

    runtime {
        docker: "quay.io/biocontainers/tabix:" + dockerTag
    }
}

task Index {
    input {
        File bamFile
        String bamIndexPath

        String dockerTag = "1.8--h46bd0b3_5"
    }

    command {
        samtools index ~{bamFile} ~{bamIndexPath}
    }

    output {
        IndexedBamFile outputBam = object {
          file: bamFile,
          index: bamIndexPath
        }
    }

    runtime {
        docker: "quay.io/biocontainers/samtools:" + dockerTag
    }
}

task Merge {
    input {
        Array[File]+ bamFiles
        String outputBamPath = "merged.bam"
        Boolean force = true

        String dockerTag = "1.8--h46bd0b3_5"
    }
    String indexPath = sub(outputBamPath, "\.bam$",".bai")

    command {
        set -e
        mkdir -p $(dirname ~{outputBamPath})
        samtools merge ~{true="-f" false="" force} ~{outputBamPath} ~{sep=' ' bamFiles}
        samtools index ~{outputBamPath} ~{indexPath}
    }

    output {
        File outputBam = outputBamPath
        File outputBamIndex = indexPath
    }

    runtime {
        docker: "quay.io/biocontainers/samtools:" + dockerTag
    }
}

task Markdup {
    input {
        File inputBam
        String outputBamPath

        String dockerTag = "1.8--h46bd0b3_5"
    }

    command {
        samtools markdup ~{inputBam} ~{outputBamPath}
    }

    output {
        File outputBam = outputBamPath
    }

    runtime {
        docker: "quay.io/biocontainers/samtools:" + dockerTag
    }
}

task Flagstat {
    input {
        File inputBam
        String outputPath

        String dockerTag = "1.8--h46bd0b3_5"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPath})
        samtools flagstat ~{inputBam} > ~{outputPath}
    }

    output {
        File flagstat = outputPath
    }

    runtime {
        docker: "quay.io/biocontainers/samtools:" + dockerTag
    }
}

task Fastq {
    input {
        File inputBam
        String outputRead1
        String? outputRead2
        String? outputRead0
        Int? includeFilter
        Int? excludeFilter
        Int? excludeSpecificFilter
        Boolean? appendReadNumber
        Boolean? outputQuality
        Int? compressionLevel

        Int threads = 1
        Int memory = 1
        String dockerTag = "1.8--h46bd0b3_5"
    }

    command {
        samtools fastq \
        ~{true="-1" false="-s" defined(outputRead2)} ~{outputRead1} \
        ~{"-2 " + outputRead2} \
        ~{"-0 " + outputRead0} \
        ~{"-f " + includeFilter} \
        ~{"-F " + excludeFilter} \
        ~{"-G " + excludeSpecificFilter} \
        ~{true="-N" false="-n" appendReadNumber} \
        ~{true="-O" false="" outputQuality} \
        ~{"-c " + compressionLevel} \
        ~{"--threads " + threads} \
        ~{inputBam}
    }

    output {
        File read1 = outputRead1
        File? read2 = outputRead2
        File? read0 = outputRead0
    }

    runtime {
        cpu: threads
        memory: memory
        docker: "quay.io/biocontainers/samtools:" + dockerTag
    }

    parameter_meta {
        inputBam: "The bam file to process."
        outputRead1: "If only outputRead1 is given '-s' flag is assumed. Else '-1'."
        includeFilter: "Include reads with ALL of these flags. Corresponds to '-f'"
        excludeFilter: "Exclude reads with ONE OR MORE of these flags. Corresponds to '-F'"
        excludeSpecificFilter: "Exclude reads with ALL of these flags. Corresponds to '-G'"
        appendReadNumber: "Append /1 and /2 to the read name, or don't. Corresponds to '-n/N"

    }
}

task Tabix {
    input {
        File inputFile
        String outputFilePath = "indexed.vcf.gz"
        String type = "vcf"
        String dockerTag = "0.2.6--ha92aebf_0"
    }
    # FIXME: It is better to do the indexing on VCF creation. Not in a separate task. With file localization this gets hairy fast.
    command {
        set -e
        mkdir -p $(dirname ~{outputFilePath})
        if [ ! -f ~{outputFilePath} ]
        then
            ln ~{inputFile} ~{outputFilePath}
        fi
        tabix ~{outputFilePath} -p ~{type}
    }

    output {
        File indexedFile = outputFilePath
        File index = outputFilePath + ".tbi"
    }

    runtime {
        docker: "quay.io/biocontainers/tabix:" + dockerTag
    }
}

task View {
    input {
        File inFile
        File? referenceFasta
        String outputFileName
        Boolean? outputBam
        Boolean? uncompressedBamOutput
        Int? includeFilter
        Int? excludeFilter
        Int? excludeSpecificFilter
        Int? MAPQthreshold

        Int threads = 1
        Int memory = 1
        String dockerTag = "1.8--h46bd0b3_5"
    }

    command {
        samtools view \
        ~{"-T " + referenceFasta} \
        ~{"-o " + outputFileName} \
        ~{true="-b " false="" outputBam} \
        ~{true="-u " false="" uncompressedBamOutput} \
        ~{"-f " + includeFilter} \
        ~{"-F " + excludeFilter} \
        ~{"-G " + excludeSpecificFilter} \
        ~{"-q " + MAPQthreshold} \
        ~{"--threads " + (threads - 1)} \
        ~{inFile}
    }

    output {
        File outputFile = outputFileName
    }

    runtime {
        cpu: threads
        memory: memory
        docker: "quay.io/biocontainers/samtools:" + dockerTag
    }
}
