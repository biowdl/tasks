version 1.0

import "common.wdl"

task BgzipAndIndex {
    input {
        File inputFile
        String outputDir
        String type = "vcf"
    }

    String outputGz = outputDir + "/" + basename(inputFile) + ".gz"

    command {
        bgzip -c ~{inputFile} > ~{outputGz}
        tabix ~{outputGz} -p ~{type}
    }

    output {
        File compressed = outputGz
        File index = outputGz + ".tbi"
    }
}

task Index {
    input {
        String? preCommand
        File bamFile
        String bamIndexPath

        String dockerTag = "1.8--h46bd0b3_5"
    }

    command {
        set -e -o pipefail
        ~{preCommand}
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
        String? preCommand
        Array[File]+ bamFiles
        String outputBamPath
    }

    command {
        set -e -o pipefail
        ~{preCommand}
        samtools merge ~{outputBamPath} ~{sep=' ' bamFiles}
    }

    output {
        File outputBam = outputBamPath
    }
}

task Markdup {
    input {
        String? preCommand
        File inputBam
        String outputBamPath
    }

    command {
        set -e -o pipefail
        ~{preCommand}
        samtools markdup ~{inputBam} ~{outputBamPath}
    }

    output {
        File outputBam = outputBamPath
    }
}

task Flagstat {
    input {
        String? preCommand
        File inputBam
        String outputPath
    }

    command {
        set -e -o pipefail
        ~{preCommand}
        mkdir -p $(dirname ~{outputPath})
        samtools flagstat ~{inputBam} > ~{outputPath}
    }

    output {
        File flagstat = outputPath
    }
}

task Fastq {
    input {
        String? preCommand
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
    }

    command {
        ~{preCommand}
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
    }

    parameter_meta {
        preCommand: "A command that is run before the task. Can be used to activate environments"
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
        String inputFile
        String type = "vcf"
    }

    command {
        tabix ~{inputFile} -p ~{type}
    }

    output {
        File index = inputFile + ".tbi"
    }
}

task View {
    input {
        String? preCommand
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
    }

    command {
        set -e -o pipefail
        ~{preCommand}
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
    }
}
