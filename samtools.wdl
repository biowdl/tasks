task Index {
    String? preCommand
    File bamFilePath
    String bamIndexPath

    command {
        set -e -o pipefail
        ${preCommand}
        samtools index ${bamFilePath} ${bamIndexPath}
    }

    output {
        File indexFile = bamIndexPath
    }
}

task Merge {
    String? preCommand
    Array[File]+ bamFiles
    String outputBamPath

    command {
        set -e -o pipefail
        ${preCommand}
        if [ ${length(bamFiles)} -gt 1 ]
          then
            samtools merge ${outputBamPath} ${sep=' ' bamFiles}
          else
            ln -sf ${sep=' ' bamFiles} ${outputBamPath}
        fi
    }

    output {
        File outputBam = outputBamPath
    }
}

task Markdup {
    String? preCommand
    File inputBam
    String outputBamPath

    command {
        set -e -o pipefail
        ${preCommand}
        samtools markdup ${inputBam} ${outputBamPath}
    }

    output {
        File outputBam = outputBamPath
    }
}

task Flagstat {
    String? preCommand
    File inputBam
    String outputPath

    command {
        set -e -o pipefail
        ${preCommand}
        samtools flagstat ${inputBam} > ${outputPath}
    }

    output {
        File flagstat = outputPath
    }
}

task fastq {
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
    Int? threads
    Int? memory
    Int totalThreads = select_first([threads, 1])

    command {
        ${preCommand}
        samtools fastq \
        ${true="-1" false="-s" defined(outputRead2)} ${outputRead1} \
        ${"-2 " + outputRead2} \
        ${"-0 " + outputRead0} \
        ${"-f " + includeFilter} \
        ${"-F " + excludeFilter} \
        ${"-G " + excludeSpecificFilter} \
        ${true="-N" false="-n" appendReadNumber} \
        ${true="-O" false="" outputQuality} \
        ${"-c " + compressionLevel} \
        ${"--threads " + totalThreads} \
        ${inputBam}
    }
    output {
        File read1 = outputRead1
        File? read2 = outputRead2
        File? read0 = outputRead0
    }
    runtime {
        cpu: totalThreads
        memory: select_first([memory, 1])
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

task view {
    String? preCommand
    File inFile
    File? referenceFasta
    String outputFileName
    Boolean? outputBam
    Boolean? uncompressedBamOutput
    Int? includeFilter
    Int? excludeFilter
    Int? excludeSpecificFilter
    Int? threads
    Int? memory

    command {
    set -e -o pipefail
    ${preCommand}
    samtools view \
    ${"-T " + referenceFasta} \
    ${"-o " + outputFileName} \
    ${true="-b " false="" outputBam} \
    ${true="-u " false="" uncompressedBamOutput} \
    ${"-f " + includeFilter} \
    ${"-F " + excludeFilter} \
    ${"-G " + excludeSpecificFilter} \
    ${"--threads " + threads - 1} \
    ${inFile}
    }

    output {
        File outputFile = outputFileName
    }
    runtime {
        cpu: select_first([threads, 1])
        memory: select_first([memory, 1])
    }
}
