version 1.0

task Stringtie {
    input {
        File bam
        File bamIndex
        File? referenceGtf
        Boolean skipNovelTranscripts = false
        String assembledTranscriptsFile
        Boolean? firstStranded
        Boolean? secondStranded
        String? geneAbundanceFile

        Int threads = 1
        String memory = "10G"
        String dockerImage = "quay.io/biocontainers/stringtie:1.3.4--py35_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{assembledTranscriptsFile})"
        stringtie \
        ~{"-p " + threads} \
        ~{"-G " + referenceGtf} \
        ~{true="-e" false="" skipNovelTranscripts} \
        ~{true="--rf" false="" firstStranded} \
        ~{true="--fr" false="" secondStranded} \
        -o ~{assembledTranscriptsFile} \
        ~{"-A " + geneAbundanceFile} \
        ~{bam}
    }

    output {
        File assembledTranscripts = assembledTranscriptsFile
        File? geneAbundance = geneAbundanceFile
    }

    runtime {
        cpu: threads
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        bam: {
            description: "The input BAM file.",
            category: "required"
        }
        bamIndex: {
            description: "The input BAM file's index.",
            category: "required"
        }
        referenceGtf: {
            description: "A reference GTF file to be used as guide.",
            category: "common"
        }
        skipNovelTranscripts: {
            description: "Whether new transcripts should be assembled or not.",
            category: "common"
        }
        assembledTranscriptsFile: {
            description: "Where the output of the assembly should be written.",
            category: "required"
        }
        firstStranded: {
            description: "Equivalent to the --rf flag of stringtie.",
            category: "required"
        }
        secondStranded: {
            description: "Equivalent to the --fr flag of stringtie.",
            category: "required"
        }
        geneAbundanceFile: {
            description: "Where the abundance file should be written.",
            category: "common"
        }
        threads: {
            description: "The number of threads to use.",
            category: "advanced"
        }
        memory: {
            description: "The amount of memory needed for this task in GB.",
            category: "advanced"
        }
        dockerImage: {
            description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
            category: "advanced"
        }
    }
}

task Merge {
    input {
        Array[File]+ gtfFiles
        String outputGtfPath
        File? guideGtf
        Int? minimumLength
        Float? minimumCoverage
        Float? minimumFPKM
        Float? minimumTPM
        Float? minimumIsoformFraction
        Boolean keepMergedTranscriptsWithRetainedIntrons = false
        String? label

        String memory = "10G"
        String dockerImage = "quay.io/biocontainers/stringtie:1.3.4--py35_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputGtfPath})"
        stringtie --merge \
        -o ~{outputGtfPath} \
        ~{"-G " + guideGtf} \
        ~{"-m " + minimumLength } \
        ~{"-c " + minimumCoverage} \
        ~{"-F " + minimumFPKM} \
        ~{"-T " + minimumTPM} \
        ~{"-f " + minimumIsoformFraction} \
        ~{true="-i" false="" keepMergedTranscriptsWithRetainedIntrons} \
        ~{"-l " + label} \
        ~{sep=" " gtfFiles}
    }

    output {
        File mergedGtfFile = outputGtfPath
    }

    runtime {
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        gtfFiles: {
            description: "The GTF files produced by stringtie.",
            category: "required"
        }
        outputGtfPath: {
            description: "Where the output should be written.",
            category: "required"
        }
        guideGtf: {
            description: "Equivalent to the -G option of 'stringtie --merge'.",
            category: "advanced"
        }
        minimumLength: {
            description: "Equivalent to the -m option of 'stringtie --merge'.",
            category: "advanced"
        }
        minimumCoverage: {
            description: "Equivalent to the -c option of 'stringtie --merge'.",
            category: "advanced"
        }
        minimumFPKM: {
            description: "Equivalent to the -F option of 'stringtie --merge'.",
            category: "advanced"
        }
        minimumTPM: {
            description: "Equivalent to the -T option of 'stringtie --merge'.",
            category: "advanced"
        }
        minimumIsoformFraction: {
            description: "Equivalent to the -f option of 'stringtie --merge'.",
            category: "advanced"
        }
        keepMergedTranscriptsWithRetainedIntrons: {
            description: "Equivalent to the -i flag of 'stringtie --merge'.",
            category: "advanced"
        }
        label: {
            description: "Equivalent to the -l option of 'stringtie --merge'.",
            category: "advanced"
        }
        memory: {
            description: "The amount of memory needed for this task in GB.",
            category: "advanced"
        }
        dockerImage: {
            description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
            category: "advanced"
        }
    }
}
