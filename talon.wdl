version 1.0

# Copyright (c) 2019 Sequencing Analysis Support Core - Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

task CreateAbundanceFileFromDatabase {
    input {
        File databaseFile
        String annotationVersion
        String genomeBuild
        String outputPrefix

        File? whitelistFile
        File? datasetsFile

        Int cores = 1
        String memory = "4G"
        String dockerImage = "biocontainers/talon:v4.4.1_cv1"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        talon_abundance \
        ~{"--db=" + databaseFile} \
        ~{"-a " + annotationVersion} \
        ~{"-b " + genomeBuild} \
        ~{"--o=" + outputPrefix} \
        ~{"--whitelist=" + whitelistFile} \
        ~{"-d " + datasetsFile}
    }

    output {
        File outputAbundanceFile = outputPrefix + "_talon_abundance.tsv"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        databaseFile: {
            description: "TALON database.",
            category: "required"
        }
        annotationVersion: {
            description: "Which annotation version to use.",
            category: "required"
        }
        genomeBuild: {
            description: "Genome build to use.",
            category: "required"
        }
        outputPrefix: {
            description: "Output directory path + output file prefix.",
            category: "required"
        }
        whitelistFile: {
            description: "Whitelist file of transcripts to include in the output.",
            category: "advanced"
        }
        datasetsFile: {
            description: "A file indicating which datasets should be included.",
            category: "advanced"
        }
        outputAbundanceFile: {
            description: "Abundance for each transcript in the TALON database across datasets.",
            category: "required"
        }
    }
}

task CreateGtfFromDatabase {
    input {
        File databaseFile
        String genomeBuild
        String annotationVersion
        String outputPrefix
        Boolean observedInDataset = false

        File? whitelistFile
        File? datasetFile

        Int cores = 1
        String memory = "4G"
        String dockerImage = "biocontainers/talon:v4.4.1_cv1"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        talon_create_GTF \
        ~{"--db=" + databaseFile} \
        ~{"-b " + genomeBuild} \
        ~{"-a " + annotationVersion} \
        ~{"--o=" + outputPrefix} \
        ~{"--whitelist=" + whitelistFile} \
        ~{true="--observed" false="" observedInDataset} \
        ~{"-d " + datasetFile}
    }

    output {
        File outputGTFfile = outputPrefix + "_talon.gtf"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        databaseFile: {
            description: "TALON database.",
            category: "required"
        }
        genomeBuild: {
            description: "Genome build to use.",
            category: "required"
        }
        annotationVersion: {
            description: "Which annotation version to use.",
            category: "required"
        }
        outputPrefix: {
            description: "Output directory path + output file prefix.",
            category: "required"
        }
        observedInDataset: {
            description: "The output will only include transcripts that were observed at least once.",
            category: "advanced"
        }
        whitelistFile: {
            description: "Whitelist file of transcripts to include in the output.",
            category: "advanced"
        }
        datasetFile: {
            description: "A file indicating which datasets should be included.",
            category: "advanced"
        }
        outputGTFfile: {
            description: "The genes, transcripts, and exons stored a TALON database in GTF format.",
            category: "required"
        }
    }
}

task FilterTalonTranscripts {
    input {
        File databaseFile
        String annotationVersion
        String outputPrefix

        File? pairingsFile

        Int cores = 1
        String memory = "4G"
        String dockerImage = "biocontainers/talon:v4.4.1_cv1"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        talon_filter_transcripts \
        ~{"--db=" + databaseFile} \
        ~{"-a " + annotationVersion} \
        ~{"--o=" + outputPrefix + "_whitelist.csv"} \
        ~{"-p " + pairingsFile}
    }

    output {
        File outputTranscriptWhitelist = outputPrefix + "_whitelist.csv"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        databaseFile: {
            description: "TALON database.",
            category: "required"
        }
        annotationVersion: {
            description: "Which annotation version to use.",
            category: "required"
        }
        outputPrefix: {
            description: "Output directory path + output file prefix.",
            category: "required"
        }
        pairingsFile: {
            description: "A file indicating which datasets should be considered together.",
            category: "advanced"
        }
    }
}

task GetReadAnnotations {
    input {
        File databaseFile
        String genomeBuild
        String outputPrefix

        File? datasetFile

        Int cores = 1
        String memory = "4G"
        String dockerImage = "biocontainers/talon:v4.4.1_cv1"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        talon_fetch_reads \
        ~{"--db " + databaseFile} \
        ~{"--build " + genomeBuild} \
        ~{"--o " + outputPrefix} \
        ~{"--datasets " + datasetFile}
    }

    output {
        File outputAnnotation = outputPrefix + "_talon_read_annot.tsv"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        databaseFile: {
            description: "TALON database.",
            category: "required"
        }
        genomeBuild: {
            description: "Genome build to use.",
            category: "required"
        }
        outputPrefix: {
            description: "Output directory path + output file prefix.",
            category: "required"
        }
        datasetFile: {
            description: "A file indicating which datasets should be included.",
            category: "advanced"
        }
        outputAnnotation: {
            description: "Read-specific annotation information from a TALON database.",
            category: "required"
        }
    }
}

task InitializeTalonDatabase {
    input {
        File GTFfile
        String genomeBuild
        String annotationVersion
        Int minimumLength = 300
        String novelIDprefix = "TALON"
        Int cutoff5p = 500
        Int cutoff3p = 300
        String outputPrefix

        Int cores = 1
        String memory = "10G"
        String dockerImage = "biocontainers/talon:v4.4.1_cv1"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        talon_initialize_database \
        ~{"--f=" + GTFfile} \
        ~{"--g=" + genomeBuild} \
        ~{"--a=" + annotationVersion} \
        ~{"--l=" +  minimumLength} \
        ~{"--idprefix=" + novelIDprefix} \
        ~{"--5p=" + cutoff5p} \
        ~{"--3p=" + cutoff3p} \
        ~{"--o=" + outputPrefix}
    }

    output {
        File outputDatabase = outputPrefix + ".db"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        GTFfile: {
            description: "GTF annotation containing genes, transcripts, and edges.",
            category: "required"
        }
        genomeBuild: {
            description: "Name of genome build that the GTF file is based on (ie hg38).",
            category: "required"
        }
        annotationVersion: {
            description: "Name of supplied annotation (will be used to label data).",
            category: "required"
        }
        minimumLength: { 
            description: "Minimum required transcript length.",
            category: "common"
        }
        novelIDprefix: {
            description: "Prefix for naming novel discoveries in eventual TALON runs.",
            category: "common"
        }
        cutoff5p: { 
            description: "Maximum allowable distance (bp) at the 5' end during annotation.",
            category: "advanced"
        }
        cutoff3p: {
            description: "Maximum allowable distance (bp) at the 3' end during annotation.",
            category: "advanced"
        }
        outputPrefix: {
            description: "Output directory path + output file prefix.",
            category: "required"
        }
        outputDatabase: {
            description: "TALON database.",
            category: "required"
        }
    }
}

task ReformatGtf {
    input {
        File GTFfile

        Int cores = 1
        String memory = "4G"
        String dockerImage = "biocontainers/talon:v4.4.1_cv1"
    }

    command {
        set -e
        talon_reformat_gtf \
        ~{"-gtf " + GTFfile}
    }

    output {
        File outputReformattedGTF = GTFfile
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        GTFfile: {
            description: "GTF annotation containing genes, transcripts, and edges.",
            category: "required"
        }
    }
}

task SummarizeDatasets {
    input {
        File databaseFile
        Boolean setVerbose = false
        String outputPrefix

        File? datasetGroupsCSV

        Int cores = 1
        String memory = "4G"
        String dockerImage = "biocontainers/talon:v4.4.1_cv1"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        talon_summarize \
        ~{"--db " + databaseFile} \
        ~{true="--verbose" false="" setVerbose} \
        ~{"--o " + outputPrefix} \
        ~{"--groups " + datasetGroupsCSV}
    }

    output {
        File outputSummaryFile = outputPrefix + "_talon_summary.tsv"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        databaseFile: {
            description: "TALON database.",
            category: "required"
        }
        setVerbose: {
            description: "Print out the counts in terminal.",
            category: "advanced"
        }
        outputPrefix: {
            description: "Output directory path + output file prefix.",
            category: "required"
        }
        datasetGroupsCSV: {
            description: "File of comma-delimited dataset groups to process together.",
            category: "advanced"
        }
        outputSummaryFile: {
            description: "Tab-delimited file of gene and transcript counts for each dataset.",
            category: "required"
        }
    }
}

task Talon {
    input {
        File configFile
        File databaseFile
        String genomeBuild
        Float minimumCoverage = 0.9
        Int minimumIdentity = 0
        String outputPrefix
        String configFileName = basename(configFile)

        Int cores = 4
        String memory = "25G"
        String dockerImage = "biocontainers/talon:v4.4.1_cv1"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        mv ${configFile} ./${configFileName}
        export TMPDIR=/tmp
        talon \
        ~{"--f " + configFileName} \
        ~{"--db " + databaseFile} \
        ~{"--build " + genomeBuild} \
        ~{"--threads " + cores} \
        ~{"--cov " + minimumCoverage} \
        ~{"--identity " + minimumIdentity} \
        ~{"--o " + outputPrefix + "/run"}
    }

    output {
        File outputUpdatedDatabase = databaseFile
        File outputLog = outputPrefix + "/run_QC.log"
        File outputAnnot = outputPrefix + "/run_talon_read_annot.tsv"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        configFile: {
            description: "Dataset config file (comma-delimited).",
            category: "required"
        }
        databaseFile: {
            description: "TALON database. Created using initialize_talon_database.py.",
            category: "required"
        }
        genomeBuild: {
            description: "Genome build (i.e. hg38) to use.",
            category: "required"
        }
        minimumCoverage: {
            description: "Minimum alignment coverage in order to use a SAM entry.",
            category: "common"
        }
        minimumIdentity: {
            description: "Minimum alignment identity in order to use a SAM entry.",
            category: "common" 
        }
        outputPrefix: {
            description: "Output directory path + output file prefix.",
            category: "required"
        }
        outputUpdatedDatabase: {
            description: "Updated TALON database.",
            category: "required"
        }
        outputLog: {
            description: "Log file from TALON run.",
            category: "required"
        }
        outputAnnot: {
            description: "Read annotation file from TALON run.",
            category: "required"
        }
    }
}
