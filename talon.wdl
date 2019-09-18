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
        String outputPrefix
        String genomeBuild
        String annotationVersion
        Boolean filterTranscripts = false

        File? filterPairingsFile

        Int cores = 1
        Int memory = 4
        String dockerImage = "biocontainers/talon:v4.2_cv2"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        create_abundance_file_from_database \
        ~{"--db=" + databaseFile} \
        ~{"--o=" + outputPrefix} \
        ~{"-b " + genomeBuild} \
        ~{"-a " + annotationVersion} \
        ~{true="--filter" false="" filterTranscripts} \
        ~{"-p " + filterPairingsFile}
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
        databaseFile: "TALON database."
        outputPrefix: "Output directory path + output file prefix."
        genomeBuild: "Genome build to use."
        annotationVersion: "Which annotation version to use."
        filterTranscripts: "The transcripts in the database will be filtered prior to GTF creation."
        filterPairingsFile: "A file indicating which datasets should be considered together."

        outputAbundanceFile: "Abundance for each transcript in the TALON database across datasets."
    }
}

task CreateGtfAbundanceFromDatabase {
    input {
        File databaseFile
        String outputPrefix
        String genomeBuild
        String annotationVersion
        Boolean filterTranscripts = false

        File? filterPairingsFile

        Int cores = 1
        Int memory = 4
        String dockerImage = "biocontainers/talon:v4.2_cv2"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        create_GTF_abundance_from_database \
        ~{"--db=" + databaseFile} \
        ~{"--o=" + outputPrefix} \
        ~{"-b " + genomeBuild} \
        ~{"-a " + annotationVersion} \
        ~{true="--filter" false="" filterTranscripts} \
        ~{"-p " + filterPairingsFile}
    }

    output {
        File outputGTFfile = outputPrefix + "_talon_observedOnly.gtf"
        File outputAbundanceFile = outputPrefix + "_talon_abundance.tsv"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        databaseFile: "TALON database."
        outputPrefix: "Output directory path + output file prefix."
        genomeBuild: "Genome build to use."
        annotationVersion: "Which annotation version to use."
        filterTranscripts: "The transcripts in the database will be filtered prior to GTF creation."
        filterPairingsFile: "A file indicating which datasets should be considered together."

        outputGTFfile: "The genes, transcripts, and exons stored a TALON database in GTF format."
        outputAbundanceFile: "Abundance for each transcript in the TALON database across datasets."
    }
}

task CreateGtfFromDatabase {
    input {
        File databaseFile
        String outputPrefix
        String genomeBuild
        String annotationVersion
        Boolean observedInDataset = false

        File? whitelistFile
        File? datasetFile

        Int cores = 1
        Int memory = 4
        String dockerImage = "biocontainers/talon:v4.2_cv2"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        create_GTF_from_database \
        ~{"--db=" + databaseFile} \
        ~{"--o=" + outputPrefix} \
        ~{"-b " + genomeBuild} \
        ~{"-a " + annotationVersion} \
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
        databaseFile: "TALON database."
        outputPrefix: "Output directory path + output file prefix."
        genomeBuild: "Genome build to use."
        annotationVersion: "Which annotation version to use."
        observedInDataset: "Output only includes transcripts that were observed at least once."
        whitelistFile: "Whitelist file of transcripts to include in the output."
        datasetFile: "A file indicating which datasets should be included."

        outputGTFfile: "The genes, transcripts, and exons stored a TALON database in GTF format."
    }
}

task InitializeTalonDatabase {
    input {
        File GTFfile
        String outputPrefix
        String genomeBuild
        String annotationVersion
        Int minimumLength = 300
        String novelIDprefix = "TALON"
        Int cutoff5p = 500
        Int cutoff3p = 300

        Int cores = 1
        Int memory = 10
        String dockerImage = "biocontainers/talon:v4.2_cv2"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        initialize_talon_database \
        ~{"--f=" + GTFfile} \
        ~{"--o=" + outputPrefix} \
        ~{"--g=" + genomeBuild} \
        ~{"--a=" + annotationVersion} \
        ~{"--l=" +  minimumLength} \
        ~{"--idprefix=" + novelIDprefix} \
        ~{"--5p=" + cutoff5p} \
        ~{"--3p=" + cutoff3p}
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
        GTFfile: "GTF annotation containing genes, transcripts, and edges."
        outputPrefix: "Output directory path + output file prefix."
        genomeBuild: "Name of genome build that the GTF file is based on (ie hg38)."
        annotationVersion: "Name of supplied annotation (will be used to label data)."
        minimumLength: "Minimum required transcript length."
        novelIDprefix: "Prefix for naming novel discoveries in eventual TALON runs."
        cutoff5p: "Maximum allowable distance (bp) at the 5' end during annotation."
        cutoff3p: "Maximum allowable distance (bp) at the 3' end during annotation."

        outputDatabase: "TALON database."
    }
}

task MapAntisenseGenesToSense {
    input {
        File databaseFile
        String outputPrefix
        String annotationVersion

        Int cores = 1
        Int memory = 4
        String dockerImage = "biocontainers/talon:v4.2_cv2"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        map_antisense_genes_to_sense \
        ~{"--db=" + databaseFile} \
        ~{"--o=" + outputPrefix} \
        ~{"-a " + annotationVersion}
    }

    output {
        File outputAntisenseMapFile = outputPrefix + "_antisense_mapping.gtf"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        databaseFile: "TALON database."
        outputPrefix: "Output directory path + output file prefix."
        annotationVersion: "Which annotation version to use."

        outputAntisenseMapFile: "IDs of the sense gene for every antisense gene in the database."
    }
}

task SummarizeDatasets {
    input {
        File databaseFile
        String outputPrefix

        File? datasetGroupsCSV

        Int cores = 1
        Int memory = 4
        String dockerImage = "biocontainers/talon:v4.2_cv2"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        summarize_datasets \
        ~{"--db " + databaseFile} \
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
        databaseFile: "TALON database."
        outputPrefix: "Output directory path + output file prefix."
        datasetGroupsCSV: "File of comma-delimited dataset groups to process together."

        outputSummaryFile: "Tab-delimited file of gene and transcript counts for each dataset."
    }
}

task Talon {
    input {
        File SAMfile
        File configFile
        File databaseFile
        String outputPrefix
        String genomeBuild
        String configFileName = basename(configFile)
        String SAMfileName = basename(SAMfile)
        Float minimumCoverage = 0.9
        Int minimumIdentity = 0

        Int cores = 1
        Int memory = 20
        String dockerImage = "biocontainers/talon:v4.2_cv2"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPrefix})
        mv ${configFile} ./${configFileName}
        mv ${SAMfile} ./${SAMfileName}
        talon \
        ~{"--f " + configFileName} \
        ~{"--db " + databaseFile} \
        ~{"--o " + outputPrefix} \
        ~{"--build " + genomeBuild} \
        ~{"--cov " + minimumCoverage} \
        ~{"--identity " + minimumIdentity}
    }

    output {
        File outputUpdatedDatabase = databaseFile
        File outputLog = outputPrefix + "_talon_QC.log"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        SAMfile: "Input SAM file, same one as described in configFile."
        configFile: "Dataset config file."
        databaseFile: "TALON database. Created using initialize_talon_database.py."
        outputPrefix: "Output directory path + output file prefix."
        genomeBuild: "Genome build (i.e. hg38) to use."
        minimumCoverage: "Minimum alignment coverage in order to use a SAM entry."
        minimumIdentity: "Minimum alignment identity in order to use a SAM entry."

        outputUpdatedDatabase: "Updated TALON database."
        outputLog: "Log file from TALON run."
    }
}
