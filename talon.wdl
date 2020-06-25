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

        String memory = "4G"
        Int timeMinutes = 30
        String dockerImage = "biocontainers/talon:v5.0_cv1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        talon_abundance \
        --db=~{databaseFile} \
        -a ~{annotationVersion} \
        -b ~{genomeBuild} \
        --o=~{outputPrefix} \
        ~{"--whitelist=" + whitelistFile} \
        ~{"-d " + datasetsFile}
    }

    output {
        File outputAbundanceFile = outputPrefix + "_talon_abundance.tsv"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        databaseFile: {description: "TALON database.", category: "required"}
        annotationVersion: {description: "Which annotation version to use.", category: "required"}
        genomeBuild: {description: "Genome build to use.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        whitelistFile: {description: "Whitelist file of transcripts to include in the output.", category: "advanced"}
        datasetsFile: {description: "A file indicating which datasets should be included.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputAbundanceFile: {description: "Abundance for each transcript in the TALON database across datasets."}

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

        String memory = "4G"
        Int timeMinutes = 30
        String dockerImage = "biocontainers/talon:v5.0_cv1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        talon_create_GTF \
        --db=~{databaseFile} \
        -b ~{genomeBuild} \
        -a ~{annotationVersion} \
        --o=~{outputPrefix} \
        ~{true="--observed" false="" observedInDataset} \
        ~{"--whitelist=" + whitelistFile} \
        ~{"-d " + datasetFile}
    }

    output {
        File outputGTFfile = outputPrefix + "_talon.gtf"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        databaseFile: {description: "TALON database.", category: "required"}
        genomeBuild: {description: "Genome build to use.", category: "required"}
        annotationVersion: {description: "Which annotation version to use.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        observedInDataset: {description: "The output will only include transcripts that were observed at least once.", category: "advanced"}
        whitelistFile: {description: "Whitelist file of transcripts to include in the output.", category: "advanced"}
        datasetFile: {description: "A file indicating which datasets should be included.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputGTFfile: {description: "The genes, transcripts, and exons stored a TALON database in GTF format."}
    }
}

task FilterTalonTranscripts {
    input {
        File databaseFile
        String annotationVersion
        String outputPrefix
        Float maxFracA = 0.5
        Int minCount = 5
        Boolean allowGenomic = false

        File? datasetsFile
        Int? minDatasets

        String memory = "4G"
        Int timeMinutes = 30
        String dockerImage = "biocontainers/talon:v5.0_cv1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        talon_filter_transcripts \
        --db=~{databaseFile} \
        -a ~{annotationVersion} \
        ~{"--o=" + outputPrefix + "_whitelist.csv"} \
        --maxFracA=~{maxFracA} \
        --minCount=~{minCount} \
        ~{true="--allowGenomic" false="" allowGenomic} \
        --datasets=~{datasetsFile} \
        --minDatasets=~{minDatasets}
    }

    output {
        File outputTranscriptWhitelist = outputPrefix + "_whitelist.csv"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        databaseFile: {description: "TALON database.", category: "required"}
        annotationVersion: {description: "Which annotation version to use.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        maxFracA: {description: "Maximum fraction of As to allow in the window located immediately after any read assigned to a novel transcript.", category: "advanced"}
        minCount: {description: "Number of minimum occurrences required for a novel transcript PER dataset.", category: "advanced"}
        allowGenomic: {description: "If this option is set, transcripts from the Genomic novelty category will be permitted in the output.", category: "advanced"}
        datasetsFile: {description: "Datasets to include.", category: "advanced"}
        minDatasets: {description: "Minimum number of datasets novel transcripts must be found in.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputTranscriptWhitelist: {description: "A transcript whitelist produced from the TALON database."}
    }
}

task GetReadAnnotations {
    input {
        File databaseFile
        String genomeBuild
        String outputPrefix

        File? datasetFile

        String memory = "4G"
        Int timeMinutes = 30
        String dockerImage = "biocontainers/talon:v5.0_cv1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        talon_fetch_reads \
        --db ~{databaseFile} \
        --build ~{genomeBuild} \
        --o ~{outputPrefix} \
        ~{"--datasets " + datasetFile}
    }

    output {
        File outputAnnotation = outputPrefix + "_talon_read_annot.tsv"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        databaseFile: { description: "TALON database.", category: "required"}
        genomeBuild: {description: "Genome build to use.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        datasetFile: {description: "A file indicating which datasets should be included.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputAnnotation: {description: "Read-specific annotation information from a TALON database."}
    }
}

task GetSpliceJunctions {
    input {
        File GTFfile
        File databaseFile
        File referenceGTF
        String runMode = "intron"
        String outputPrefix

        String memory = "4G"
        Int timeMinutes = 30
        String dockerImage = "biocontainers/talon:v5.0_cv1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        talon_get_sjs \
        --gtf ~{GTFfile} \
        --db ~{databaseFile} \
        --ref ~{referenceGTF} \
        --mode ~{runMode} \
        --outprefix ~{outputPrefix}
    }

    output {
        
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        GTFfile: {description: "TALON GTF file from which to extract exons/introns.", category: "required"}
        databaseFile: { description: "TALON database.", category: "required"}
        referenceGTF: {description: "GTF reference file (ie GENCODE).", category: "required"}
        runMode: {description: "Determines whether to include introns or exons in the output.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
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

        String memory = "10G"
        Int timeMinutes = 60
        String dockerImage = "biocontainers/talon:v5.0_cv1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        talon_initialize_database \
        --f=~{GTFfile} \
        --g=~{genomeBuild} \
        --a=~{annotationVersion} \
        --l=~{minimumLength} \
        --idprefix=~{novelIDprefix} \
        --5p=~{cutoff5p} \
        --3p=~{cutoff3p} \
        --o=~{outputPrefix}
    }

    output {
        File outputDatabase = outputPrefix + ".db"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        GTFfile: {description: "GTF annotation containing genes, transcripts, and edges.", category: "required"}
        genomeBuild: {description: "Name of genome build that the GTF file is based on (ie hg38).", category: "required"}
        annotationVersion: {description: "Name of supplied annotation (will be used to label data).", category: "required"}
        minimumLength: { description: "Minimum required transcript length.", category: "common"}
        novelIDprefix: {description: "Prefix for naming novel discoveries in eventual TALON runs.", category: "common"}
        cutoff5p: { description: "Maximum allowable distance (bp) at the 5' end during annotation.", category: "advanced"}
        cutoff3p: {description: "Maximum allowable distance (bp) at the 3' end during annotation.", category: "advanced"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputDatabase: {description: "TALON database."}
    }
}

task ReformatGtf {
    input {
        File GTFfile

        String memory = "4G"
        Int timeMinutes = 30
        String dockerImage = "biocontainers/talon:v5.0_cv1"
    }

    command {
        set -e
        talon_reformat_gtf \
        -gtf ~{GTFfile}
    }

    output {
        File outputReformattedGTF = GTFfile
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        GTFfile: {description: "GTF annotation containing genes, transcripts, and edges.", category: "required"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputReformattedGTF: {description: "Reformatted GTF file."}
    }
}

task SummarizeDatasets {
    input {
        File databaseFile
        Boolean setVerbose = false
        String outputPrefix

        File? datasetGroupsCSV

        String memory = "4G"
        Int timeMinutes = 50
        String dockerImage = "biocontainers/talon:v5.0_cv1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        talon_summarize \
        --db ~{databaseFile} \
        ~{true="--verbose" false="" setVerbose} \
        --o ~{outputPrefix} \
        ~{"--groups " + datasetGroupsCSV}
    }

    output {
        File outputSummaryFile = outputPrefix + "_talon_summary.tsv"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        databaseFile: {description: "TALON database.", category: "required"}
        setVerbose: {description: "Print out the counts in terminal.", category: "advanced"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        datasetGroupsCSV: {description: "File of comma-delimited dataset groups to process together.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputSummaryFile: {description: "Tab-delimited file of gene and transcript counts for each dataset."}
    }
}

task Talon {
    input {
        Array[File] SAMfiles
        String organism
        String sequencingPlatform = "PacBio-RS-II"
        File databaseFile
        String genomeBuild
        Float minimumCoverage = 0.9
        Float minimumIdentity = 0.8
        String outputPrefix

        Int cores = 4
        String memory = "25G"
        Int timeMinutes = 2880
        String dockerImage = "biocontainers/talon:v5.0_cv1"
    }

    command <<<
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        mkdir -p $PWD/tmp #Standard /tmp fills up which makes the SQLite process crash.
        ln -s $PWD/tmp /tmp/sqltmp #Multiprocessing will crash if the absolute path is too long.
        export TMPDIR=/tmp/sqltmp
        printf "" > ~{outputPrefix}/talonConfigFile.csv #File needs to be emptied when task is rerun.
        for file in ~{sep=" " SAMfiles}
        do
            configFileLine="$(basename ${file%.*}),~{organism},~{sequencingPlatform},${file}"
            echo ${configFileLine} >> ~{outputPrefix}/talonConfigFile.csv
        done
        talon \
        ~{"--f " + outputPrefix + "/talonConfigFile.csv"} \
        --db ~{databaseFile} \
        --build ~{genomeBuild} \
        --threads ~{cores} \
        --cov ~{minimumCoverage} \
        --identity ~{minimumIdentity} \
        ~{"--o " + outputPrefix + "/run"}
    >>>

    output {
        File outputUpdatedDatabase = databaseFile
        File outputLog = outputPrefix + "/run_QC.log"
        File outputAnnot = outputPrefix + "/run_talon_read_annot.tsv"
        File outputConfigFile = outputPrefix + "/talonConfigFile.csv"
    }

    runtime {
        cpu: cores
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        SAMfiles: {description: "Input SAM files.", category: "required"}
        organism: {description: "The name of the organism from which the samples originated.", category: "required"}
        sequencingPlatform: {description: "The sequencing platform used to generate long reads.", category: "required"}
        databaseFile: {description: "TALON database. Created using initialize_talon_database.py.", category: "required"}
        genomeBuild: {description: "Genome build (i.e. hg38) to use.", category: "required"}
        minimumCoverage: {description: "Minimum alignment coverage in order to use a SAM entry.", category: "common"}
        minimumIdentity: {description: "Minimum alignment identity in order to use a SAM entry.", category: "common" }
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        cores: {description: "The number of cores to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputUpdatedDatabase: {description: "Updated TALON database."}
        outputLog: {description: "Log file from TALON run."}
        outputAnnot: {description: "Read annotation file from TALON run."}
        outputConfigFile: {description: "The TALON configuration file."}
    }
}
