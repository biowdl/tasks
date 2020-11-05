version 1.0

# Copyright (c) 2019 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
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
        File abundanceFile = outputPrefix + "_talon_abundance.tsv"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        databaseFile: {description: "Talon database.", category: "required"}
        annotationVersion: {description: "Which annotation version to use.", category: "required"}
        genomeBuild: {description: "Genome build to use.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        whitelistFile: {description: "Whitelist file of transcripts to include in the output.", category: "advanced"}
        datasetsFile: {description: "A file indicating which datasets should be included.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        abundanceFile: {description: "Abundance for each transcript in the talon database across datasets."}

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
        File gtfFile = outputPrefix + "_talon.gtf"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        databaseFile: {description: "Talon database.", category: "required"}
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
        gtfFile: {description: "The genes, transcripts, and exons stored a talon database in gtf format."}
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
        File transcriptWhitelist = outputPrefix + "_whitelist.csv"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        databaseFile: {description: "Talon database.", category: "required"}
        annotationVersion: {description: "Which annotation version to use.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        maxFracA: {description: "Maximum fraction of As to allow in the window located immediately after any read assigned to a novel transcript.", category: "advanced"}
        minCount: {description: "Number of minimum occurrences required for a novel transcript per dataset.", category: "advanced"}
        allowGenomic: {description: "If this option is set, transcripts from the Genomic novelty category will be permitted in the output.", category: "advanced"}
        datasetsFile: {description: "Datasets to include.", category: "advanced"}
        minDatasets: {description: "Minimum number of datasets novel transcripts must be found in.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        transcriptWhitelist: {description: "Transcript whitelist produced from the talon database."}
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
        File readAnnotations = outputPrefix + "_talon_read_annot.tsv"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        databaseFile: { description: "Talon database.", category: "required"}
        genomeBuild: {description: "Genome build to use.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        datasetFile: {description: "A file indicating which datasets should be included.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        readAnnotations: {description: "Read-specific annotation information from a talon database."}
    }
}

task GetSpliceJunctions {
    input {
        File sjInformationFile
        String inputFileType = "db"
        File referenceGtf
        String runMode = "intron"
        String outputPrefix

        String memory = "4G"
        Int timeMinutes = 30
        String dockerImage = "biocontainers/talon:v5.0_cv1"
    }

    Map[String, String] SJfileType = {"db": "--db", "gtf": "--gtf"}

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        talon_get_sjs \
        ~{SJfileType[inputFileType] + sjInformationFile} \
        --ref ~{referenceGtf} \
        --mode ~{runMode} \
        --outprefix ~{outputPrefix}
    }

    output {
        File spliceJunctions = outputPrefix + "_" + runMode + "s.tsv"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        sjInformationFile: {description: "Talon gtf file or database from which to extract exons/introns.", category: "required"}
        inputFileType: {description: "The file type of sjInformationFile.", category: "common"}
        referenceGtf: {description: "Gtf reference file (ie gencode).", category: "required"}
        runMode: {description: "Determines whether to include introns or exons in the output.", category: "common"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        spliceJunctions: {description: "File containing locations, novelty and transcript assignments of exons/introns."}
    }
}

task InitializeTalonDatabase {
    input {
        File gtfFile
        String genomeBuild
        String annotationVersion
        Int minimumLength = 300
        String novelPrefix = "TALON"
        Int cutOff5p = 500
        Int cutOff3p = 300
        String outputPrefix

        String memory = "10G"
        Int timeMinutes = 60
        String dockerImage = "biocontainers/talon:v5.0_cv1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        talon_initialize_database \
        --f=~{gtfFile} \
        --g=~{genomeBuild} \
        --a=~{annotationVersion} \
        --l=~{minimumLength} \
        --idprefix=~{novelPrefix} \
        --5p=~{cutOff5p} \
        --3p=~{cutOff3p} \
        --o=~{outputPrefix}
    }

    output {
        File databaseFile = outputPrefix + ".db"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        gtfFile: {description: "Gtf annotation containing genes, transcripts, and edges.", category: "required"}
        genomeBuild: {description: "Name of genome build that the gtf file is based on (ie hg38).", category: "required"}
        annotationVersion: {description: "Name of supplied annotation (will be used to label data).", category: "required"}
        minimumLength: { description: "Minimum required transcript length.", category: "common"}
        novelPrefix: {description: "Prefix for naming novel discoveries in eventual talon runs.", category: "common"}
        cutOff5p: { description: "Maximum allowable distance (bp) at the 5' end during annotation.", category: "advanced"}
        cutOff3p: {description: "Maximum allowable distance (bp) at the 3' end during annotation.", category: "advanced"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        databaseFile: {description: "Talon database."}
    }
}

task LabelReads {
    input {
        File inputSam
        File referenceGenome
        Int fracaRangeSize = 20
        String tmpDir = "./tmp_label_reads"
        Boolean deleteTmp = true
        String outputPrefix

        Int threads = 4
        String memory = "25G"
        Int timeMinutes = 2880
        String dockerImage = "biocontainers/talon:v5.0_cv1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        talon_label_reads \
        --f=~{inputSam} \
        --g=~{referenceGenome} \
        --t=~{threads} \
        --ar=~{fracaRangeSize} \
        --tmpDir=~{tmpDir} \
        ~{true="--deleteTmp" false="" deleteTmp} \
        --o=~{outputPrefix}
    }

    output {
        File labeledSam = outputPrefix + "_labeled.sam"
        File readLabels = outputPrefix + "_read_labels.tsv"
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputSam: {description: "Sam file of transcripts.", category: "required"}
        referenceGenome: {description: "Reference genome fasta file.", category: "required"}
        fracaRangeSize: {description: "Size of post-transcript interval to compute fraction.", category: "common"}
        tmpDir: {description: "Path to directory for tmp files.", category: "advanced"}
        deleteTmp: {description: "If set, tmp dir will be removed.", category: "advanced"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        threads: {description: "The number of threads to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        labeledSam: {description: "Sam file with labeled transcripts."}
        readLabels: {description: "Tabular file with fraction description per read."}
    }
}

task ReformatGtf {
    input {
        File gtfFile

        String memory = "4G"
        Int timeMinutes = 30
        String dockerImage = "biocontainers/talon:v5.0_cv1"
    }

    command {
        set -e
        talon_reformat_gtf \
        -gtf ~{gtfFile}
    }

    output {
        File reformattedGtf = gtfFile
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        gtfFile: {description: "Gtf annotation containing genes, transcripts, and edges.", category: "required"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        reformattedGtf: {description: "Reformatted gtf file."}
    }
}

task SummarizeDatasets {
    input {
        File databaseFile
        Boolean setVerbose = false
        String outputPrefix

        File? datasetGroupsCsv

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
        ~{"--groups " + datasetGroupsCsv}
    }

    output {
        File summaryFile = outputPrefix + "_talon_summary.tsv"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        databaseFile: {description: "Talon database.", category: "required"}
        setVerbose: {description: "Print out the counts in terminal.", category: "advanced"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        datasetGroupsCsv: {description: "File of comma-delimited dataset groups to process together.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        summaryFile: {description: "Tab-delimited file of gene and transcript counts for each dataset."}
    }
}

task Talon {
    input {
        Array[File] samFiles
        String organism
        String sequencingPlatform = "PacBio-RS-II"
        File databaseFile
        String genomeBuild
        Float minimumCoverage = 0.9
        Float minimumIdentity = 0.8
        String outputPrefix

        Int threads = 4
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
        for file in ~{sep=" " samFiles}
        do
            configFileLine="$(basename ${file%.*}),~{organism},~{sequencingPlatform},${file}"
            echo ${configFileLine} >> ~{outputPrefix}/talonConfigFile.csv
        done
        talon \
        ~{"--f " + outputPrefix + "/talonConfigFile.csv"} \
        --db ~{databaseFile} \
        --build ~{genomeBuild} \
        --threads ~{threads} \
        --cov ~{minimumCoverage} \
        --identity ~{minimumIdentity} \
        ~{"--o " + outputPrefix + "/run"}
    >>>

    output {
        File updatedDatabase = databaseFile
        File talonLog = outputPrefix + "/run_QC.log"
        File talonAnnotation = outputPrefix + "/run_talon_read_annot.tsv"
        File talonConfigFile = outputPrefix + "/talonConfigFile.csv"
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        samFiles: {description: "Input sam files.", category: "required"}
        organism: {description: "The name of the organism from which the samples originated.", category: "required"}
        sequencingPlatform: {description: "The sequencing platform used to generate long reads.", category: "required"}
        databaseFile: {description: "Talon database. Created using initialize_talon_database.py.", category: "required"}
        genomeBuild: {description: "Genome build (i.e. hg38) to use.", category: "required"}
        minimumCoverage: {description: "Minimum alignment coverage in order to use a sam entry.", category: "common"}
        minimumIdentity: {description: "Minimum alignment identity in order to use a sam entry.", category: "common" }
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        threads: {description: "The number of threads to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        updatedDatabase: {description: "Updated talon database."}
        talonLog: {description: "Log file from talon run."}
        talonAnnotation: {description: "Read annotation file from talon run."}
        talonConfigFile: {description: "The talon configuration file."}
    }
}
