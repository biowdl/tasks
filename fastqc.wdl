version 1.0

task Fastqc {
    input {
        File seqFile
        String outdirPath
        Boolean casava = false
        Boolean nano = false
        Boolean noFilter = false
        Boolean extract = false
        Boolean nogroup = false
        Int? minLength
        String? format
        File? contaminants
        File? adapters
        File? limits
        Int? kmers
        String? dir

        Int threads = 1
        String memory = "4G"
        String dockerImage = "quay.io/biocontainers/fastqc:0.11.7--4"
        Array[File]? NoneArray
        File? NoneFile
    }

    # Chops of the .gz extension if present.
    # The Basename needs to be taken here. Otherwise paths might differ between similar jobs.
    String name = basename(sub(seqFile, "\.gz$",""))
    # This regex chops of the extension and replaces it with _fastqc for the reportdir.
    # Just as fastqc does it.
    String reportDir = outdirPath + "/" + sub(name, "\.[^\.]*$", "_fastqc")

    command {
        set -e
        mkdir -p ~{outdirPath}
        fastqc \
        ~{"--outdir " + outdirPath} \
        ~{true="--casava" false="" casava} \
        ~{true="--nano" false="" nano} \
        ~{true="--nofilter" false="" noFilter} \
        ~{true="--extract" false="" extract} \
        ~{true="--nogroup" false="" nogroup} \
        ~{"--min_length " + minLength } \
        ~{"--format " + format} \
        ~{"--threads " + threads} \
        ~{"--contaminants " + contaminants} \
        ~{"--adapters " + adapters} \
        ~{"--limits " + limits} \
        ~{"--kmers " + kmers} \
        ~{"--dir " + dir} \
        ~{seqFile}
    }

    output {
        File? rawReport = if extract then reportDir + "/fastqc_data.txt" else NoneFile
        File htmlReport = reportDir + ".html"
        File reportZip = reportDir + ".zip"
        File? summary = if extract then reportDir + "/summary.txt" else NoneFile
        Array[File]? images = if extract then glob(reportDir + "/Images/*.png") else NoneArray
    }

    runtime {
        cpu: threads
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        seqFile: {description: "A fastq file.", category: "required"}
        outdirPath: {description: "The path to write the output to", catgory: "required"}
        casava: {description: "Equivalent to fastqc's --casava flag.", category: "advanced"}
        nano: {description: "Equivalent to fastqc's --nano flag.", category: "advanced"}
        noFilter: {description: "Equivalent to fastqc's --nofilter flag.", category: "advanced"}
        extract: {description: "Equivalent to fastqc's --extract flag.", category: "advanced"}
        nogroup: {description: "Equivalent to fastqc's --nogroup flag.", category: "advanced"}
        minLength: {description: "Equivalent to fastqc's --min_length option.", category: "advanced"}
        format: {description: "Equivalent to fastqc's --format option.", category: "advanced"}
        contaminants: {description: "Equivalent to fastqc's --contaminants option.", category: "advanced"}
        adapters: {description: "Equivalent to fastqc's --adapters option.", category: "advanced"}
        limits: {description: "Equivalent to fastqc's --limits option.", category: "advanced"}
        kmers: {description: "Equivalent to fastqc's --kmers option.", category: "advanced"}
        dir: {description: "Equivalent to fastqc's --dir option.", category: "advanced"}
        threads: {description: "The number of cores to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }

    meta {
        WDL_AID: {
            exclude: ["NoneFile", "NoneArray"]
        }
    }
}

task GetConfiguration {
    input {
        String dockerImage = "quay.io/biocontainers/fastqc:0.11.7--4"
    }

    command <<<
        set -e
        fastqcExe="$(command -v fastqc)"
        fastqcPath="$(readlink -f $fastqcExe)"
        fastqcDir="$(dirname $fastqcPath)"
        mkdir Configuration
        cp ${fastqcDir}/Configuration/adapter_list.txt Configuration/adapter_list.txt
        cp ${fastqcDir}/Configuration/contaminant_list.txt Configuration/contaminant_list.txt
        cp ${fastqcDir}/Configuration/limits.txt Configuration/limits.txt
    >>>

    output {
        File adapterList = "Configuration/adapter_list.txt"
        File contaminantList = "Configuration/contaminant_list.txt"
        File limits = "Configuration/limits.txt"
    }

    runtime {
        memory: "2G" # Needs more than 1 to pull the docker image
        docker: dockerImage
    }

    parameter_meta {
        dockerImage: {
            description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
            category: "advanced"
        }
    }
}
