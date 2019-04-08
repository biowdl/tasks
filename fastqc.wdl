version 1.0

task Fastqc {
    input {
        File seqFile
        String outdirPath
        Boolean? casava
        Boolean? nano
        Boolean? noFilter
        Boolean extract = true
        Boolean? nogroup
        Int? minLength
        String? format
        File? contaminants
        File? adapters
        File? limits
        Int? kmers
        String? dir

        Int threads = 1
        String dockerTag = "0.11.7--4"
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
        File rawReport = reportDir + "/fastqc_data.txt"
        File htmlReport = reportDir + "/fastqc_report.html"
        File summary = reportDir + "/summary.txt"
        Array[File] images = glob(reportDir + "/Images/*.png")
    }

    runtime {
        cpu: threads
        docker: "quay.io/biocontainers/fastqc:" + dockerTag
    }
}

task GetConfiguration {
    input {
        String dockerTag = "0.11.7--4"
    }

    command <<<
        set -e
        fastqcDir=$(dirname $(readlink -f $(which fastqc)))
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
        memory: 2 # Needs more than 1 to pull the docker image
        docker: "quay.io/biocontainers/fastqc:" + dockerTag
    }
}
