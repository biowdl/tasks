version 1.0

task Fastqc {
    input {
        File seqFile
        String outdirPath
        String? preCommand
        Boolean? casava
        Boolean? nano
        Boolean? noFilter
        Boolean extract = true
        Boolean? nogroup
        Int? minLength
        String? format
        Int threads = 1
        File? contaminants
        File? adapters
        File? limits
        Int? kmers
        String? dir
    }

    # Chops of the .gz extension if present.
    String name = sub(seqFile, "\\.gz$","")
    # This regex chops of the extension and replaces it with _fastqc for the reportdir.
    # Just as fastqc does it.
    String reportDir = outdirPath + "/" + sub(basename(name), "\\.[^\\.]*$", "_fastqc")

    command {
        set -e -o pipefail
        ~{preCommand}
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
    }
}

task GetConfiguration {
    input {
        String? preCommand
        String fastqcDirFile = "fastqcDir.txt"
    }

    command {
        set -e -o pipefail
        ~{preCommand}
        echo $(dirname $(readlink -f $(which fastqc))) > ~{fastqcDirFile}
    }

    output {
        String fastqcDir = read_string(fastqcDirFile)
        File adapterList = fastqcDir + "/Configuration/adapter_list.txt"
        File contaminantList = fastqcDir + "/Configuration/contaminant_list.txt"
        File limits = fastqcDir + "/Configuration/limits.txt"
    }

    runtime {
        memory: 1
    }
}
