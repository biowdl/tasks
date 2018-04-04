task fastqc {
    File seqFile
    String outdirPath
    String? preCommand
    Boolean? casava
    Boolean? nano
    Boolean? noFilter
    Boolean? extract = true
    Boolean? nogroup
    Int? minLength
    String? format
    Int? threads = 1
    File? contaminants
    File? adapters
    File? limits
    Int? kmers
    String? dir

    command {
    set -e -o pipefail
    ${preCommand}
    mkdir -p ${outdirPath}
    fastqc \
    ${"--outdir " + outdirPath} \
    ${true="--casava" false="" casava} \
    ${true="--nano" false="" nano} \
    ${true="--nofilter" false="" noFilter} \
    ${true="--extract" false="" extract} \
    ${true="--nogroup" false="" nogroup} \
    ${"--min_length " + minLength } \
    ${"--format " + format} \
    ${"--threads " + threads} \
    ${"--contaminants " + contaminants} \
    ${"--adapters " + adapters} \
    ${"--limits " + limits} \
    ${"--kmers " + kmers} \
    ${"--dir " + dir} \
    ${seqFile}

    }

    output {
        File rawReport = select_first(glob(outdirPath + "/*/fastqc_data.txt"))
        File htmlReport = select_first(glob(outdirPath + "/*/fastqc_report.html"))
        File summary = select_first(glob(outdirPath + "/*/summary.txt"))
        Array[File] images = glob(outdirPath + "/*/Images/*.png")
    }

    runtime {
        cpu: select_first([threads])
    }
}

task extractAdapters {
    File extractAdaptersFastqcJar
    File inputFile
    String outputDir
    String? adapterOutputFilePath = outputDir + "/adapter.list"
    String? contamsOutputFilePath = outputDir + "/contaminations.list"
    Boolean? skipContams
    File? knownContamFile
    File? knownAdapterFile
    Float? adapterCutoff
    Boolean? outputAsFasta
    command {
    set -e
    mkdir -p ${outputDir}
    java -Xmx4G -jar ${extractAdaptersFastqcJar} \
    --inputFile ${inputFile} \
    ${"--adapterOutputFile " + adapterOutputFilePath } \
    ${"--contamsOutputFile " + contamsOutputFilePath } \
    ${"--knownContamFile " + knownContamFile} \
    ${"--knownAdapterFile " + knownAdapterFile} \
    ${"--adapterCutoff " + adapterCutoff} \
    ${true="--skipContams" false="" skipContams} \
    ${true="--outputAsFasta" false="" outputAsFasta}
    }

    output {
        File adapterOutputFile = select_first([adapterOutputFilePath])
        File contamsOutputFile = select_first([contamsOutputFilePath])
        Array[String] adapterList = read_lines(select_first([adapterOutputFilePath]))
        Array[String] contamsList = read_lines(select_first([contamsOutputFilePath]))
    }

    runtime {
        memory: 6
    }
}

task getConfiguration {
    String? preCommand
    String? fastqcDirFile = "fastqcDir.txt"

    command {
        set -e -o pipefail
        ${preCommand}
        echo $(dirname $(readlink -f $(which fastqc))) > ${fastqcDirFile}
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