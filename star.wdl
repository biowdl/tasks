version 1.0

task Star {
    input {
        Array[File]+ inputR1
        Array[File]? inputR2
        Array[File]+ indexFiles
        String outFileNamePrefix
        String outSAMtype = "BAM SortedByCoordinate"
        String readFilesCommand = "zcat"
        String? outStd
        String? twopassMode = "Basic"
        Array[String]? outSAMattrRGline
        String? outSAMunmapped = "Within KeepPairs"
        Int? limitBAMsortRAM

        Int runThreadN = 4
        String memory = "48G"
        String dockerImage = "quay.io/biocontainers/star:2.6.0c--0"
    }

    #TODO Needs to be extended for all possible output extensions
    Map[String, String] samOutputNames = {"BAM SortedByCoordinate": "sortedByCoord.out.bam"}

    command {
        set -e
        mkdir -p $(dirname ~{outFileNamePrefix})
        STAR \
        --readFilesIn ~{sep=',' inputR1} ~{sep="," inputR2} \
        --outFileNamePrefix ~{outFileNamePrefix} \
        --genomeDir ~{sub(indexFiles[0], basename(indexFiles[0]), "")} \
        --outSAMtype ~{outSAMtype} \
        --readFilesCommand ~{readFilesCommand} \
        ~{"--outSAMunmapped " + outSAMunmapped} \
        ~{"--runThreadN " + runThreadN} \
        ~{"--outStd " + outStd} \
        ~{"--twopassMode " + twopassMode} \
        ~{"--limitBAMsortRAM " + limitBAMsortRAM} \
        ~{true="--outSAMattrRGline " false="" defined(outSAMattrRGline)} ~{sep=" , " outSAMattrRGline}
    }

    output {
        File bamFile = outFileNamePrefix + "Aligned." +  samOutputNames[outSAMtype]
    }

    runtime {
        cpu: runThreadN
        memory: memory
        docker: dockerImage
    }
}

task MakeStarRGline {
    input {
        String sample
        String library
        String platform = "ILLUMINA"
        String readgroup
    }

    command {
        printf '"ID:~{readgroup}" "LB:~{library}" "PL:~{platform}" "SM:~{sample}"'
    }

    output {
        String rgLine = read_string(stdout())
    }
}
