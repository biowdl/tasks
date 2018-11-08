version 1.0

task Star {
    input {
        String? preCommand

        Array[File] inputR1
        Array[File]? inputR2
        String genomeDir
        String outFileNamePrefix

        String outSAMtype = "BAM SortedByCoordinate"
        String readFilesCommand = "zcat"
        Int runThreadN = 1
        String? outStd
        String? twopassMode
        Array[String]? outSAMattrRGline
        Int? limitBAMsortRAM

        Int memory = 10
    }

    # Needs to be extended for all possible output extensions
    Map[String, String] samOutputNames = {"BAM SortedByCoordinate": "sortedByCoord.out.bam"}

    command {
        set -e -o pipefail
        mkdir -p ~{sub(outFileNamePrefix, basename(outFileNamePrefix) + "$", "")}
        ~{preCommand}
        STAR \
        --readFilesIn ~{sep=',' inputR1} ~{sep="," inputR2} \
        --outFileNamePrefix ~{outFileNamePrefix} \
        --genomeDir ~{genomeDir} \
        --outSAMtype ~{outSAMtype} \
        --readFilesCommand ~{readFilesCommand} \
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
