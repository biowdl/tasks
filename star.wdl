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
        Int runThreadN = 4
        String? outStd
        String? twopassMode
        Array[String]? outSAMattrRGline
        String? outSAMunmapped = "Within KeepPairs"
        Int? limitBAMsortRAM


        Int memory = 48

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
        # Return memory per CPU here due to SGE backend.
        # Can also work with slurms mem-per-cpu flag
        memory: (memory / runThreadN) + 1
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
