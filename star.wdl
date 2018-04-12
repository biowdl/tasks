task Star {
    String? preCommand

    Array[File] inputR1
    Array[File?] inputR2
    String genomeDir
    String outFileNamePrefix

    String? outSAMtype
    String? readFilesCommand
    Int? runThreadN
    String? outStd
    String? twopassMode
    Array[String]? outSAMattrRGline

    Int? memory

    #TODO needs to be extended for all possible output extensions
    Map[String, String] samOutputNames = {"BAM SortedByCoordinate": "sortedByCoord.out.bam"}

    # converts String? to String for use as key (for the Map above) in output
    String key = select_first([outSAMtype, "BAM SortedByCoordinate"])

    command {
        set -e -o pipefail
        mkdir -p ${sub(outFileNamePrefix, basename(outFileNamePrefix) + "$", "")}
        ${preCommand}
        STAR \
        --readFilesIn ${sep=',' inputR1} ${sep="," inputR2} \
        --outFileNamePrefix ${outFileNamePrefix} \
        --genomeDir ${genomeDir} \
        --outSAMtype ${default="BAM SortedByCoordinate" outSAMtype} \
        --readFilesCommand ${default="zcat" readFilesCommand} \
        ${"--runThreadN " + runThreadN} \
        ${"--outStd " + outStd} \
        ${"--twopassMode " + twopassMode} \
        ${true="--outSAMattrRGline " false="" defined(outSAMattrRGline)} ${sep=" , " outSAMattrRGline}
    }

    output {
        File bamFile = outFileNamePrefix + "Aligned." +  samOutputNames[key]
    }

    runtime {
        cpu: select_first([runThreadN, 1])
        memory: select_first([memory, 10])
    }
}