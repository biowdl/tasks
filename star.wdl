task Star {
    String? preCommand

    Array[File] inputR1
    Array[File]? inputR2
    String genomeDir
    String outFileNamePrefix

    String? outSAMtype = "BAM SortedByCoordinate"
    String? readFilesCommand = "zcat"
    Int? runThreadN
    String? outStd
    String? twopassMode
    Array[String]? outSAMattrRGline

    #TODO needs to be extended for all possible output extensions
    Map[String, String] samOutputNames = {"BAM SortedByCoordinate": "sortedByCoord.out.bam"}

    command {
        set -e -o pipefail
        mkdir -p ${sub(outFileNamePrefix, basename(outFileNamePrefix) + "$", "")}
        ${preCommand}
        STAR \
        --readFilesIn ${sep=',' inputR1} ${sep="," inputR2} \
        --outFileNamePrefix ${outFileNamePrefix} \
        --genomeDir ${genomeDir} \
        ${"--readFilesCommand " + readFilesCommand} \
        ${"--outSAMtype " + outSAMtype} \
        ${"--runThreadN " + runThreadN} \
        ${"--outStd " + outStd} \
        ${"--twopassMode " + twopassMode} \
        ${true="--outSAMattrRGline " false="" defined(outSAMattrRGline)} ${sep=" , " outSAMattrRGline}
    }

    output {
        File bamFile = outFileNamePrefix + "Aligned." +  samOutputNames["${outSAMtype}"]
    }

    runtime {
        threads: runThreadN
    }
}