task Star {
    String? preCommand

    File inputR1
    File? inputR2
    String genomeDir
    String outFileNamePrefix

    String? outSAMtype = "BAM SortedByCoordinate"
    String? readFilesCommand = "zcat"
    Int? runThreadN
    String? outStd
    String? twopassMode
    String? outSAMattrRGline

    Map[String, String] samOutputNames = {"BAM SortedByCoordinate": "sortedByCoord.out.bam"} #needs to be extended for all possible values

    command {
        set -e -o pipefail
        mkdir -p ${sub(outFileNamePrefix, basename(outFileNamePrefix) + "$", "")}
        ${preCommand}
        STAR \
        --readFilesIn ${inputR1} ${inputR2} \
        --outFileNamePrefix ${outFileNamePrefix} \
        --genomeDir ${genomeDir} \
        ${"--readFilesCommand " + readFilesCommand} \
        ${"--outSAMtype " + outSAMtype} \
        ${"--runThreadN " + runThreadN} \
        ${"--outStd " + outStd} \
        ${"--twopassMode " + twopassMode} \
        ${"--outSAMattrRGline " + outSAMattrRGline}
    }

    output {
        File bamFile = outFileNamePrefix + "Aligned." +  samOutputNames["${outSAMtype}"]
    }
}