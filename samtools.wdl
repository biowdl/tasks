task getContigMostCoverage {
    String? preCommand
    File bamFile
    String? outputFilePath = "mostCoveredContig.bam"
    command <<<
        set -e -o pipefail
        ${preCommand}
        samtools view -hb wgs1.bam $( \
            samtools idxstats wgs1.bam \
            | grep -ve "^\*" \
            | awk '{ print $1"\t"($3/$2)}' \
            | sort -rnk2 \
            | head -n 1 \
            | cut -f1) \
        > ${outputFilePath}
        samtools index ${outputFilePath}
    >>>

    output {
        File contig = outputFilePath
    }

    parameter_meta {
        preCommand: "This command is run before running the samtools command. Can be used to set up environments such as conda."
        bamFile: "Must be a sorted BAM file with an index"
        outputFilePath: "where the output file is stored"
    }
}