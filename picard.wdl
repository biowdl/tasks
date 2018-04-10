task ScatterIntervalList {
    String? preCommand
    File interval_list
    Int scatter_count
    String picard_jar

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        mkdir scatter_list
        java -Xmx${mem}G -jar ${picard_jar} \
          IntervalListTools \
          SCATTER_COUNT=${scatter_count} \
          SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
          UNIQUE=true \
          SORT=true \
          INPUT=${interval_list} \
          OUTPUT=scatter_list
    }

    output {
        Array[File] out = glob("scatter_list/*/*.interval_list")
        Int interval_count = read_int(stdout())
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 1.5]))
    }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
    String? preCommand
    Array[File]+ input_bams
    String output_bam_path
    Int? compression_level
    String picard_jar

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java ${"-Dsamjdk.compression_level=" + compression_level} \
        -Xmx${mem}G -jar ${picard_jar} \
          GatherBamFiles \
          INPUT=${sep=' INPUT=' input_bams} \
          OUTPUT=${output_bam_path} \
          CREATE_INDEX=true \
          CREATE_MD5_FILE=true
    }

    output {
        File output_bam = "${output_bam_path}"
        File output_bam_index = sub(output_bam_path, ".bam$", ".bai")
        File output_bam_md5 = "${output_bam_path}.md5"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 1.5]))
    }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
    String? preCommand
    Array[File] input_bams
    String output_bam_path
    String metrics_path
    Int? compression_level
    String picard_jar

    Float? memory
    Float? memoryMultiplier

    # The program default for READ_NAME_REGEX is appropriate in nearly every case.
    # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
    # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
    String? read_name_regex

    # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
    # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
    # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        mkdir -p $(dirname ${output_bam_path})
        java ${"-Dsamjdk.compression_level=" + compression_level} \
        -Xmx${mem}G -jar ${picard_jar} \
          MarkDuplicates \
          INPUT=${sep=' INPUT=' input_bams} \
          OUTPUT=${output_bam_path} \
          METRICS_FILE=${metrics_path} \
          VALIDATION_STRINGENCY=SILENT \
          ${"READ_NAME_REGEX=" + read_name_regex} \
          OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
          CLEAR_DT="false" \
          CREATE_INDEX=true \
          ADD_PG_TAG_TO_READS=false
    }

    output {
        File output_bam = output_bam_path
        File output_bam_index = sub(output_bam_path, ".bam$", ".bai")
        File duplicate_metrics = metrics_path
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 1.5]))
    }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
    String? preCommand
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_path
    Int? compression_level
    String picard_jar

    Float? memory
    Float? memoryMultiplier

    # Using MergeVcfs instead of GatherVcfs so we can create indices
    # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
    Int mem = ceil(select_first([memory, 4.0]))
    command {
        set -e -o pipefail
        ${preCommand}
        java ${"-Dsamjdk.compression_level=" + compression_level} \
        -Xmx${mem}G -jar ${picard_jar} \
          MergeVcfs \
          INPUT=${sep=' INPUT=' input_vcfs} \
          OUTPUT=${output_vcf_path}
    }

    output {
        File output_vcf = output_vcf_path
        File output_vcf_index = output_vcf_path + ".tbi"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 1.5]))
    }
}