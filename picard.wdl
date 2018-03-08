task ScatterIntervalList {
    File interval_list
    Int scatter_count
    String picard_jar

    command <<<
        set -e
        mkdir scatter_list
        java -Xmx4G -jar ${picard_jar} \
          IntervalListTools \
          SCATTER_COUNT=${scatter_count} \
          SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
          UNIQUE=true \
          SORT=true \
          INPUT=${interval_list} \
          OUTPUT=scatter_list
    >>>
    output {
        Array[File] out = glob("scatter_list/*/*.interval_list")
        Int interval_count = read_int(stdout())
    }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  Array[File]+ input_bams
  String output_bam_path
  Int? compression_level
  String picard_jar

  command {
    java ${"-Dsamjdk.compression_level=" + compression_level} -Xmx6000 -jar ${picard_jar} \
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
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  Array[File] input_bams
  String output_bam_path
  String metrics_path
  Int? compression_level
  String picard_jar

  # The program default for READ_NAME_REGEX is appropriate in nearly every case.
  # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
  # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
  String? read_name_regex

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    java ${"-Dsamjdk.compression_level=" + compression_level} -Xmx4000 -jar ${picard_jar} \
      MarkDuplicates \
      INPUT=${sep=' INPUT=' input_bams} \
      OUTPUT=${output_bam_path} \
      METRICS_FILE=${metrics_path} \
      VALIDATION_STRINGENCY=SILENT \
      ${"READ_NAME_REGEX=" + read_name_regex} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT="false" \
      CREATE_INDEX=true \
      ADD_PG_TAG_TO_READS=false
  }
  output {
    File output_bam = output_bam_path
    File output_bam_index = sub(output_bam_path, ".bam$", ".bai")
    File duplicate_metrics = metrics_path
  }
}
