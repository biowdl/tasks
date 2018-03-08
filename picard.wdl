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
