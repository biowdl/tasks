task ScatterIntervalList {
    String? preCommand
    File interval_list
    Int scatter_count
    String? picardJar

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))

    String toolCommand = if defined(picardJar)
    then "java -Xmx" + mem + "G -jar " + picardJar
    else "picard -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        ${preCommand}
        mkdir scatter_list
        ${toolCommand} \
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
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
    String? preCommand
    Array[File]+ input_bams
    String output_bam_path
    Int? compression_level
    String? picardJar

    Float? memory
    Float? memoryMultiplier

    Int mem = ceil(select_first([memory, 4.0]))

    String toolCommand = if defined(picardJar)
    then "java -Xmx" + mem + "G -jar " + picardJar
    else "picard -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        ${preCommand}
        ${toolCommand} \
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
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
    String? preCommand
    Array[File] input_bams
    String output_bam_path
    String metrics_path
    Int? compression_level
    String? picardJar

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

    String toolCommand = if defined(picardJar)
    then "java -Xmx" + mem + "G -jar " + picardJar
    else "picard -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        ${preCommand}
        mkdir -p $(dirname ${output_bam_path})
        ${toolCommand} \
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
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
    String? preCommand
    Array[File] inputVCFs
    Array[File] inputVCFsIndexes
    String outputVCFpath
    Int? compressionLevel
    String? picardJar

    Float? memory
    Float? memoryMultiplier

    # Using MergeVcfs instead of GatherVcfs so we can create indices
    # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
    Int mem = ceil(select_first([memory, 4.0]))

    String toolCommand = if defined(picardJar)
    then "java -Xmx" + mem + "G -jar " + picardJar
    else "picard -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        ${preCommand}
        ${toolCommand} \
          MergeVcfs \
          INPUT=${sep=' INPUT=' inputVCFs} \
          OUTPUT=${outputVCFpath}
    }

    output {
        File outputVCF = outputVCFpath
        File outputVCFindex = outputVCFpath + ".tbi"
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}

task SamToFastq {
    String? preCommand
    File inputBam
    String outputRead1
    String? outputRead2
    String? outputUnpaired
    String? picardJar
    Float? memory
    Float? memoryMultiplier
    Int mem = ceil(select_first([memory, 16.0])) # High memory default to avoid crashes.

    String toolCommand = if defined(picardJar)
    then "java -Xmx" + mem + "G -jar " + picardJar
    else "picard -Xmx" + mem + "G"

    command {
        set -e -o pipefail
        ${preCommand}
        ${toolCommand} \
        SamToFastq \
        I=${inputBam} \
        ${"FASTQ=" + outputRead1} \
        ${"SECOND_END_FASTQ=" + outputRead2} \
        ${"UNPAIRED_FASTQ=" + outputUnpaired}
    }

    output {
        File read1 = outputRead1
        File? read2 = outputRead2
        File? unpairedRead = outputUnpaired
    }

    runtime {
        memory: ceil(mem * select_first([memoryMultiplier, 3.0]))
    }
}