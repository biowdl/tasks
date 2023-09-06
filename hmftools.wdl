version 1.0

# Copyright (c) 2020 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

task Amber {
    input {
        String? referenceName
        File? referenceBam
        File? referenceBamIndex
        String tumorName
        File tumorBam
        File tumorBamIndex
        String outputDir = "./amber"
        File loci
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String refGenomeVersion

        Int? tumorOnlyMinDepth

        Int threads = 2
        String memory = "85GiB"
        String javaXmx = "80G"
        Int timeMinutes = 480
        String dockerImage = "quay.io/biocontainers/hmftools-amber:3.9--hdfd78af_1"
    }

    command {
        AMBER -Xmx~{javaXmx} \
        ~{"-reference " + referenceName} \
        ~{"-reference_bam " + referenceBam} \
        -tumor ~{tumorName} \
        -tumor_bam ~{tumorBam} \
        -output_dir ~{outputDir} \
        -threads ~{threads} \
        -ref_genome ~{referenceFasta} \
        -ref_genome_version ~{refGenomeVersion} \
        -loci ~{loci} \
        ~{"-tumor-only-min-depth " + tumorOnlyMinDepth}
    }

    output {
        File version = "~{outputDir}/amber.version"
        File tumorBafPcf = "~{outputDir}/~{tumorName}.amber.baf.pcf"
        File tumorBafTsv = "~{outputDir}/~{tumorName}.amber.baf.tsv.gz"
        File tumorContaminationVcf = "~{outputDir}/~{tumorName}.amber.contamination.vcf.gz"
        File tumorContaminationVcfIndex = "~{outputDir}/~{tumorName}.amber.contamination.vcf.gz.tbi"
        File tumorContaminationTsv = "~{outputDir}/~{tumorName}.amber.contamination.tsv"
        File tumorQc = "~{outputDir}/~{tumorName}.amber.qc"
        File normalHomozygousregionsTsv = "~{outputDir}/~{referenceName}.amber.homozygousregion.tsv"
        File normalSnpVcf = "~{outputDir}/~{referenceName}.amber.snp.vcf.gz"
        File normalSnpVcfIndex = "~{outputDir}/~{referenceName}.amber.snp.vcf.gz.tbi"
        Array[File] outputs = [version, tumorBafPcf, tumorBafTsv, tumorContaminationVcf,
            tumorContaminationVcfIndex, tumorContaminationTsv, tumorQc, normalHomozygousregionsTsv,
            normalSnpVcf, normalSnpVcfIndex]
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        cpu: threads
    }

    parameter_meta {
        referenceName: {description: "the name of the normal sample.", category: "required"}
        referenceBam: {description: "The normal BAM file.", category: "required"}
        referenceBamIndex: {description: "The index for the normal BAM file.", category: "required"}
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        tumorBam: {description: "The tumor BAM file.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor BAM file.", category: "required"}
        outputDir: {description: "The path to the output directory.", category: "common"}
        loci: {description: "A VCF file containing likely heterozygous sites.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        refGenomeVersion: {description: "The version of the reference genome: 37 or 38.", category: "required"}
        threads: {description: "The number of threads the program will use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Cobalt {
    input {
        String? referenceName
        File? referenceBam
        File? referenceBamIndex
        String tumorName
        File tumorBam
        File tumorBamIndex
        String outputDir = "./cobalt"
        File gcProfile
        File refGenomeFile

        File? tumorOnlyDiploidBed
        File? targetRegionsNormalisationTsv
        Int? pcfGamma

        Int threads = 1
        String memory = "5GiB"
        String javaXmx = "4G"
        Int timeMinutes = 960
        String dockerImage = "quay.io/biocontainers/hmftools-cobalt:1.13--hdfd78af_1"
    }

    command {
        COBALT -Xmx~{javaXmx} \
        ~{"-reference " + referenceName} \
        ~{"-reference_bam " + referenceBam} \
        -tumor ~{tumorName} \
        -tumor_bam ~{tumorBam} \
        -output_dir ~{outputDir} \
        -threads ~{threads} \
        -gc_profile ~{gcProfile} \
        -ref_genome ~{refGenomeFile} \
        ~{"-tumor_only_diploid_bed " + tumorOnlyDiploidBed} \
        ~{"-target_region " + targetRegionsNormalisationTsv} \
        ~{"-pcf_gamma" + pcfGamma}
    }

    output {
        File version = "~{outputDir}/cobalt.version"
        File normalGcMedianTsv = "~{outputDir}/~{referenceName}.cobalt.gc.median.tsv"
        File normalRationMedianTsv = "~{outputDir}/~{referenceName}.cobalt.ratio.median.tsv"
        File normalRationPcf = "~{outputDir}/~{referenceName}.cobalt.ratio.pcf"
        File tumorGcMedianTsv = "~{outputDir}/~{tumorName}.cobalt.gc.median.tsv"
        File tumorRatioPcf = "~{outputDir}/~{tumorName}.cobalt.ratio.pcf"
        File tumorRatioTsv = "~{outputDir}/~{tumorName}.cobalt.ratio.tsv.gz"
        Array[File] outputs = [version, normalGcMedianTsv, normalRationMedianTsv,
            normalRationPcf, tumorGcMedianTsv, tumorRatioPcf, tumorRatioTsv]
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        cpu: threads
    }

    parameter_meta {
        referenceName: {description: "the name of the normal sample.", category: "required"}
        referenceBam: {description: "The normal BAM file.", category: "required"}
        referenceBamIndex: {description: "The index for the normal BAM file.", category: "required"}
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        tumorBam: {description: "The tumor BAM file.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor BAM file.", category: "required"}
        outputDir: {description: "The path to the output directory.", category: "common"}
        gcProfile: {description: "A file describing the GC profile of the reference genome.", category: "required"}
        refGenomeFile: {description: "The reference genome fasta file.", category: "required"}
        threads: {description: "The number of threads the program will use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task CupGenerateReport {
    input {
        String sampleName
        File cupData
        String outputDir = "./cuppa"

        String memory = "5GiB"
        Int timeMinutes = 10
        String dockerImage = "quay.io/biowdl/cuppa@sha256:e76d367a3226068967fb64ad6adaa889cbdcc01397075b0cbc382bbba4350b98"
    }

    # This script writes to the directory that the input is located in.
    # Giving the input directly will cause the script to write in the
    # locallized input dir, which may cause issues with write permissions
    # in certain execution engines or backends. We, therefore, make links
    # to a working directory, and give that directory as input instead.
    # We can't just use the outputDir directly. This could be an
    # absolute path in which case the linking might fail due to name
    # collisions. Outputs are copied to the given output dir afterwards.
    command {
        set -e
        mkdir -p ./workdir ~{outputDir}
        ln -s -t workdir ~{cupData}
        CupGenerateReport \
        ~{sampleName} \
        workdir/
        mv -t ~{outputDir} \
        ./workdir/~{sampleName}.cup.report.summary.png \
        ./workdir/~{sampleName}_cup_report.pdf
        if [ -f ./workdir/~{sampleName}.cup.report.features.png ]
          then
            mv -t ~{outputDir} \
            ./workdir/~{sampleName}.cup.report.features.png
        fi
    }

    output {
        File summaryPng = "~{outputDir}/~{sampleName}.cup.report.summary.png"
        File? featuresPng = "~{outputDir}/~{sampleName}.cup.report.features.png"
        File reportPdf = "~{outputDir}/~{sampleName}_cup_report.pdf"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        sampleName: {description: "The sample id.", category: "required"}
        cupData: {description: "The output produced by cuppa.", category: "required"}
        outputDir: {description: "The directory the ouput will be placed in.", category: "common"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Cuppa {
    input {
        Array[File]+ linxOutput
        Array[File]+ purpleOutput
        File virusInterpreterOutput
        String sampleName
        Array[String]+ categories = ["DNA"]
        Array[File]+ referenceData
        String outputDir = "./cuppa"

        String javaXmx = "4G"
        String memory = "5GiB"
        Int timeMinutes = 10
        String dockerImage = "quay.io/biowdl/cuppa@sha256:e76d367a3226068967fb64ad6adaa889cbdcc01397075b0cbc382bbba4350b98"
    }

    command {
        set -e
        mkdir -p sampleData ~{outputDir}
        ln -s -t sampleData ~{sep=" " linxOutput} ~{sep=" " purpleOutput}
        ln -s -t sampleData ~{virusInterpreterOutput}
        cuppa -Xmx~{javaXmx} \
        -output_dir ~{outputDir} \
        -categories '~{sep="," categories}' \
        -ref_data_dir ~{sub(referenceData[0], basename(referenceData[0]), "")} \
        -sample_data_dir sampleData \
        -sample_data ~{sampleName}
    }

    output {
        File cupData = "~{outputDir}/~{sampleName}.cup.data.csv"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        linxOutput: {description: "The files produced by linx.", category: "required"}
        purpleOutput: {description: "The files produced by purple.", category: "required"}
        sampleName: {description: "The name of the sample.", category: "required"}
        categories: {description: "The classifiers to use.", category: "advanced"}
        referenceData : {description: "The reference data.", category: "required"}
        outputDir: {description: "The directory the ouput will be placed in.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task CuppaChart {
    input {
        String sampleName
        File cupData
        String outputDir = "./cuppa"

        String memory = "4GiB"
        Int timeMinutes = 5
        String dockerImage = "quay.io/biowdl/cuppa@sha256:e76d367a3226068967fb64ad6adaa889cbdcc01397075b0cbc382bbba4350b98"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        cuppa-chart \
        -sample ~{sampleName} \
        -sample_data ~{cupData} \
        -output_dir ~{outputDir}
    }

    output {
        File cuppaChart = "~{outputDir}/~{sampleName}.cuppa.chart.png"
        File cuppaConclusion = "~{outputDir}/~{sampleName}.cuppa.conclusion.txt"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        sampleName: {description: "The name of the sample.", category:"common"}
        cupData: {description: "The cuppa output.", category: "required"}
        outputDir: {description: "The directory the output will be written to.", category:"common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Gripss {
    input {
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        File knownFusionPairBedpe
        File breakendPon
        File breakpointPon
        File repeatMaskFile
        String? referenceName
        String sampleName
        File vcf
        File vcfIndex
        String outputId
        String outputDir = "./"
        Boolean hg38 = false
        Int? hardMinTumorQual
        Int? minQualBreakPoint
        Int? minQualBreakEnd
        Boolean filterSgls = false
        Boolean germline = false

        String memory = "17GiB"
        String javaXmx = "16G"
        Int timeMinutes = 50
        String dockerImage = "quay.io/biocontainers/hmftools-gripss:2.3.2--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        gripss -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -ref_genome ~{referenceFasta} \
        -ref_genome_version ~{if hg38 then "38" else "37"} \
        -known_hotspot_file ~{knownFusionPairBedpe} \
        -pon_sgl_file ~{breakendPon} \
        -pon_sv_file ~{breakpointPon} \
        -repeat_mask_file ~{repeatMaskFile} \
        ~{"-reference " + referenceName} \
        -sample ~{sampleName} \
        -vcf ~{vcf} \
        -output_dir ~{outputDir} \
        -output_id ~{outputId} \
        ~{if filterSgls then "-filter_sgls" else ""} \
        ~{"-hard_min_tumor_qual " + hardMinTumorQual} \
        ~{"-min_qual_break_point " + minQualBreakPoint} \
        ~{"-min_qual_break_end " + minQualBreakEnd} \
        ~{if germline then "-germline" else ""}
    }

    String suffix = if defined(referenceName) then "somatic" else "germline"

    output {
        File fullVcf = "~{outputDir}/~{sampleName}.gripss.~{suffix}.vcf.gz"
        File fullVcfIndex = "~{outputDir}/~{sampleName}.gripss.~{suffix}.vcf.gz.tbi"
        File filteredVcf = "~{outputDir}/~{sampleName}.gripss.filtered.~{suffix}.vcf.gz"
        File filteredVcfIndex = "~{outputDir}/~{sampleName}.gripss.filtered.~{suffix}.vcf.gz.tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        knownFusionPairBedpe: {description: "Equivalent to the `-known_hotspot_file` option.", category: "required"}
        breakendPon: {description: "Equivalent to the `-pon_sgl_file` option.", category: "required"}
        breakpointPon: {description: "Equivalent to the `-pon_sv_file` option.", category: "required"}
        sampleName: {description: "The name of the tumor sample.", category: "required"}
        referenceName: {description: "The name of the normal sample.", category: "required"}
        vcf: {description: "The input VCF.", category: "required"}
        vcfIndex: {description: "The index for the input VCF.", category: "required"}
        outputDir: {description: "The path the output will be written to.", category:"required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task HealthChecker {
    input {
        String outputDir = "."
        String referenceName
        File referenceFlagstats
        File referenceMetrics
        String tumorName
        File tumorFlagstats
        File tumorMetrics
        Array[File]+ purpleOutput

        String javaXmx = "2G"
        String memory = "3GiB"
        Int timeMinutes = 1
        String dockerImage = "quay.io/biowdl/health-checker:3.4"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        health-checker -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -reference ~{referenceName} \
        -ref_flagstat_file ~{referenceFlagstats} \
        -ref_wgs_metrics_file ~{referenceMetrics} \
        -tumor ~{tumorName} \
        -tum_flagstat_file ~{tumorFlagstats} \
        -tum_wgs_metrics_file ~{tumorMetrics} \
        -purple_dir ~{sub(purpleOutput[0], basename(purpleOutput[0]), "")} \
        -output_dir ~{outputDir}
        if [ -e '~{outputDir}/~{tumorName}.HealthCheckSucceeded' ]
          then
            echo 'true' > '~{outputDir}/succeeded'
        fi
        if [ -e '~{outputDir}/~{tumorName}.HealthCheckFailed' ]
          then
            echo 'false' > '~{outputDir}/succeeded'
        fi
    }

    output {
        Boolean succeeded = read_boolean("succeeded")
        File outputFile = if succeeded
                          then "~{outputDir}/~{tumorName}.HealthCheckSucceeded"
                          else "~{outputDir}/~{tumorName}.HealthCheckFailed"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        outputDir: {description: "The path the output will be written to.", category:"required"}
        referenceName: {description: "The name of the normal sample.", category: "required"}
        referenceFlagstats: {description: "The flagstats for the normal sample.", category: "required"}
        referenceMetrics: {description: "The picard WGS metrics for the normal sample.", category: "required"}
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        tumorFlagstats: {description: "The flagstats for the tumor sample.", category: "required"}
        tumorMetrics: {description: "The picard WGS metrics for the tumor sample.", category: "required"}
        purpleOutput: {description: "The files from purple's output directory.", category: "required"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Isofox {
    input {
        String sampleName
        File neoepitopeFile
        File bamFile
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String refGenomeVersion
        File expCountsFile
        File expGcRatiosFile

        String outputDir = "./isofox"
        Int readLength = 151

        #The following should be in the same directory.
        File geneDataCsv
        File proteinFeaturesCsv
        File transExonDataCsv
        File transSpliceDataCsv

        Int threads = 10
        String javaXmx = "12G"
        String memory = "13GiB"
        Int timeMinutes = 120
        String dockerImage = "quay.io/biocontainers/hmftools-isofox:1.6.2--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        sed 's/\t/,/g' ~{neoepitopeFile} > tmp.neo_data.csv
        isofox -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -sample ~{sampleName} \
        -functions 'NEO_EPITOPES;TRANSCRIPT_COUNTS;ALT_SPLICE_JUNCTIONS;FUSIONS' \
        -neoepitope_file tmp.neo_data.csv \
        -bam_file ~{bamFile} \
        -ref_genome ~{referenceFasta} \
        -ref_genome_version ~{refGenomeVersion} \
        -ensembl_data_dir ~{sub(geneDataCsv, basename(geneDataCsv), "")} \
        -output_dir ~{outputDir} \
        -log_debug \
        -threads ~{threads}
    }

    output {
        File neoepitopeTsv = "~{outputDir}/~{sampleName}.isf.neoepitope.tsv"
        Array[File] outputs = [neoepitopeTsv]
        #TODO
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        sampleName: {description: "The name of the sample.", category: "required"}
        neoepitopeFile: {description: "Neo's data file.", category: "required"}
        bamFile: {description: "Input rna BAM file.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}  
        refGenomeVersion: {description: "The version of the genome assembly used for alignment. Either \"37\" or \"38\".", category: "required"}
        expCountsFile: {description: "Isofox reference file.", category: "required"}
        expGcRatiosFile: {description: "Isofox reference file.", category: "required"}
        outputDir: {description: "The directory the outputs will be written to.", category: "required"}
        geneDataCsv: {description: "A  CSV file containing gene information, must be in the same directory as `proteinFeaturesCsv`, `transExonDataCsv` and `transSpliceDataCsv`.", category: "required"}
        proteinFeaturesCsv: {description: "A  CSV file containing protein feature information, must be in the same directory as `geneDataCsv`, `transExonDataCsv` and `transSpliceDataCsv`.", category: "required"}
        transExonDataCsv: {description: "A  CSV file containing transcript exon information, must be in the same directory as `geneDataCsv`, `proteinFeaturesCsv` and `transSpliceDataCsv`.", category: "required"}
        transSpliceDataCsv: {description: "A  CSV file containing transcript splicing information, must be in the same directory as `geneDataCsv`, `proteinFeaturesCsv` and `transExonDataCsv`.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Lilac {
    input {
        String sampleName
        File referenceBam
        File referenceBamIndex
        File? tumorBam
        File? tumorBamIndex
        String refGenomeVersion
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        File? geneCopyNumberFile
        File? somaticVariantsFile
        File? somaticVariantsFileIndex
        String outputDir = "./lilac"

        #The following need to be in the same directory
        File hlaRefAminoacidSequencesCsv
        File hlaRefNucleotideSequencesCsv
        File lilacAlleleFrequenciesCsv

        String javaXmx = "15G"
        String memory = "16GiB"
        Int timeMinutes = 1440 #FIXME
        Int threads = 1
        String dockerImage = "quay.io/biocontainers/hmftools-lilac:1.4.2--hdfd78af_0"
    }

    command {
        LILAC -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -sample ~{sampleName} \
        -reference_bam ~{referenceBam} \
        -ref_genome ~{referenceFasta} \
        -ref_genome_version ~{refGenomeVersion} \
        -resource_dir ~{sub(hlaRefAminoacidSequencesCsv, basename(hlaRefAminoacidSequencesCsv), "")} \
        -output_dir ~{outputDir} \
        -threads ~{threads} \
        ~{"-tumor_bam " + tumorBam} \
        ~{"-gene_copy_number " + geneCopyNumberFile} \
        ~{"-somatic_vcf " + somaticVariantsFile}
    }

    output {
        File lilacCsv = "~{outputDir}/~{sampleName}.lilac.csv"
        File lilacQcCsv = "~{outputDir}/~{sampleName}.lilac.qc.csv"
        File candidatesCoverageCsv = "~{outputDir}/~{sampleName}.candidates.coverage.csv"
        Array[File] outputs = [lilacCsv, lilacQcCsv, candidatesCoverageCsv]
    }

    runtime {
        memory: memory
        cpu: threads
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        sampleName: {description: "The name of the sample.", category: "required"}
        referenceBam: {description: "The bam file for the reference sample.", category: "required"}
        referenceBamIndex: {description: "The index for the reference sample's bam file.", category: "required"}
        tumorBam: {description: "The bam file for the tumor sample.", category: "common"}
        tumorBamIndex: {description: "The index for the tumor sample's bam file.", category: "required"}
        refGenomeVersion: {description: "The version of the genome assembly used for alignment. Either \"37\" or \"38\".", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        geneCopyNumberFile: {description: "Gene copy number file produced by purple.", category: "common"}
        somaticVariantsFile: {description: "Somatic variant VCF produced by purple.", category: "common"}
        somaticVariantsFileIndex: {description: "Index for the somatic variant VCf produced by purple.", category: "common"}
        outputDir: {description: "The directory the outputs will be written to.", category: "required"}
        hlaRefAminoacidSequencesCsv: {description: "LILAC reference file.", category: "required"}
        hlaRefNucleotideSequencesCsv: {description: "LILAC reference file.", category: "required"}
        lilacAlleleFrequenciesCsv: {description: "LILAC reference file.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        threads: {description: "The number of threads to use", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}


task Linx {
    input {
        String sampleName
        File svVcf
        File svVcfIndex
        Array[File] purpleOutput = []
        String refGenomeVersion
        String outputDir = "./linx"
        File? fragileSiteCsv
        File lineElementCsv
        File? knownFusionCsv
        File driverGenePanel
        Boolean writeAllVisFusions = false
        Boolean germline = false
        Boolean checkFusions = true
        Boolean checkDrivers = true
        Boolean writeVisData = true
        Boolean writeNeoEpitopes = false
        #The following should be in the same directory.
        File geneDataCsv
        File proteinFeaturesCsv
        File transExonDataCsv
        File transSpliceDataCsv

        String memory = "9GiB"
        String javaXmx = "8G"
        Int timeMinutes = 10
        String dockerImage = "quay.io/biocontainers/hmftools-linx:1.22.1--hdfd78af_0"

        String? DONOTDEFINE
    }

    String? purpleDir = if length(purpleOutput) > 0
        then sub(purpleOutput[0], basename(purpleOutput[0]), "")
        else DONOTDEFINE

    command {
        linx -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -sample ~{sampleName} \
        -sv_vcf ~{svVcf} \
        ~{"-purple_dir " + purpleDir} \
        -ref_genome_version ~{refGenomeVersion} \
        -output_dir ~{outputDir} \
        ~{"-fragile_site_file " + fragileSiteCsv} \
        -line_element_file ~{lineElementCsv} \
        -ensembl_data_dir ~{sub(geneDataCsv, basename(geneDataCsv), "")} \
        ~{if checkFusions then "-check_fusions" else ""} \
        ~{"-known_fusion_file " + knownFusionCsv} \
        ~{if checkDrivers then "-check_drivers" else ""} \
        -driver_gene_panel ~{driverGenePanel} \
        ~{if writeVisData then "-write_vis_data" else ""} \
        ~{if writeAllVisFusions then "-write_all_vis_fusions" else ""} \
        ~{if writeNeoEpitopes then "-write_neo_epitopes" else ""} \
        ~{if germline then "-germline" else ""}
    }

    String prefix = if germline then "~{sampleName}.linx.germline" else "~{sampleName}.linx"

    output {
        File driverCatalog = "~{outputDir}/~{prefix}.driver.catalog.tsv"
        File linxClusters = "~{outputDir}/~{prefix}.clusters.tsv"
        File linxLinks = "~{outputDir}/~{prefix}.links.tsv"
        File linxSvs = "~{outputDir}/~{prefix}.svs.tsv"
        File? linxBreakend = "~{outputDir}/~{prefix}.breakend.tsv"
        File? linxDrivers = "~{outputDir}/~{prefix}.drivers.tsv"
        File? linxFusion = "~{outputDir}/~{prefix}.fusion.tsv"
        File? linxVisCopyNumber = "~{outputDir}/~{prefix}.vis_copy_number.tsv"
        File? linxVisFusion = "~{outputDir}/~{prefix}.vis_fusion.tsv"
        File? linxVisGeneExon = "~{outputDir}/~{prefix}.vis_gene_exon.tsv"
        File? linxVisProteinDomain = "~{outputDir}/~{prefix}.vis_protein_domain.tsv"
        File? linxVisSegments = "~{outputDir}/~{prefix}.vis_segments.tsv"
        File? linxVisSvData = "~{outputDir}/~{prefix}.vis_sv_data.tsv"
        File? linxDisruptionTsv = "~{outputDir}/~{prefix}.disruption.tsv"
        File? linxNeoepitopeTsv = "~{outputDir}/~{prefix}.neoepitope.tsv"
        File linxVersion = "~{outputDir}/linx.version"
        Array[File] outputs = select_all([driverCatalog, linxBreakend, linxClusters, linxDrivers, linxFusion,
                               linxLinks, linxSvs, linxVisCopyNumber, linxVisFusion,
                               linxVisGeneExon, linxVisProteinDomain, linxVisSegments, linxVisSvData,
                               linxDisruptionTsv, linxNeoepitopeTsv, linxVersion])
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        sampleName: {description: "The name of the sample.", category: "required"}
        svVcf: {description: "A VCF file containing structural variants, produced using GRIDSS, annotated for viral insertions and postprocessed with GRIPSS.", category: "required"}
        svVcfIndex: {description: "Index for the structural variants VCf file.", category: "required"}
        purpleOutput: {description: "The files produced by PURPLE.", category: "required"}
        refGenomeVersion: {description: "The version of the genome assembly used for alignment. Either \"37\" or \"38\".", category: "required"}
        outputDir: {description: "The directory the outputs will be written to.", category: "required"}
        fragileSiteCsv: {description: "A list of known fragile sites.", category: "required"}
        lineElementCsv: {description: "A list of known LINE source regions.", category: "required"}
        knownFusionCsv: {description: "A CSV file describing known fusions.", category: "required"}
        driverGenePanel: {description: "A TSV file describing the driver gene panel.", category: "required"}
        writeAllVisFusions: {description: "Equivalent to the -write_all_vis_fusions flag.", category: "advanced"}
        writeNeoEpitopes: {description: "Equivalent to the -write_neo_epitopes flag.", category: "advanced"}
        geneDataCsv: {description: "A  CSV file containing gene information, must be in the same directory as `proteinFeaturesCsv`, `transExonDataCsv` and `transSpliceDataCsv`.", category: "required"}
        proteinFeaturesCsv: {description: "A  CSV file containing protein feature information, must be in the same directory as `geneDataCsv`, `transExonDataCsv` and `transSpliceDataCsv`.", category: "required"}
        transExonDataCsv: {description: "A  CSV file containing transcript exon information, must be in the same directory as `geneDataCsv`, `proteinFeaturesCsv` and `transSpliceDataCsv`.", category: "required"}
        transSpliceDataCsv: {description: "A  CSV file containing transcript splicing information, must be in the same directory as `geneDataCsv`, `proteinFeaturesCsv` and `transExonDataCsv`.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task LinxVisualisations {
    input {
        String outputDir = "./linx"
        String sample
        String refGenomeVersion
        Array[File]+ linxOutput
        Boolean plotReportable = true

        String memory = "9GiB"
        String javaXmx = "8G"
        Int timeMinutes = 1440
        String dockerImage = "quay.io/biocontainers/hmftools-linx:1.22.1--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        java -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -cp /usr/local/share/hmftools-linx-1.22.1-0/linx.jar \
        com.hartwig.hmftools.linx.visualiser.SvVisualiser \
        -sample ~{sample} \
        -ref_genome_version ~{refGenomeVersion} \
        -circos /usr/local/bin/circos \
        -vis_file_dir ~{sub(linxOutput[0], basename(linxOutput[0]), "")} \
        -data_out ~{outputDir}/circos \
        -plot_out ~{outputDir}/plots \
        ~{if plotReportable then "-plot_reportable" else ""}
    }

    output {
        Array[File] circos = glob("~{outputDir}/circos/*")
        Array[File] plots = glob("~{outputDir}/plots/*")
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        outputDir: {description: "The directory the outputs will be written to.", category: "required"}
        sample: {description: "The sample's name.", category: "required"}
        refGenomeVersion: {description: "The version of the genome assembly used for alignment. Either \"37\" or \"38\".", category: "required"}
        linxOutput: {description: "The directory containing the linx output.", category: "required"}
        plotReportable: {description: "Equivalent to the -plot_reportable flag.", category: "advanced"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Neo {
    input {
        String sampleId
        File somaticVcf
        File somaticVcfIndex
        Array[File]+ linxOutput
        String refGenomeVersion
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String outputDir = "./neo"
        #The following should be in the same directory.
        File geneDataCsv
        File proteinFeaturesCsv
        File transExonDataCsv
        File transSpliceDataCsv

        Int reqAminoAcids = 15

        String memory = "9GiB"
        String javaXmx = "8G"
        Int timeMinutes = 1440
        String dockerImage = "quay.io/biocontainers/hmftools-neo:1.0.1--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        neo -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -sample ~{sampleId} \
        -ref_genome_version ~{refGenomeVersion} \
        -ref_genome ~{referenceFasta} \
        -ensembl_data_dir ~{sub(geneDataCsv, basename(geneDataCsv), "")} \
        -linx_dir ~{sub(linxOutput[0], basename(linxOutput[0]), "")} \
        -somatic_vcf ~{somaticVcf} \
        -req_amino_acids ~{reqAminoAcids} \
        -output_dir ~{outputDir}
    }

    output {
        File neoData = "~{outputDir}/~{sampleId}.neo.neo_data.tsv"
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        sampleId: {description: "The name/id of the sample.", category: "required"}
        somaticVcf: {description: "The vcf containing the samples's somatic variants.", category: "required"}
        somaticVcfIndex: {description: "The vcf containing the samples's somatic variants.", category: "required"}
        linxOutput: {description: "The directory containing the linx output.", category: "required"}
        refGenomeVersion: {description: "The version of the genome assembly used for alignment. Either \"37\" or \"38\".", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        outputDir: {description: "The directory the outputs will be written to.", category: "common"}
        geneDataCsv: {description: "A  CSV file containing gene information, must be in the same directory as `proteinFeaturesCsv`, `transExonDataCsv` and `transSpliceDataCsv`.", category: "required"}
        proteinFeaturesCsv: {description: "A  CSV file containing protein feature information, must be in the same directory as `geneDataCsv`, `transExonDataCsv` and `transSpliceDataCsv`.", category: "required"}
        transExonDataCsv: {description: "A  CSV file containing transcript exon information, must be in the same directory as `geneDataCsv`, `proteinFeaturesCsv` and `transSpliceDataCsv`.", category: "required"}
        transSpliceDataCsv: {description: "A  CSV file containing transcript splicing information, must be in the same directory as `geneDataCsv`, `proteinFeaturesCsv` and `transExonDataCsv`.", category: "required"}
        reqAminoAcids: {description: "Equivalent to neo's -req_amino_acids option.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task NeoScorer {
    input {
        String sampleId
        Array[File]+ neoBindingFiles
        String neoBindingFileId
        File cancerTpmMedians
        File neoData
        Array[File]+ lilacOutput
        Array[File]+ purpleOutput
        String outputDir = "./neo"

        #The following should be in the same directory.
        File geneDataCsv
        File proteinFeaturesCsv
        File transExonDataCsv
        File transSpliceDataCsv

        String? cancerType
        Array[File]? isofoxOutput
        File? rnaSomaticVcf
        File? rnaSomaticVcfIndex

        String memory = "9GiB"
        String javaXmx = "8G"
        Int timeMinutes = 1440
        String dockerImage = "quay.io/biocontainers/hmftools-neo:1.0.1--hdfd78af_0"
    }

    String isofoxDir = sub(select_first([isofoxOutput, [""]])[0], basename(select_first([isofoxOutput, [""]])[0]), "")

    command {
        set -e
        mkdir -p ~{outputDir}
        neo com.hartwig.hmftools.neo.scorer.NeoScorer Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -sample ~{sampleId} \
        ~{"-cancer_type " + cancerType} \
        -ensembl_data_dir ~{sub(geneDataCsv, basename(geneDataCsv), "")} \
        -score_file_dir ~{sub(neoBindingFiles[0], basename(neoBindingFiles[0]), "")} \
        -score_file_id ~{neoBindingFileId} \
        -cancer_tpm_medians_file ~{cancerTpmMedians} \
        -neo_dir ~{sub(neoData, basename(neoData), "")} \
        ~{if defined(isofoxOutput) then "-isofox_dir " + isofoxDir else ""} \
        -lilac_dir ~{sub(lilacOutput[0], basename(lilacOutput[0]), "")} \
        -purple_dir ~{sub(purpleOutput[0], basename(purpleOutput[0]), "")} \
        ~{"-rna_somatic_vcf " + rnaSomaticVcf} \
        -output_dir ~{outputDir}
    }

    output {
        File neoepitopes = "~{outputDir}/~{sampleId}.neo.neoepitope.tsv"
        File peptideScores = "~{outputDir}/~{sampleId}.neo.peptide_scores.tsv"
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        sampleId: {description: "The name/id of the sample.", category: "required"}
        neoBindingFiles: {description: "The neo binding reference files.", category: "required"}
        neoBindingFileId: {description: "The neo binding reference file version id.", category: "required"}
        cancerTpmMedians: {description: "HMF RNA cohort transcript median TPM file.", category: "required"}
        neoData: {description: "Data file produced by neo.", category: "required"}
        lilacOutput: {description: "The output produced by lilac.", category: "required"}
        purpleOutput: {description: "The output produced by purple.", category: "required"}
        outputDir: {description: "The directory the outputs will be written to.", category: "required"}
        geneDataCsv: {description: "A  CSV file containing gene information, must be in the same directory as `proteinFeaturesCsv`, `transExonDataCsv` and `transSpliceDataCsv`.", category: "required"}
        proteinFeaturesCsv: {description: "A  CSV file containing protein feature information, must be in the same directory as `geneDataCsv`, `transExonDataCsv` and `transSpliceDataCsv`.", category: "required"}
        transExonDataCsv: {description: "A  CSV file containing transcript exon information, must be in the same directory as `geneDataCsv`, `proteinFeaturesCsv` and `transSpliceDataCsv`.", category: "required"}
        transSpliceDataCsv: {description: "A  CSV file containing transcript splicing information, must be in the same directory as `geneDataCsv`, `proteinFeaturesCsv` and `transExonDataCsv`.", category: "required"}
        cancerType: {description: "The cancer type.", category: "common"}
        isofoxOutput: {description: "The output produced by isofox.", category: "common"}
        rnaSomaticVcf: {description: "SageAppend produced rna somatic VCF file.", category: "common"}
        rnaSomaticVcfIndex: {description: "Index for the rna somatic VCF file.", category: "common"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Orange {
    input {
        String outputDir = "./orange"
        File doidJson
        Array[String] sampleDoids
        String tumorName
        String referenceName
        File referenceMetrics
        File tumorMetrics
        File referenceFlagstats
        File tumorFlagstats
        File sageGermlineGeneCoverageTsv
        File sageSomaticRefSampleBqrPlot
        File sageSomaticTumorSampleBqrPlot
        File purpleGeneCopyNumberTsv
        File purpleGermlineDriverCatalogTsv
        File purpleGermlineDeletionTsv
        File purpleGermlineVariantVcf
        File purpleGermlineVariantVcfIndex
        Array[File]+ purplePlots
        File purplePurityTsv
        File purpleQcFile
        File purpleSomaticCopyNumberFile
        File purpleSomaticDriverCatalogTsv
        File purpleSomaticVariantVcf
        File purpleSomaticVariantVcfIndex
        File lilacQcCsv
        File lilacResultCsv
        File linxFusionTsv
        File linxBreakendTsv
        File linxDriverCatalogTsv
        File linxDriverTsv
        File linxGermlineDisruptionTsv
        Array[File]+ linxPlots
        File linxStructuralVariantTsv
        File cuppaResultCsv
        File cuppaSummaryPlot
        File? cuppaFeaturePlot
        File chordPredictionTxt
        File peachGenotypeTsv
        File protectEvidenceTsv
        File annotatedVirusTsv
        #File pipelineVersionFile
        File cohortMappingTsv
        File cohortPercentilesTsv
        Boolean hg38 = false
        File driverGenePanel
        File knownFusionFile

        String memory = "17GiB"
        String javaXmx = "16G"
        Int timeMinutes = 10
        String dockerImage = "quay.io/biocontainers/hmftools-orange:1.10.2--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        export JAVA_TOOL_OPTIONS='--add-opens=java.base/java.time=ALL-UNNAMED'
        orange -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -output_dir ~{outputDir} \
        -ref_genome_version ~{if hg38 then "38" else "37"} \
        -doid_json ~{doidJson} \
        -primary_tumor_doids '~{sep=";" sampleDoids}' \
        -max_evidence_level C \
        -tumor_sample_id ~{tumorName} \
        -reference_sample_id ~{referenceName} \
        -ref_sample_wgs_metrics_file ~{referenceMetrics} \
        -tumor_sample_wgs_metrics_file ~{tumorMetrics} \
        -ref_sample_flagstat_file ~{referenceFlagstats} \
        -tumor_sample_flagstat_file ~{tumorFlagstats} \
        -sage_germline_gene_coverage_tsv ~{sageGermlineGeneCoverageTsv} \
        -sage_somatic_ref_sample_bqr_plot ~{sageSomaticRefSampleBqrPlot} \
        -sage_somatic_tumor_sample_bqr_plot ~{sageSomaticTumorSampleBqrPlot} \
        -purple_gene_copy_number_tsv ~{purpleGeneCopyNumberTsv} \
        -purple_germline_driver_catalog_tsv ~{purpleGermlineDriverCatalogTsv} \
        -purple_germline_deletion_tsv ~{purpleGermlineDeletionTsv} \
        -purple_germline_variant_vcf ~{purpleGermlineVariantVcf} \
        -purple_plot_directory ~{sub(purplePlots[0], basename(purplePlots[0]), "")} \
        -purple_purity_tsv ~{purplePurityTsv} \
        -purple_qc_file ~{purpleQcFile} \
        -purple_somatic_copy_number_tsv ~{purpleSomaticCopyNumberFile} \
        -purple_somatic_driver_catalog_tsv ~{purpleSomaticDriverCatalogTsv} \
        -purple_somatic_variant_vcf ~{purpleSomaticVariantVcf} \
        -lilac_qc_csv ~{lilacQcCsv} \
        -lilac_result_csv ~{lilacResultCsv} \
        -linx_fusion_tsv ~{linxFusionTsv} \
        -linx_breakend_tsv ~{linxBreakendTsv} \
        -linx_driver_catalog_tsv ~{linxDriverCatalogTsv} \
        -linx_driver_tsv ~{linxDriverTsv} \
        -linx_germline_disruption_tsv ~{linxGermlineDisruptionTsv} \
        -linx_plot_directory ~{sub(linxPlots[0], basename(linxPlots[0]), "")} \
        -linx_structural_variant_tsv ~{linxStructuralVariantTsv} \
        -cuppa_result_csv ~{cuppaResultCsv} \
        -cuppa_summary_plot ~{cuppaSummaryPlot} \
        ~{"-cuppa_feature_plot " + cuppaFeaturePlot} \
        -chord_prediction_txt ~{chordPredictionTxt} \
        -peach_genotype_tsv ~{peachGenotypeTsv} \
        -protect_evidence_tsv ~{protectEvidenceTsv} \
        -annotated_virus_tsv ~{annotatedVirusTsv} \
        -cohort_mapping_tsv ~{cohortMappingTsv} \
        -cohort_percentiles_tsv ~{cohortPercentilesTsv} \
        -driver_gene_panel_tsv ~{driverGenePanel} \
        -known_fusion_file ~{knownFusionFile}
    }

    output {
        File orangeJson = "~{outputDir}/~{tumorName}.orange.json"
        File orangePdf = "~{outputDir}/~{tumorName}.orange.pdf"
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        outputDir: {description: "The directory the outputs will be written to.", category: "common"}
        doidJson: {description: "A json with the DOID (Human Disease Ontology) tree.", category: "required"}
        sampleDoids: {description: "The DOIDs (Human Disease Ontology) for the primary tumor.", category: "required"}
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        referenceName: {description: "The name of the normal sample.", category: "required"}
        referenceMetrics: {description: "The picard WGS metrics for the normal sample.", category: "required"}
        tumorMetrics: {description: "The picard WGS metrics for the tumor sample.", category: "required"}
        referenceFlagstats: {description: "The flagstats for the normal sample.", category: "required"}
        tumorFlagstats: {description: "The flagstats for the tumor sample.", category: "required"}
        sageGermlineGeneCoverageTsv: {description: "Gene coverage file produced by the germline sage run.", category: "required"}
        sageSomaticRefSampleBqrPlot: {description: "The reference bqr plot produced by the somatic sage run.", category: "required"}
        sageSomaticTumorSampleBqrPlot: {description: "The reference bqr plot produced by the somatic sage run.", category: "required"}
        purpleGeneCopyNumberTsv: {description: "Copy number tsv produced by purple.", category: "required"}
        purpleGermlineDriverCatalogTsv: {description: "Germline driver catalog produced by purple.", category: "required"}
        purpleGermlineVariantVcf: {description: "Germline variant vcf produced by purple.", category: "required"}
        purplePlots: {description: "The plots generated by purple.", category: "required"}
        purplePurityTsv: {description: "The purity file produced by purple.", category: "required"}
        purpleQcFile: {description: "The qc file produced by purple.", category: "required"}
        purpleSomaticDriverCatalogTsv: {description: "Somatic driver catalog produced by purple.", category: "required"}
        purpleSomaticVariantVcf: {description: "Somatic variant vcf produced by purple.", category: "required"}
        linxFusionTsv: {description: "The fusions tsv produced by linx.", category: "required"}
        linxBreakendTsv: {description: "The breakend tsv produced by linx.", category: "required"}
        linxDriverCatalogTsv: {description: "The driver catalog produced by linx.", category: "required"}
        linxDriverTsv: {description: "The driver tsv produced by linx.", category: "required"}
        linxPlots: {description: "The plots generated by linx.", category: "required"}
        cuppaResultCsv: {description: "The cuppa results csv.", category: "required"}
        cuppaSummaryPlot: {description: "The cuppa summary plot.", category: "required"}
        cuppaFeaturePlot: {description: "The cuppa feature plot.", category: "common"}
        chordPredictionTxt: {description: "Chord prediction results.", category: "required"}
        peachGenotypeTsv: {description: "Genotype tsv produced by peach.", category: "required"}
        protectEvidenceTsv: {description: "Evidence tsv produced by protect.", category: "required"}
        annotatedVirusTsv: {description: "Annotated virus tsv produced by virus-interpreter.", category: "required"}
        #pipelineVersionFile: {description: "", category: "required"}
        cohortMappingTsv: {description: "Cohort mapping file from the HMFTools resources.", category: "required"}
        cohortPercentilesTsv: {description: "Cohort percentile file from the HMFTools resources.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Pave {
    input {
        String outputDir = "./"
        String sampleName
        File vcfFile
        File vcfFileIndex
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String refGenomeVersion
        File driverGenePanel
        File mappabilityBed
        Array[File] gnomadFreqFiles = []

        #The following should be in the same directory.
        File geneDataCsv
        File proteinFeaturesCsv
        File transExonDataCsv
        File transSpliceDataCsv

        File? ponFile
        File? ponArtefactFile
        String? ponFilters
        File? clinvarVcf
        File? clinvarVcfIndex
        File? blacklistVcf
        File? blacklistBed
        File? blacklistVcfIndex
        Boolean writePassOnly = false

        Int timeMinutes = 50
        String javaXmx = "8G"
        String memory = "9GiB"
        String dockerImage = "quay.io/biocontainers/hmftools-pave:1.4.1--hdfd78af_0"

        String? DONOTDEFINE
    }

    String? gnomadFreqDir = if length(gnomadFreqFiles) > 0
        then sub(gnomadFreqFiles[0], basename(gnomadFreqFiles[0]), "")
        else DONOTDEFINE

    command {
        set -e
        mkdir -p ~{outputDir}
        pave -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -sample ~{sampleName} \
        -vcf_file ~{vcfFile} \
        -output_dir ~{outputDir} \
        -ensembl_data_dir ~{sub(geneDataCsv, basename(geneDataCsv), "")} \
        -ref_genome ~{referenceFasta} \
        -ref_genome_version ~{refGenomeVersion} \
        -driver_gene_panel ~{driverGenePanel} \
        -read_pass_only \
        -mappability_bed ~{mappabilityBed} \
        ~{"-pon_file " + ponFile} \
        ~{"-pon_artefact_file " + ponArtefactFile} \
        ~{if defined(ponFilters) then ("-pon_filters '" + ponFilters + "'") else ""} \
        ~{"-gnomad_freq_dir " + gnomadFreqDir} \
        ~{if defined(gnomadFreqDir) then "-gnomad_load_chr_on_demand" else ""} \
        ~{"-clinvar_vcf " + clinvarVcf} \
        ~{"-blacklist_bed " + blacklistBed} \
        ~{"-blacklist_vcf " + blacklistVcf} \
        ~{if writePassOnly then "-write_pass_only" else ""}
    }

    output {
        File outputVcf = "~{outputDir}/~{sub(basename(vcfFile), 'vcf.gz$', 'pave.vcf.gz')}"
        File outputVcfIndex = "~{outputDir}/~{sub(basename(vcfFile), 'vcf.gz$', 'pave.vcf.gz.tbi')}"
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        outputDir: {description: "The directory the outputs will be written to.", category: "required"}
        sampleName: {description: "The name of the sample.", category: "required"}
        vcfFile: {description: "The input VCF file.", category: "required"}
        vcfFileIndex: {description: "The index for the input vcf file.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        refGenomeVersion: {description: "The version of the genome assembly used for alignment. Either \"HG19\" or \"HG38\".", category: "required"}
        driverGenePanel: {description: "A TSV file describing the driver gene panel.", category: "required"}
        geneDataCsv: {description: "A  CSV file containing gene information, must be in the same directory as `proteinFeaturesCsv`, `transExonDataCsv` and `transSpliceDataCsv`.", category: "required"}
        proteinFeaturesCsv: {description: "A  CSV file containing protein feature information, must be in the same directory as `geneDataCsv`, `transExonDataCsv` and `transSpliceDataCsv`.", category: "required"}
        transExonDataCsv: {description: "A  CSV file containing transcript exon information, must be in the same directory as `geneDataCsv`, `proteinFeaturesCsv` and `transSpliceDataCsv`.", category: "required"}
        transSpliceDataCsv: {description: "A  CSV file containing transcript splicing information, must be in the same directory as `geneDataCsv`, `proteinFeaturesCsv` and `transExonDataCsv`.", category: "required"}
        mappabilityBed: {description: "A bed file with mappability information.", category: "required"}
        ponFile: {description: "A panel of normals files.", category: "common"}
        ponArtefactFile: {description: "A panel of normals artefact file.", category: "common"}
        ponFilters: {description: "Filters to be applied based on the panel of normals.", category: "common"}
        gnomadFreqFiles: {description: "A directory with gnomad frequency information.", category: "common"}
        clinvarVcf: {description: "A clinvar VCF file.", category: "common"}
        clinvarVcfIndex: {description: "The index for the clinvar VCF file.", category: "common"}
        blacklistVcf: {description: "A blacklist VCF file.", category: "common"}
        blacklistBed: {description: "A blacklist bed file.", category: "common"}
        blacklistVcfIndex: {description: "The index for the blacklist vcf file.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Protect {
    input {
        String refGenomeVersion
        String tumorName
        String referenceName
        Array[String]+ sampleDoids
        String outputDir = "./protect"
        Array[File]+ serveActionability
        File doidJson
        File purplePurity
        File purpleQc
        File purpleDriverCatalogSomatic
        File purpleDriverCatalogGermline
        File purpleSomaticVariants
        File purpleSomaticVariantsIndex
        File purpleGermlineVariants
        File purpleGermlineVariantsIndex
        File purpleGeneCopyNumber
        File linxFusion
        File linxBreakend
        File linxDriversCatalog
        File chordPrediction
        File annotatedVirus
        File lilacResultCsv
        File lilacQcCsv
        File driverGeneTsv

        String memory = "9GiB"
        String javaXmx = "8G"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/hmftools-protect:2.3--hdfd78af_0"
    }

    command {
        protect -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -ref_genome_version ~{refGenomeVersion} \
        -tumor_sample_id ~{tumorName} \
        -reference_sample_id ~{referenceName} \
        -primary_tumor_doids '~{sep=";" sampleDoids}' \
        -output_dir ~{outputDir} \
        -serve_actionability_dir ~{sub(serveActionability[0], basename(serveActionability[0]), "")} \
        -driver_gene_tsv ~{driverGeneTsv} \
        -doid_json ~{doidJson} \
        -purple_purity_tsv ~{purplePurity} \
        -purple_qc_file ~{purpleQc} \
        -purple_somatic_driver_catalog_tsv ~{purpleDriverCatalogSomatic} \
        -purple_germline_driver_catalog_tsv ~{purpleDriverCatalogGermline} \
        -purple_somatic_variant_vcf ~{purpleSomaticVariants} \
        -purple_germline_variant_vcf ~{purpleGermlineVariants} \
        -purple_gene_copy_number_tsv ~{purpleGeneCopyNumber} \
        -linx_fusion_tsv ~{linxFusion} \
        -linx_breakend_tsv ~{linxBreakend} \
        -linx_driver_catalog_tsv ~{linxDriversCatalog} \
        -chord_prediction_txt ~{chordPrediction} \
        -annotated_virus_tsv ~{annotatedVirus} \
        -lilac_result_csv ~{lilacResultCsv} \
        -lilac_qc_csv ~{lilacQcCsv}
    }

    output {
        File protectTsv = "~{outputDir}/~{tumorName}.protect.tsv"
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        refGenomeVersion: {description: "The version of the genome assembly used for alignment. Either \"37\" or \"38\".", category: "required"}
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        referenceName: {description: "The name of the normal sample.", category: "required"}
        sampleDoids: {description: "The DOIDs (Human Disease Ontology) for the primary tumor.", category: "required"}
        outputDir: {description: "The directory the outputs will be written to.", category: "required"}
        serveActionability: {description: "The actionability files generated by hmftools' serve.", category: "required"}
        doidJson: {description: "A json with the DOID (Human Disease Ontology) tree.", category: "required"}
        purplePurity: {description: "The purity file generated by purple.", category: "required"}
        purpleQc: {description: "The QC file generated by purple.", category: "required"}
        purpleDriverCatalogSomatic: {description: "The somatic driver catalog generated by purple.", category: "required"}
        purpleDriverCatalogGermline: {description: "The germline driver catalog generated by purple.", category: "required"}
        purpleSomaticVariants: {description: "The somatic VCF generated by purple.", category: "required"}
        purpleSomaticVariantsIndex: {description: "The index for the somatic VCF generated by purple.", category: "required"}
        purpleGermlineVariants: {description: "The germline VCF generated by purple.", category: "required"}
        purpleGermlineVariantsIndex: {description: "The index of the germline VCF generated by purple.", category: "required"}
        purpleGeneCopyNumber: {description: "The gene copy number file generated by purple.", category: "required"}
        linxFusion: {description: "The fusion file generated by linx.", category: "required"}
        linxBreakend: {description: "The breakend file generated by linx.", category: "required"}
        linxDriversCatalog: {description: "The driver catalog generated generated by linx.", category: "required"}
        chordPrediction: {description: "The chord prediction file.", category: "required"}
        annotatedVirus: {description: "The virus-interpreter output.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Purple {
    input {
        String? referenceName
        String tumorName
        String outputDir = "./purple"
        Array[File]+ amberOutput
        Array[File]+ cobaltOutput
        File gcProfile
        File somaticVcf
        File? germlineVcf
        File filteredSvVcf
        File filteredSvVcfIndex
        File fullSvVcf
        File fullSvVcfIndex
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String refGenomeVersion
        File driverGenePanel
        File somaticHotspots
        File? germlineHotspots
        File? germlineDelFreqFile
        Float? highlyDiploidPercentage
        Float? somaticMinPuritySpread
        File? targetRegionsBed
        File? targetRegionsRatios
        File? targetRegionsMsiIndels
        Int? minDiploidTumorRatioCount
        Int? minDiploidTumorRatioCountCentromere
        #The following should be in the same directory.
        File geneDataCsv
        File proteinFeaturesCsv
        File transExonDataCsv
        File transSpliceDataCsv

        Int threads = 1
        Int timeMinutes = 30
        String memory = "9GiB"
        String javaXmx = "8G"
        # clone of quay.io/biocontainers/hmftools-purple:3.2--hdfd78af_0 with 'ln -s /usr/local/lib/libwebp.so.7 /usr/local/lib/libwebp.so.6'
        #String dockerImage = "quay.io/biowdl/hmftools-purple:3.2" FIXME see if biocontainer works
        String dockerImage = "quay.io/biocontainers/hmftools-purple:3.7.1--hdfd78af_0"
    }

    command {
        PURPLE -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        ~{"-reference " + referenceName} \
        ~{"-germline_vcf " + germlineVcf} \
        ~{"-germline_hotspots " + germlineHotspots} \
        ~{"-germline_del_freq_file " + germlineDelFreqFile} \
        -tumor ~{tumorName} \
        -output_dir ~{outputDir} \
        -amber ~{sub(amberOutput[0], basename(amberOutput[0]), "")} \
        -cobalt ~{sub(cobaltOutput[0], basename(cobaltOutput[0]), "")} \
        -gc_profile ~{gcProfile} \
        -somatic_vcf ~{somaticVcf} \
        -structural_vcf ~{filteredSvVcf} \
        -sv_recovery_vcf ~{fullSvVcf} \
        -circos /usr/local/bin/circos \
        -ref_genome ~{referenceFasta} \
        -ref_genome_version ~{refGenomeVersion} \
        -ensembl_data_dir ~{sub(geneDataCsv, basename(geneDataCsv), "")} \
        -run_drivers \
        -somatic_hotspots ~{somaticHotspots} \
        -driver_gene_panel ~{driverGenePanel} \
        ~{"-highly_diploid_percentage " + highlyDiploidPercentage} \
        ~{"-somatic_min_purity_spread " + somaticMinPuritySpread} \
        ~{"-target_regions_bed " + targetRegionsBed} \
        ~{"-target_regions_ratios " + targetRegionsRatios} \
        ~{"-target_regions_msi_indels " + targetRegionsMsiIndels} \
        ~{"-min_diploid_tumor_ratio_count " + minDiploidTumorRatioCount} \
        ~{"-min_diploid_tumor_ratio_count_centromere" + minDiploidTumorRatioCountCentromere} \
        -threads ~{threads}
    }

    output {
        File driverCatalogGermlineTsv = "~{outputDir}/~{tumorName}.driver.catalog.germline.tsv"
        File driverCatalogSomaticTsv = "~{outputDir}/~{tumorName}.driver.catalog.somatic.tsv"
        File purpleCnvGeneTsv = "~{outputDir}/~{tumorName}.purple.cnv.gene.tsv"
        File purpleCnvSomaticTsv = "~{outputDir}/~{tumorName}.purple.cnv.somatic.tsv"
        File purpleGermlineDeletionTsv = "~{outputDir}/~{tumorName}.purple.germline.deletion.tsv"
        File purpleGermlineVcf = "~{outputDir}/~{tumorName}.purple.germline.vcf.gz"
        File purpleGermlineVcfIndex = "~{outputDir}/~{tumorName}.purple.germline.vcf.gz.tbi"
        File purplePurityRangeTsv = "~{outputDir}/~{tumorName}.purple.purity.range.tsv"
        File purplePurityTsv = "~{outputDir}/~{tumorName}.purple.purity.tsv"
        File purpleQc = "~{outputDir}/~{tumorName}.purple.qc"
        File purpleSegmentTsv = "~{outputDir}/~{tumorName}.purple.segment.tsv"
        File purpleSomaticClonalityTsv = "~{outputDir}/~{tumorName}.purple.somatic.clonality.tsv"
        File purpleSomaticHistTsv = "~{outputDir}/~{tumorName}.purple.somatic.hist.tsv"
        File purpleSomaticVcf = "~{outputDir}/~{tumorName}.purple.somatic.vcf.gz"
        File purpleSomaticVcfIndex = "~{outputDir}/~{tumorName}.purple.somatic.vcf.gz.tbi"
        File purpleSvVcf = "~{outputDir}/~{tumorName}.purple.sv.vcf.gz"
        File purpleSvVcfIndex = "~{outputDir}/~{tumorName}.purple.sv.vcf.gz.tbi"
        File purpleVersion = "~{outputDir}/purple.version"
        File circosPlot = "~{outputDir}/plot/~{tumorName}.circos.png"
        File copynumberPlot = "~{outputDir}/plot/~{tumorName}.copynumber.png"
        File inputPlot = "~{outputDir}/plot/~{tumorName}.input.png"
        File mapPlot = "~{outputDir}/plot/~{tumorName}.map.png"
        File purityRangePlot = "~{outputDir}/plot/~{tumorName}.purity.range.png"
        File segmentPlot = "~{outputDir}/plot/~{tumorName}.segment.png"
        File somaticClonalityPlot = "~{outputDir}/plot/~{tumorName}.somatic.clonality.png"
        File somaticPlot = "~{outputDir}/plot/~{tumorName}.somatic.png"
        File? somaticRainfallPlot = "~{outputDir}/plot/~{tumorName}.somatic.rainfall.png"
        File circosNormalRatio = "~{outputDir}/circos/~{referenceName}.ratio.circos"
        File circosBaf = "~{outputDir}/circos/~{tumorName}.baf.circos"
        File circosConf = "~{outputDir}/circos/~{tumorName}.circos.conf"
        File circosCnv = "~{outputDir}/circos/~{tumorName}.cnv.circos"
        File circosIndel = "~{outputDir}/circos/~{tumorName}.indel.circos"
        File circosInputConf = "~{outputDir}/circos/~{tumorName}.input.conf"
        File circosLink = "~{outputDir}/circos/~{tumorName}.link.circos"
        File circosMap = "~{outputDir}/circos/~{tumorName}.map.circos"
        File circosTumorRatio = "~{outputDir}/circos/~{tumorName}.ratio.circos"
        File circosSnp = "~{outputDir}/circos/~{tumorName}.snp.circos"
        File circosGaps = "~{outputDir}/circos/gaps.txt"
        Array[File] outputs = [driverCatalogSomaticTsv, purpleCnvGeneTsv,
            purpleCnvSomaticTsv, purplePurityRangeTsv, purplePurityTsv, purpleQc,
            purpleSegmentTsv, purpleSomaticClonalityTsv, purpleSomaticHistTsv,
            purpleSomaticVcf, purpleSomaticVcfIndex, purpleSvVcf, purpleSvVcfIndex,
            purpleVersion, purpleGermlineVcf, purpleGermlineVcfIndex, driverCatalogGermlineTsv,
            purpleGermlineDeletionTsv]
        Array[File] plots = select_all([circosPlot, copynumberPlot, inputPlot, mapPlot, purityRangePlot,
            segmentPlot, somaticClonalityPlot, somaticPlot, somaticRainfallPlot])
        Array[File] circos = [circosNormalRatio, circosConf, circosIndel, circosLink,
            circosTumorRatio, circosGaps, circosBaf, circosCnv, circosInputConf, circosMap,
            circosSnp]
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        cpu: threads
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        referenceName: {description: "the name of the normal sample.", category: "required"}
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        outputDir: {description: "The path to the output directory.", category: "common"}
        amberOutput: {description: "The output files of hmftools amber.", category: "required"}
        cobaltOutput: {description: "The output files of hmftools cobalt", category: "required"}
        gcProfile: {description: "A file describing the GC profile of the reference genome.", category: "required"}
        somaticVcf: {description: "The somatic variant calling results.", category: "required"}
        germlineVcf: {description: "The germline variant calling results.", category: "required"}
        filteredSvVcf: {description: "The filtered structural variant calling results.", category: "required"}
        fullSvVcf: {description: "The unfiltered structural variant calling results.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        driverGenePanel: {description: "A TSV file describing the driver gene panel.", category: "required"}
        somaticHotspots: {description: "A vcf file with hotspot somatic variant sites.", category: "required"}
        germlineHotspots: {description: "A vcf file with hotspot germline variant sites.", category: "required"}
        highlyDiploidPercentage: {description: "Equivalent to PURPLE's `-highly_diploid_percentage` option.", category: "advanced"}
        somaticMinPuritySpread: {description: "Equivalent to PURPLE's `-somatic_min_purity_spread` option.", category: "advanced"}
        geneDataCsv: {description: "A  CSV file containing gene information, must be in the same directory as `proteinFeaturesCsv`, `transExonDataCsv` and `transSpliceDataCsv`.", category: "required"}
        proteinFeaturesCsv: {description: "A  CSV file containing protein feature information, must be in the same directory as `geneDataCsv`, `transExonDataCsv` and `transSpliceDataCsv`.", category: "required"}
        transExonDataCsv: {description: "A  CSV file containing transcript exon information, must be in the same directory as `geneDataCsv`, `proteinFeaturesCsv` and `transSpliceDataCsv`.", category: "required"}
        transSpliceDataCsv: {description: "A  CSV file containing transcript splicing information, must be in the same directory as `geneDataCsv`, `proteinFeaturesCsv` and `transExonDataCsv`.", category: "required"}


        threads: {description: "The number of threads the program will use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Rose {
    input {
        File actionabilityDatabaseTsv
        Boolean hg38 = false
        File driverGeneTsv
        File purplePurityTsv
        File purpleQc
        File purpleGeneCopyNumberTsv
        File purpleSomaticDriverCatalogTsv
        File purpleGermlineDriverCatalogTsv
        File purpleSomaticVcf
        File purpleSomaticVcfIndex
        File purpleGermlineVcf
        File purpleGermlineVcfIndex
        File linxFusionTsv
        File linxBreakendTsv
        File linxDriverCatalogTsv
        File annotatedVirusTsv
        File chordPredictionTxt
        File cuppaResultCsv
        String outputDir = "./rose"
        String tumorName
        String referenceName

        String memory = "9GiB"
        String javaXmx = "8G"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/hmftools-rose:1.3--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        rose -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -actionability_database_tsv ~{actionabilityDatabaseTsv} \
        -ref_genome_version ~{if hg38 then "38" else "37"} \
        -driver_gene_tsv ~{driverGeneTsv} \
        -purple_purity_tsv ~{purplePurityTsv} \
        -purple_qc_file ~{purpleQc} \
        -purple_gene_copy_number_tsv ~{purpleGeneCopyNumberTsv} \
        -purple_somatic_driver_catalog_tsv ~{purpleSomaticDriverCatalogTsv} \
        -purple_germline_driver_catalog_tsv ~{purpleGermlineDriverCatalogTsv} \
        -purple_somatic_variant_vcf ~{purpleSomaticVcf} \
        -purple_germline_variant_vcf ~{purpleGermlineVcf} \
        -linx_fusion_tsv ~{linxFusionTsv} \
        -linx_breakend_tsv ~{linxBreakendTsv} \
        -linx_driver_catalog_tsv ~{linxDriverCatalogTsv} \
        -annotated_virus_tsv ~{annotatedVirusTsv} \
        -chord_prediction_txt ~{chordPredictionTxt} \
        -cuppa_result_csv ~{cuppaResultCsv} \
        -output_dir ~{outputDir} \
        -tumor_sample_id ~{tumorName} \
        -ref_sample_id ~{referenceName} \
        -patient_id not_used_because_primary_tumor_tsv_has_only_headers
    }

    output {
        File roseTsv = "~{outputDir}/~{tumorName}.rose.tsv"
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        memory: memory
    }

    parameter_meta {

    }
}

task Sage {
    input {
        Array[String]+ tumorName
        Array[File]+ tumorBam
        Array[File]+ tumorBamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        File hotspots
        File panelBed
        File highConfidenceBed
        File coverageBed
        Boolean hg38 = false
        Boolean panelOnly = false
        String outputPath = "./sage.vcf.gz"
        #The following should be in the same directory.
        File geneDataCsv
        File proteinFeaturesCsv
        File transExonDataCsv
        File transSpliceDataCsv

        Array[String] referenceName = []
        Array[File] referenceBam = []
        Array[File] referenceBamIndex = []
        Int? hotspotMinTumorQual
        Int? panelMinTumorQual
        Int? hotspotMaxGermlineVaf
        Int? hotspotMaxGermlineRelRawBaseQual
        Int? panelMaxGermlineVaf
        Int? panelMaxGermlineRelRawBaseQual
        Int? refSampleCount
        Float? hotspotMinTumorVaf
        Int? highConfidenceMinTumorQual
        Int? lowConfidenceMinTumorQual

        Int threads = 32
        String javaXmx = "16G"
        String memory = "20GiB"
        Int timeMinutes = 720
        String dockerImage = "quay.io/biocontainers/hmftools-sage:3.2.3--hdfd78af_0"
    }

    command {
        SAGE -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -tumor ~{sep="," tumorName} \
        -tumor_bam ~{sep="," tumorBam} \
        ~{if length(referenceName) > 0 then "-reference" else ""} ~{sep="," referenceName} \
        ~{if length(referenceBam) > 0 then "-reference_bam" else ""}  ~{sep="," referenceBam} \
        -hotspots ~{hotspots} \
        ~{"-hotspot_min_tumor_qual " + hotspotMinTumorQual} \
        -high_confidence_bed ~{highConfidenceBed} \
        -panel_bed ~{panelBed} \
        -coverage_bed ~{coverageBed} \
        -ref_genome ~{referenceFasta} \
        -ref_genome_version ~{true="38" false="37" hg38} \
        -ensembl_data_dir ~{sub(geneDataCsv, basename(geneDataCsv), "")} \
        -write_bqr_data \
        -write_bqr_plot \
        -out ~{outputPath} \
        -threads ~{threads} \
        ~{"-panel_min_tumor_qual " + panelMinTumorQual} \
        ~{"-hotspot_max_germline_vaf " + hotspotMaxGermlineVaf} \
        ~{"-hotspot_max_germline_rel_raw_base_qual " + hotspotMaxGermlineRelRawBaseQual} \
        ~{"-panel_max_germline_vaf " + panelMaxGermlineVaf} \
        ~{"-panel_max_germline_rel_raw_base_qual " + panelMaxGermlineRelRawBaseQual} \
        ~{true="-panel_only" false="" panelOnly} \
        ~{"-ref_sample_count " + refSampleCount} \
        ~{"-hotspot_min_tumor_vaf " + hotspotMinTumorVaf} \
        ~{"-high_confidence_min_tumor_qual " + highConfidenceMinTumorQual} \
        ~{"-low_confidence_min_tumor_qual " + lowConfidenceMinTumorQual}
    }

    String outputDir = sub(outputPath, basename(outputPath), "")

    output { #FIXME does it produce multiple plots/tsvs if multiple samples are given?
        File outputVcf = outputPath
        File outputVcfIndex = outputPath + ".tbi"
        File? referenceSageBqrPng = "~{outputDir}/~{referenceName[0]}.sage.bqr.png"
        File? referenceSageBqrTsv = "~{outputDir}/~{referenceName[0]}.sage.bqr.tsv"
        File tumorSageBqrPng = "~{outputDir}/~{tumorName[0]}.sage.bqr.png"
        File tumorSageBqrTsv = "~{outputDir}/~{tumorName[0]}.sage.bqr.tsv"
        File sageGeneCoverageTsv = "~{outputDir}/~{tumorName[0]}.sage.gene.coverage.tsv"
        File referenceSageExonMediansTsv = "~{outputDir}/~{tumorName[0]}.sage.exon.medians.tsv"
        Array[File] outputs = select_all([outputVcf, outputVcfIndex, referenceSageBqrPng,
                                          referenceSageBqrTsv, tumorSageBqrPng, tumorSageBqrTsv,
                                          sageGeneCoverageTsv, referenceSageExonMediansTsv])
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        cpu: threads
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        tumorBam: {description: "The BAM file for the tumor sample.", category: "required"}
        tumorBamIndex: {description: "The index of the BAM file for the tumor sample.", category: "required"}
        referenceName: {description: "The name of the normal/reference sample.", category: "common"}
        referenceBam: {description: "The BAM file for the normal sample.", category: "common"}
        referenceBamIndex: {description: "The index of the BAM file for the normal sample.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        hotspots: {description: "A vcf file with hotspot variant sites.", category: "required"}
        panelBed: {description: "A bed file describing coding regions to search for in frame indels.", category: "required"}
        highConfidenceBed: {description: "A bed files describing high confidence mapping regions.", category: "required"}
        hotspotMinTumorQual: {description: "Equivalent to sage's `hotspot_min_tumor_qual` option.", category: "advanced"}
        panelMinTumorQual: {description: "Equivalent to sage's `panel_min_tumor_qual` option.", category: "advanced"}
        hotspotMaxGermlineVaf: {description: "Equivalent to sage's `hotspot_max_germline_vaf` option.", category: "advanced"}
        hotspotMaxGermlineRelRawBaseQual: {description: "Equivalent to sage's `hotspot_max_germline_rel_raw_base_qual` option.", category: "advanced"}
        panelMaxGermlineVaf: {description: "Equivalent to sage's `panel_max_germline_vaf` option.", category: "advanced"}
        panelMaxGermlineRelRawBaseQual: {description: "Equivalent to sage's `panel_max_germline_vaf` option.", category: "advanced"}
        refSampleCount: {description: "Equivalent to sage's `ref_sample_count` option.", category: "advanced"}
        hg38: {description: "Whether or not the refernce genome is HG18, if false HG19 is assumed.", category: "common"}

        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task SageAppend {
    input {
        String sampleName
        File bamFile
        File? bamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        File sageVcf
        String outPath = "./sage_append.vcf"

        Int threads = 2
        String javaXmx = "32G"
        String memory = "33GiB"
        Int timeMinutes = 720
        String dockerImage = "quay.io/biocontainers/hmftools-sage:3.2.3--hdfd78af_0"
    }

    command {
        java -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -cp /usr/local/share/hmftools-sage-3.2.3-0/sage.jar \
        com.hartwig.hmftools.sage.append.SageAppendApplication \
        -reference ~{sampleName} \
        -reference_bam ~{bamFile} \
        -ref_genome ~{referenceFasta} \
        -input_vcf ~{sageVcf} \
        -out ~{outPath} \
        -threads ~{threads}
    }

    output {
        File vcf = outPath
        File index = "~{outPath}.tbi"
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        cpu: threads
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        sampleName: {description: "The sample id.", category: "required"}
        bamFile: {description: "The input BAM file.", category: "required"}
        bamIndex: {description: "Index for the input BAM file", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        sageVcf: {description: "A VCF file from Sage or Purple.", category: "required"}
        outPath: {description: "Location to write the output to.", category: "required"}

        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Sigs {
    input {
        String sampleName
        File signaturesFile
        File somaticVcfFile
        File somaticVcfIndex
        String outputDir = "./sigs"

        String javaXmx = "4G"
        String memory = "5GiB"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/hmftools-sigs:1.1--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        sigs -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -sample ~{sampleName} \
        -signatures_file ~{signaturesFile} \
        -somatic_vcf_file ~{somaticVcfFile} \
        -output_dir ~{outputDir}
    }

    output {
        File sigAllocationTsv = "~{outputDir}/~{sampleName}.sig.allocation.tsv"
        File sigSnvCountsCsv = "~{outputDir}/~{sampleName}.sig.snv_counts.csv"
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        memory: memory
    }

    parameter_meta {

    }
}

task SvPrep {
    # for ref also add tumorJunctionFile
    input {
        String sampleName
        File bamFile
        File bamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        File blacklistBed
        File knownFusionBed
        String outputDir = "."

        File? existingJunctionFile
        Boolean hg38 = false

        Int threads = 10
        String javaXmx = "48G"
        String memory = "50GiB"
        Int timeMinutes = 120
        String dockerImage = "quay.io/biocontainers/hmftools-sv-prep:1.1--hdfd78af_1"
    }

    command {
        set -e
        SvPrep -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -sample ~{sampleName} \
        -bam_file ~{bamFile} \
        -ref_genome ~{referenceFasta} \
        -ref_genome_version ~{true="38" false="37" hg38} \
        -blacklist_bed ~{blacklistBed} \
        -known_fusion_bed ~{knownFusionBed} \
        ~{"-existing_junction_file " + existingJunctionFile} \
        -write_types "JUNCTIONS;BAM;FRAGMENT_LENGTH_DIST" \
        -output_dir ~{outputDir} \
        -threads ~{threads}
        samtools sort -O bam ~{outputDir}/~{sampleName}.sv_prep.bam -o ~{outputDir}/~{sampleName}.sv_prep.sorted.bam
        samtools index ~{outputDir}/~{sampleName}.sv_prep.sorted.bam
    }

    output {
        File preppedBam = "~{outputDir}/~{sampleName}.sv_prep.sorted.bam"
        File preppedBamIndex = "~{outputDir}/~{sampleName}.sv_prep.sorted.bam.bai"
        File junctions = "~{outputDir}/~{sampleName}.sv_prep.junctions.csv"
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        cpu: threads
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        sampleName: {description: "The name of the sample.", category: "required"}
        bamFile: {description: "The BAM file to prepare for SV calling with GRIDSS.", category: "required"}
        bamIndex: {description: "The index for the BAM file.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        blacklistBed: {description: "Blacklist bed file.", category: "required"}
        knownFusionBed: {description: "Bed file with known fusion sites", category: "required"}
        outputDir: {description: "Path to the output directory.", category: "common"}
        existingJunctionFile: {description: "Junctions file generated by an earlier run of this tool, eg. from a paired sample.", category: "common"}
        hg38: {description: "Whether or not the refernce genome is HG18, if false HG19 is assumed.", category: "common"}

        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task SvPrepDepthAnnotator {
    input {
        File inputVcf
        File inputVcfIndex
        Array[File]+ bamFiles
        Array[File]+ bamIndexes
        Array[String]+ samples
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        Boolean hg38 = false
        String outputVcf = "gridss.depth_annotated.vcf.gz"

        Int threads = 10
        String javaXmx = "48G"
        String memory = "50GiB"
        Int timeMinutes = 240
        String dockerImage = "quay.io/biocontainers/hmftools-sv-prep:1.1--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputVcf})"
        java -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -cp /usr/local/share/hmftools-sv-prep-1.1-0/sv-prep.jar \
        com.hartwig.hmftools.svprep.depth.DepthAnnotator \
        -input_vcf ~{inputVcf} \
        -output_vcf ~{outputVcf} \
        -samples ~{sep="," samples} \
        -bam_files ~{sep="," bamFiles} \
        -ref_genome ~{referenceFasta} \
        -ref_genome_version ~{if hg38 then "38" else "37"} \
        -threads ~{threads}
    }

    output {
        File vcf = outputVcf
        File vcfIndex = outputVcf + ".tbi"
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        cpu: threads
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        samples: {description: "The names of the samples.", category: "required"}
        bamFiles: {description: "The BAM files.", category: "required"}
        bamIndexes: {description: "The indexes for the BAM files.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        hg38: {description: "Whether or not the refernce genome is HG18, if false HG19 is assumed.", category: "common"}
        outputVcf: {description: "The path for the output VCF.", category: "common"}

        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task VirusInterpreter {
    input {
        String sampleId
        File purplePurityTsv
        File prupleQcFile
        File tumorSampleWgsMetricsFile
        File virusBreakendTsv
        File taxonomyDbTsv
        File virusReportingDbTsv
        String outputDir = "."

        String memory = "3GiB"
        String javaXmx = "2G"
        Int timeMinutes = 15
        String dockerImage = "quay.io/biowdl/virus-interpreter:1.2"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        virus-interpreter -Xmx~{javaXmx} -XX:ParallelGCThreads=1  \
        -sample_id ~{sampleId} \
        -purple_purity_tsv ~{purplePurityTsv} \
        -purple_qc_file ~{prupleQcFile} \
        -tumor_sample_wgs_metrics_file ~{tumorSampleWgsMetricsFile} \
        -virus_breakend_tsv ~{virusBreakendTsv} \
        -taxonomy_db_tsv ~{taxonomyDbTsv} \
        -virus_reporting_db_tsv ~{virusReportingDbTsv} \
        -output_dir ~{outputDir}
    }

    output {
        File virusAnnotatedTsv = "~{outputDir}/~{sampleId}.virus.annotated.tsv"
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        sampleId: {description: "The name of the sample.", category: "required"}
        purplePurityTsv: {description: "The purity file produced by purple.", category: "required"}
        prupleQcFile: {description: "The QC file produced by purple.", category: "required"}
        tumorSampleWgsMetricsFile: {description: "The picard WGS metrics file for this sample.", category: "required"}
        virusBreakendTsv: {description: "The TSV output from virusbreakend.", category: "required"}
        taxonomyDbTsv: {description: "A taxonomy database tsv.", category: "required"}
        virusReportingDbTsv: {description: "A virus reporting tsv.", category: "required"}
        outputDir: {description: "The directory the output will be written to.", category: "required"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
