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
        String normalName
        File normalBam
        File normalBamIndex
        String tumorName
        File tumorBam
        File tumorBamIndex
        String outputDir = "./amber"
        File loci
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict

        Int threads = 2
        String memory = "33G"
        String javaXmx = "32G"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/hmftools-amber:3.5--0"
    }

    command {
        AMBER -Xmx~{javaXmx} \
        -reference ~{normalName} \
        -reference_bam ~{normalBam} \
        -tumor ~{tumorName} \
        -tumor_bam ~{tumorBam} \
        -output_dir ~{outputDir} \
        -threads ~{threads} \
        -ref_genome ~{referenceFasta} \
        -loci ~{loci}
    }

    output {
        File version = "~{outputDir}/amber.version"
        File tumorBafPcf = "~{outputDir}/~{tumorName}.amber.baf.pcf"
        File tumorBafTsv = "~{outputDir}/~{tumorName}.amber.baf.tsv"
        File tumorBafVcf = "~{outputDir}/~{tumorName}.amber.baf.vcf.gz"
        File tumorBafVcfIndex = "~{outputDir}/~{tumorName}.amber.baf.vcf.gz.tbi"
        File tumorContaminationVcf = "~{outputDir}/~{tumorName}.amber.contamination.vcf.gz"
        File tumorContaminationVcfIndex = "~{outputDir}/~{tumorName}.amber.contamination.vcf.gz.tbi"
        File tumorContaminationTsv = "~{outputDir}/~{tumorName}.amber.contamination.tsv"
        File tumorQc = "~{outputDir}/~{tumorName}.amber.qc"
        File normalSnpVcf = "~{outputDir}/~{normalName}.amber.snp.vcf.gz"
        File normalSnpVcfIndex = "~{outputDir}/~{normalName}.amber.snp.vcf.gz.tbi"
        Array[File] outputs = [version, tumorBafPcf, tumorBafTsv, tumorBafVcf, tumorBafVcfIndex, 
            tumorContaminationVcf, tumorContaminationVcfIndex, tumorContaminationTsv, tumorQc, 
            normalSnpVcf, normalSnpVcfIndex]
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        cpu: threads
    }

    parameter_meta {
        normalName: {description: "the name of the normal sample.", category: "required"}
        normalBam: {description: "The normal BAM file.", category: "required"}
        normalBamIndex: {description: "The index for the normal BAM file.", category: "required"}
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        tumorBam: {description: "The tumor BAM file.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor BAM file.", category: "required"}
        outputDir: {description: "The path to the output directory.", category: "common"}
        loci: {description: "A VCF file containing likely heterozygous sites.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
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
        String normalName
        File normalBam
        File normalBamIndex
        String tumorName
        File tumorBam
        File tumorBamIndex
        String outputDir = "./cobalt"
        File gcProfile
        
        Int threads = 1
        String memory = "9G"
        String javaXmx = "8G"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/hmftools-cobalt:1.10--0"
    }

    command {
        COBALT -Xmx~{javaXmx} \
        -reference ~{normalName} \
        -reference_bam ~{normalBam} \
        -tumor ~{tumorName} \
        -tumor_bam ~{tumorBam} \
        -output_dir ~{outputDir} \
        -threads ~{threads} \
        -gc_profile ~{gcProfile}
    }

    output {
        File version = "~{outputDir}/cobalt.version"
        File normalGcMedianTsv = "~{outputDir}/~{normalName}.cobalt.gc.median.tsv"
        File normalRationMedianTsv = "~{outputDir}/~{normalName}.cobalt.ratio.median.tsv"
        File normalRationPcf = "~{outputDir}/~{normalName}.cobalt.ratio.pcf"
        File tumorGcMedianTsv = "~{outputDir}/~{tumorName}.cobalt.gc.median.tsv"
        File tumorRatioPcf = "~{outputDir}/~{tumorName}.cobalt.ratio.pcf"
        File tumorRatioTsv = "~{outputDir}/~{tumorName}.cobalt.ratio.tsv"
        File tumorChrLen = "~{outputDir}/~{tumorName}.chr.len"
        Array[File] outputs = [version, normalGcMedianTsv, normalRationMedianTsv,
            normalRationPcf, tumorGcMedianTsv, tumorRatioPcf, tumorRatioTsv, tumorChrLen]
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
        cpu: threads
    }

    parameter_meta {
        normalName: {description: "the name of the normal sample.", category: "required"}
        normalBam: {description: "The normal BAM file.", category: "required"}
        normalBamIndex: {description: "The index for the normal BAM file.", category: "required"}
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        tumorBam: {description: "The tumor BAM file.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor BAM file.", category: "required"}
        outputDir: {description: "The path to the output directory.", category: "common"}
        gcProfile: {description: "A file describing the GC profile of the reference genome.", category: "required"}
        threads: {description: "The number of threads the program will use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task GripssApplicationKt {
    input {
        File inputVcf
        String outputPath = "gripss.vcf.gz"
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        File breakpointHotspot
        File breakendPon
        File breakpointPon

        String memory = "25G"
        String javaXmx = "24G"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/hmftools-gripss:1.8--0"
    }

    command {
        java -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -cp /usr/local/share/hmftools-gripss-1.8-0/gripss.jar \
        com.hartwig.hmftools.gripss.GripssApplicationKt \
        -ref_genome ~{referenceFasta} \
        -breakpoint_hotspot ~{breakpointHotspot} \
        -breakend_pon ~{breakendPon} \
        -breakpoint_pon ~{breakpointPon} \
        -input_vcf ~{inputVcf} \
        -output_vcf ~{outputPath} 
    }

    output {
        File outputVcf = outputPath
        File outputVcfIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        inputVcf: {description: "The input VCF.", category: "required"}
        outputPath: {description: "The path where th eoutput VCF will be written.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        breakpointHotspot: {description: "Equivalent to the `-breakpoint_hotspot` option.", category: "required"}
        breakendPon: {description: "Equivalent to the `-breakend_pon` option.", category: "required"}
        breakpointPon: {description: "Equivalent to the `breakpoint_pon` option.", category: "required"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task GripssHardFilterApplicationKt {
    input {
        File inputVcf
        String outputPath = "gripss_hard_filter.vcf.gz"

        String memory = "25G"
        String javaXmx = "24G"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/hmftools-gripss:1.8--0"
    }

    command {
        java -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -cp /usr/local/share/hmftools-gripss-1.8-0/gripss.jar \
        com.hartwig.hmftools.gripss.GripssHardFilterApplicationKt \
        -input_vcf ~{inputVcf} \
        -output_vcf ~{outputPath} 
    }

    output {
        File outputVcf = outputPath
        File outputVcfIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        inputVcf: {description: "The input VCF.", category: "required"}
        outputPath: {description: "The path where th eoutput VCF will be written.", category: "common"}
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
        String normalName
        String tumorName

        String javaXmx = "10G"
    }

    command {
        java -Xmx10G \
        -jar /opt/tools/health-checker/3.1/health-checker.jar \
        -reference ~{normalName} \
        -tumor ~{tumorName} \
        -metrics_dir ~{metricsPath} \
        -amber_dir ~{sub(amberOutput[0], basename(amberOutput[0]), "")} \
        -purple_dir ~{sub(purpleOutput[0], basename(purpleOutput[0]), "")} \
        -output_dir ~{outputDir}
    }

    #    super("health-checker",
    #             Versions.HEALTH_CHECKER,
    #             "health-checker.jar",
    #             "10G",
    #             Lists.newArrayList("-reference",
    #                     referenceSampleName,
    #                     "-tumor",
    #                     tumorSampleName,
    #                     "-ref_wgs_metrics_file",
    #                     referenceMetricsPath,
    #                     "-tum_wgs_metrics_file",
    #                     tumorMetricsPath,
    #                     "-ref_flagstat_file",
    #                     referenceFlagstatPath,
    #                     "-tum_flagstat_file",
    #                     tumorFlagstatPath,
    #                     "-purple_dir",
    #                     purplePath,
    #                     "-output_dir",
    #                     outputPath));    

    output {

    }


}

task Purple {
    input {
        String normalName
        String tumorName
        String outputDir = "./purple"
        Array[File]+ amberOutput
        Array[File]+ cobaltOutput
        File gcProfile
        File somaticVcf
        File filteredSvVcf
        File fullSvVcf
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        File driverGenePanel
        File hotspots
        
        Int threads = 1
        Int timeMinutes = 60
        String memory = "13G"
        String javaXmx = "12G"
        String dockerImage = "quay.io/biocontainers/hmftools-purple:2.51--1"
    }

    command {
        PURPLE -Xmx~{javaXmx} \
        -reference ~{normalName} \
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
        -driver_catalog \
        -driver_gene_panel ~{driverGenePanel} \
        -hotspots ~{hotspots} \
        -threads ~{threads}

        # TODO if shallow also the following:
        #-highly_diploid_percentage 0.88 \
        #-somatic_min_total 100 \
        #-somatic_min_purity_spread 0.1
    }

    output {
        #TODO
    }

    runtime {
        time_minutes: timeMinutes # !UnknownRuntimeKey
        cpu: threads
        docker: dockerImage
        memory: memory
    }

    parameter_meta {
        normalName: {description: "the name of the normal sample.", category: "required"}
        tumorName: {description: "The name of the tumor sample.", category: "required"}
        outputDir: {description: "The path to the output directory.", category: "common"}
        amberOutput: {description: "The output files of hmftools amber.", category: "required"}
        cobaltOutput: {description: "The output files of hmftools cobalt", category: "required"}
        gcProfile: {description: "A file describing the GC profile of the reference genome.", category: "required"}
        somaticVcf: {description: "The somatic variant calling results.", category: "required"}
        filteredSvVcf: {description: "The filtered structural variant calling results.", category: "required"}
        fullSvVcf: {description: "The unfiltered structural variant calling results.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        driverGenePanel: {description: "A bed file describing the driver gene panel.", category: "required"}
        hotspots: {description: "A vcf file with hotspot variant sites.", category: "required"}

        threads: {description: "The number of threads the program will use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Sage {
    input {
        String tumorName
        File tumorBam
        File tumorBamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        File hotspots
        File panelBed
        File highConfidenceBed
        Boolean hg38 = false
        String outputPath = "./sage.vcf.gz"

        String? normalName
        File? normalBam
        File? normalBamIndex

        Int threads = 2
        String javaXmx = "32G"
        String memory = "33G"
        Int timeMinutes = 1 + ceil(size(select_all([tumorBam, normalBam]), "G") * 10 / threads) #FIXME make sure this is enough
        String dockerImage = "quay.io/biocontainers/hmftools-sage:2.2--2"
    }

    command {
        java -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -cp /usr/local/share/hmftools-sage-2.2-2/sage.jar \
        com.hartwig.hmftools.sage.SageApplication \
        -tumor ~{tumorName} \
        -tumor_bam ~{tumorBam} \
        ~{"-reference " + normalName} \
        ~{"-reference_bam " + normalBam} \
        -ref_genome ~{referenceFasta} \
        -hotspots ~{hotspots} \
        -panel_bed ~{panelBed} \
        -high_confidence_bed ~{highConfidenceBed} \
        -assembly ~{true="hg38" false="hg19" hg38} \
        -threads ~{threads} \
        -out ~{outputPath}
    }

    output {
        File outputVcf = outputPath
        File outputVcfIndex = outputPath + ".tbi"
        # There is some plots as well, but in the current container the labels in the plots are just series of `□`s.
        # This seems to be a systemic issue with R generated plots in biocontainers...
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
        normalName: {description: "The name of the normal/reference sample.", category: "common"}
        normalBam: {description: "The BAM file for the normal sample.", category: "common"}
        normalBamIndex: {description: "The index of the BAM file for the normal sample.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.",
                             category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        hotspots: {description: "A vcf file with hotspot variant sites.", category: "required"}
        panelBed: {description: "A bed file describing coding regions to search for in frame indels.", category: "required"}
        highConfidenceBed: {description: "A bed files describing high confidence mapping regions.", category: "required"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
