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
        Int timeMinutes = 1200
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
        Int timeMinutes = 1200
        String dockerImage = "quay.io/biocontainers/hmftools-cobalt:1.11--0"
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
        String tumorName
        String normalName
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        File breakpointHotspot
        File breakendPon
        File breakpointPon

        String memory = "25G"
        String javaXmx = "24G"
        Int timeMinutes = 120
        String dockerImage = "quay.io/biocontainers/hmftools-gripss:1.9--0"
    }

    command {
        java -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -cp /usr/local/share/hmftools-gripss-1.9-0/gripss.jar \
        com.hartwig.hmftools.gripss.GripssApplicationKt \
        -tumor ~{tumorName} \
        ~reference ~{normalName} \
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
        Int timeMinutes = 120
        String dockerImage = "quay.io/biocontainers/hmftools-gripss:1.9--0"
    }

    command {
        java -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -cp /usr/local/share/hmftools-gripss-1.9-0/gripss.jar \
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
    # WIP
    input {
        String outputDir = "."
        String normalName
        File normalFlagstats
        File normalMetrics
        String tumorName
        File tumorFlagstats
        File tumorMetrics
        Array[File]+ purpleOutput

        String javaXmx = "10G"
        String memory = "11G"
        Int timeMinutes = 10
        String dockerImage = "quay.io/biowdl/health-checker:3.2"
    }

    command {
        set -e
        mkdir -p ~{outputDir}
        health-checker -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -reference ~{normalName} \
        -ref_flagstat_file ~{normalFlagstats} \
        -ref_wgs_metrics_file ~{normalMetrics} \
        -tumor ~{tumorName} \
        -tum_flagstat_file ~{tumorFlagstats} \
        -tum_wgs_metrics_file ~{tumorMetrics} \
        -purple_dir ~{sub(purpleOutput[0], basename(purpleOutput[0]), "")} \
        -output_dir ~{outputDir}   
    }
 

    output {
        File? healthCheckSucceeded = "~{outputDir}/~{tumorName}.HealthCheckSucceeded"
        File? healthCheckFailed = "~{outputDir}/~{tumorName}.HealthCheckFailed"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        outputDir: {description: "The path the output will be written to.", category:"required"}
        normalName: {description: "The name of the normal sample.", category: "required"}
        normalFlagstats: {description: "The flagstats for the normal sample.", category: "required"}
        normalMetrics: {description: "The picard WGS metrics for the normal sample.", category: "required"}
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

task Linx {
    input {
        String sampleName
        File svVcf
        File svVcfIndex
        Array[File]+ purpleOutput
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String refGenomeVersion
        String outputDir = "./linx"
        File fragileSiteCsv
        File lineElementCsv
        File replicationOriginsBed
        File viralHostsCsv
        File knownFusionCsv
        File driverGenePanel
        #The following should be in the same directory.
        File geneDataCsv
        File proteinFeaturesCsv
        File transExonDataCsv
        File transSpliceDataCsv

        String memory = "9G"
        String javaXmx = "8G"
        Int timeMinutes = 30
        String dockerImage = "quay.io/biocontainers/hmftools-linx:1.13--0"
    }

    command {
        linx -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -sample ~{sampleName} \
        -sv_vcf ~{svVcf} \
        -purple_dir ~{sub(purpleOutput[0], basename(purpleOutput[0]), "")} \
        -ref_genome ~{referenceFasta} \
        -ref_genome_version ~{refGenomeVersion} \
        -output_dir ~{outputDir} \
        -fragile_site_file ~{fragileSiteCsv} \
        -line_element_file ~{lineElementCsv} \
        -replication_origins_file ~{replicationOriginsBed} \
        -viral_hosts_file ~{viralHostsCsv} \
        -gene_transcripts_dir ~{sub(geneDataCsv, basename(geneDataCsv), "")} \
        -check_fusions \
        -known_fusion_file ~{knownFusionCsv} \
        -check_drivers \
        -driver_gene_panel ~{driverGenePanel} \
        -chaining_sv_limit 0 \
        -write_vis_data
    }

    output {
        File driverCatalog = "~{outputDir}/~{sampleName}.driver.catalog.tsv"
        File linxBreakend = "~{outputDir}/~{sampleName}.linx.breakend.tsv"
        File linxClusters = "~{outputDir}/~{sampleName}.linx.clusters.tsv"
        File linxDrivers = "~{outputDir}/~{sampleName}.linx.drivers.tsv"
        File linxFusion = "~{outputDir}/~{sampleName}.linx.fusion.tsv"
        File linxLinks = "~{outputDir}/~{sampleName}.linx.links.tsv"
        File linxSvs = "~{outputDir}/~{sampleName}.linx.svs.tsv"
        File linxViralInserts = "~{outputDir}/~{sampleName}.linx.viral_inserts.tsv"
        File linxVisCopyNumber = "~{outputDir}/~{sampleName}.linx.vis_copy_number.tsv"
        File linxVisFusion = "~{outputDir}/~{sampleName}.linx.vis_fusion.tsv"
        File linxVisGeneExon = "~{outputDir}/~{sampleName}.linx.vis_gene_exon.tsv"
        File linxVisProteinDomain = "~{outputDir}/~{sampleName}.linx.vis_protein_domain.tsv"
        File linxVisSegments = "~{outputDir}/~{sampleName}.linx.vis_segments.tsv"
        File linxVisSvData = "~{outputDir}/~{sampleName}.linx.vis_sv_data.tsv"
        File linxVersion = "~{outputDir}/linx.version"
        Array[File] outputs = [driverCatalog, linxBreakend, linxClusters, linxDrivers, linxFusion,
                               linxLinks, linxSvs, linxViralInserts, linxVisCopyNumber,
                               linxVisFusion, linxVisGeneExon, linxVisProteinDomain,
                               linxVisSegments, linxVisSvData, linxVersion]
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
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        refGenomeVersion: {description: "The version of the genome assembly used for alignment. Either \"HG19\" or \"HG38\".", category: "required"}
        outputDir: {description: "The directory the outputs will be written to.", category: "required"}
        fragileSiteCsv: {description: "A list of known fragile sites.", category: "required"}
        lineElementCsv: {description: "A list of known LINE source regions.", category: "required"}
        replicationOriginsBed: {description: "Replication timing input in BED format with replication timing as the 4th column.", category: "required"}
        viralHostsCsv: {description: "A list of the viruses which were used for annotation of the GRIDSS results.", category: "required"}
        knownFusionCsv: {description: "A CSV file describing known fusions.", category: "required"}
        driverGenePanel: {description: "A TSV file describing the driver gene panel.", category: "required"}
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
        File fullSvVcfIndex
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        File driverGenePanel
        File somaticHotspots
        
        Int threads = 1
        Int timeMinutes = 60
        String memory = "13G"
        String javaXmx = "12G"
        String dockerImage = "quay.io/biocontainers/hmftools-purple:2.52--0"
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
        -somatic_hotspots ~{somaticHotspots} \
        -threads ~{threads}
    }

    output {
        File driverCatalogSomaticTsv = "~{outputDir}/~{tumorName}.driver.catalog.somatic.tsv"
        File purpleCnvGeneTsv = "~{outputDir}/~{tumorName}.purple.cnv.gene.tsv"
        File purpleCnvGermlineTsv = "~{outputDir}/~{tumorName}.purple.cnv.germline.tsv"
        File purpleCnvSomaticTsv = "~{outputDir}/~{tumorName}.purple.cnv.somatic.tsv"
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
        File circosPlot = "~{outputDir}/plot/~{tumorName}.circos.png"
        File copynumberPlot = "~{outputDir}/plot/~{tumorName}.copynumber.png"
        File inputPlot = "~{outputDir}/plot/~{tumorName}.input.png"
        File mapPlot = "~{outputDir}/plot/~{tumorName}.map.png"
        File purityRangePlot = "~{outputDir}/plot/~{tumorName}.purity.range.png"
        File segmentPlot = "~{outputDir}/plot/~{tumorName}.segment.png"
        File somaticClonalityPlot = "~{outputDir}/plot/~{tumorName}.somatic.clonality.png"
        File somaticPlot = "~{outputDir}/plot/~{tumorName}.somatic.png"
        File somaticRainfallPlot = "~{outputDir}/plot/~{tumorName}.somatic.rainfall.png"
        File purpleVersion = "~{outputDir}/purple.version"
        File circosNormalRatio = "~{outputDir}/circos/~{normalName}.ratio.circos"
        File circosConf = "~{outputDir}/circos/~{tumorName}.circos.conf"
        File circosIndel = "~{outputDir}/circos/~{tumorName}.indel.circos"
        File circosLink = "~{outputDir}/circos/~{tumorName}.link.circos"
        File circosTumorRatio = "~{outputDir}/circos/~{tumorName}.ratio.circos"
        File circosGaps = "~{outputDir}/circos/gaps.txt"
        File circosBaf = "~{outputDir}/circos/~{tumorName}.baf.circos"
        File circosCnv = "~{outputDir}/circos/~{tumorName}.cnv.circos"
        File circosInputConf = "~{outputDir}/circos/~{tumorName}.input.conf"
        File circosMap = "~{outputDir}/circos/~{tumorName}.map.circos"
        File circosSnp = "~{outputDir}/circos/~{tumorName}.snp.circos"
        Array[File] outputs = [driverCatalogSomaticTsv, purpleCnvGeneTsv, purpleCnvGermlineTsv,
            purpleCnvSomaticTsv, purplePurityRangeTsv, purplePurityTsv, purpleQc, 
            purpleSegmentTsv, purpleSomaticClonalityTsv, purpleSomaticHistTsv, 
            purpleSomaticVcf, purpleSomaticVcfIndex, purpleSvVcf, purpleSvVcfIndex,
            purpleVersion]
        Array[File] plots = [circosPlot, copynumberPlot, inputPlot, mapPlot, purityRangePlot,
            segmentPlot, somaticClonalityPlot, somaticPlot, somaticRainfallPlot]
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
        driverGenePanel: {description: "A TSV file describing the driver gene panel.", category: "required"}
        somaticHotspots: {description: "A vcf file with hotspot variant sites.", category: "required"}

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
        Boolean panelOnly = false
        String outputPath = "./sage.vcf.gz"

        String? normalName
        File? normalBam
        File? normalBamIndex
        Int? hotspotMinTumorQual
        Int? panelMinTumorQual
        Int? hotspotMaxGermlineVaf
        Int? hotspotMaxGermlineRelRawBaseQual
        Int? panelMaxGermlineVaf
        Int? panelMaxGermlineRelRawBaseQual
        String? mnvFilterEnabled
        File? coverageBed

        Int threads = 2
        String javaXmx = "32G"
        String memory = "33G"
        Int timeMinutes = 1 + ceil(size(select_all([tumorBam, normalBam]), "G") * 10 / threads) #FIXME make sure this is enough
        String dockerImage = "quay.io/biocontainers/hmftools-sage:2.6--0"
    }

    command {
        SAGE -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -tumor ~{tumorName} \
        -tumor_bam ~{tumorBam} \
        ~{"-reference " + normalName} \
        ~{"-reference_bam " + normalBam} \
        -ref_genome ~{referenceFasta} \
        -hotspots ~{hotspots} \
        -panel_bed ~{panelBed} \
        -high_confidence_bed ~{highConfidenceBed} \
        -assembly ~{true="hg38" false="hg19" hg38} \
        ~{"-hotspot_min_tumor_qual " + hotspotMinTumorQual} \
        ~{"-panel_min_tumor_qual " + panelMinTumorQual} \
        ~{"-hotspot_max_germline_vaf " + hotspotMaxGermlineVaf} \
        ~{"-hotspot_max_germline_rel_raw_base_qual " + hotspotMaxGermlineRelRawBaseQual} \
        ~{"-panel_max_germline_vaf " + panelMaxGermlineVaf} \
        ~{"-panel_max_germline_rel_raw_base_qual " + panelMaxGermlineRelRawBaseQual} \
        ~{"-mnv_filter_enabled " + mnvFilterEnabled} \
        ~{"-coverage_bed " + coverageBed} \
        ~{true="-panel_only" false="" panelOnly} \
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
        hotspotMinTumorQual: {description: "Equivalent to sage's `hotspot_min_tumor_qual` option.", category: "advanced"}
        panelMinTumorQual: {description: "Equivalent to sage's `panel_min_tumor_qual` option.", category: "advanced"}
        hotspotMaxGermlineVaf: {description: "Equivalent to sage's `hotspot_max_germline_vaf` option.", category: "advanced"}
        hotspotMaxGermlineRelRawBaseQual: {description: "Equivalent to sage's `hotspot_max_germline_rel_raw_base_qual` option.", category: "advanced"}
        panelMaxGermlineVaf: {description: "Equivalent to sage's `panel_max_germline_vaf` option.", category: "advanced"}
        panelMaxGermlineRelRawBaseQual: {description: "Equivalent to sage's `panel_max_germline_vaf` option.", category: "advanced"}
        mnvFilterEnabled: {description: "Equivalent to sage's `mnv_filter_enabled` option.", category: "advanced"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
