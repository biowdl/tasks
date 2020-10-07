version 1.0

# Copyright (c) 2018 Leiden University Medical Center
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


task Phase {
    input {
        String outputVCF
        File? reference
        File? referenceIndex
        Boolean? no_reference
        String? tag
        File? output_read_list
        String? algorithm
        Boolean? merge_reads
        String? internal_downsampling
        String? mapping_quality
        Boolean? indels
        Boolean? ignore_read_groups
        String? sample
        String? chromosome
        String? error_rate
        String? maximum_error_rate
        String? threshold
        String? negative_threshold
        Boolean? full_genotyping
        Boolean? distrust_genotypes
        Boolean? include_homozygous
        String? default_gq
        String? gl_regularize_r
        File? changed_genotype_list
        String? ped
        File? recombination_list
        String? recomb_rate
        File? gen_map
        Boolean? no_genetic_haplo_typing
        Boolean? use_ped_samples
        File vcf
        File vcfIndex
        File phaseInput
        File phaseInputIndex

        String memory = "4G"
        Int timeMinutes = 120
        # Whatshap 1.0, tabix 0.2.5
        String dockerImage = "quay.io/biocontainers/mulled-v2-5c61fe1d8c284dd05d26238ce877aa323205bf82:89b4005d04552bdd268e8af323df83357e968d83-0"
    }

    command {
        whatshap phase \
        ~{vcf} \
        ~{phaseInput} \
        ~{if defined(outputVCF) then ("--output " +  '"' + outputVCF + '"') else ""} \
        ~{if defined(reference) then ("--reference " +  '"' + reference + '"') else ""} \
        ~{true="--no-reference" false="" no_reference} \
        ~{if defined(tag) then ("--tag " +  '"' + tag + '"') else ""} \
        ~{if defined(output_read_list) then ("--output-read-list " +  '"' + output_read_list + '"') else ""} \
        ~{if defined(algorithm) then ("--algorithm " +  '"' + algorithm + '"') else ""} \
        ~{true="--merge-reads" false="" merge_reads} \
        ~{if defined(internal_downsampling) then ("--internal-downsampling " +  '"' + internal_downsampling + '"') else ""} \
        ~{if defined(mapping_quality) then ("--mapping-quality " +  '"' + mapping_quality + '"') else ""} \
        ~{true="--indels" false="" indels} \
        ~{true="--ignore-read-groups" false="" ignore_read_groups} \
        ~{if defined(sample) then ("--sample " +  '"' + sample + '"') else ""} \
        ~{if defined(chromosome) then ("--chromosome " +  '"' + chromosome + '"') else ""} \
        ~{if defined(error_rate) then ("--error-rate " +  '"' + error_rate + '"') else ""} \
        ~{if defined(maximum_error_rate) then ("--maximum-error-rate " +  '"' + maximum_error_rate + '"') else ""} \
        ~{if defined(threshold) then ("--threshold " +  '"' + threshold + '"') else ""} \
        ~{if defined(negative_threshold) then ("--negative-threshold " +  '"' + negative_threshold + '"') else ""} \
        ~{true="--full-genotyping" false="" full_genotyping} \
        ~{true="--distrust-genotypes" false="" distrust_genotypes} \
        ~{true="--include-homozygous" false="" include_homozygous} \
        ~{if defined(default_gq) then ("--default-gq " +  '"' + default_gq + '"') else ""} \
        ~{if defined(gl_regularize_r) then ("--gl-regularizer " +  '"' + gl_regularize_r + '"') else ""} \
        ~{if defined(changed_genotype_list) then ("--changed-genotype-list " +  '"' + changed_genotype_list + '"') else ""} \
        ~{if defined(ped) then ("--ped " +  '"' + ped + '"') else ""} \
        ~{if defined(recombination_list) then ("--recombination-list " +  '"' + recombination_list + '"') else ""} \
        ~{if defined(recomb_rate) then ("--recombrate " +  '"' + recomb_rate + '"') else ""} \
        ~{if defined(gen_map) then ("--genmap " +  '"' + gen_map + '"') else ""} \
        ~{true="--no-genetic-haplotyping" false="" no_genetic_haplo_typing} \
        ~{true="--use-ped-samples" false="" use_ped_samples} && \
        tabix -p vcf ~{outputVCF}
    }

    output {
        File phasedVCF = outputVCF
        File phasedVCFIndex = outputVCF + ".tbi"
    }

    runtime {
        docker: dockerImage
        time_minutes: timeMinutes
        memory: memory
    }

    parameter_meta {
        outputVCF: {description: "Output VCF file. Add .gz to the file name to get compressed output. If omitted, use standard output.", category: "common"}
        reference: {description: "Reference file. Provide this to detect alleles through re-alignment. If no index (.fai) exists, it will be created", category: "common"}
        no_reference: {description: "Detect alleles without requiring a reference, at the expense of phasing quality (in particular for long reads)", category: "common"}
        tag: {description: "Store phasing information with PS tag (standardized) or HP tag (used by GATK ReadBackedPhasing) (default: {description: PS)", category: "common"}
        output_read_list: {description: "Write reads that have been used for phasing to FILE.", category: "advanced"}
        algorithm: {description: "Phasing algorithm to use (default: {description: whatshap)", category: "advanced"}
        merge_reads: {description: "Merge reads which are likely to come from the same haplotype (default: {description: do not merge reads)", category: "common"}
        internal_downsampling: {description: "Coverage reduction parameter in the internal core phasing algorithm. Higher values increase runtime *exponentially* while possibly improving phasing quality marginally. Avoid using this in the normal case! (default: {description: 15)", category: "advanced"}
        mapping_quality: {description: "Minimum mapping quality (default: {description: 20)", category: "common"}
        indels: {description: "Also phase indels (default: {description: do not phase indels)", category: "common"}
        ignore_read_groups: {description: "Ignore read groups in BAM/CRAM header and assume all reads come from the same sample.", category: "advanced"}
        sample: {description: "Name of a sample to phase. If not given, all samples in the input VCF are phased. Can be used multiple times.", category: "common"}
        chromosome: {description: "Name of chromosome to phase. If not given, all chromosomes in the input VCF are phased. Can be used multiple times.", category: "common"}
        error_rate: {description: "The probability that a nucleotide is wrong in read merging model (default: {description: 0.15).", category: "advanced"}
        maximum_error_rate: {description: "The maximum error rate of any edge of the read merging graph before discarding it (default: {description: 0.25).", category: "advanced"}
        threshold: {description: "The threshold of the ratio between the probabilities that a pair of reads come from the same haplotype and different haplotypes in the read merging model (default: {description: 1000000).", category: "advanced"}
        negative_threshold: {description: "The threshold of the ratio between the probabilities that a pair of reads come from different haplotypes and the same haplotype in the read merging model (default: {description: 1000).", category: "advanced"}
        full_genotyping: {description: "Completely re-genotype all variants based on read data, ignores all genotype data that might be present in the VCF (EXPERIMENTAL FEATURE).", category: "experimental"}
        distrust_genotypes: {description: "Allow switching variants from hetero- to homozygous in an optimal solution (see documentation).", category: "advanced"}
        include_homozygous: {description: "Also work on homozygous variants, which might be turned to heterozygous", category: "advanced"}
        default_gq: {description: "Default genotype quality used as cost of changing a genotype when no genotype likelihoods are available (default 30)", category: "advanced"}
        gl_regularize_r: {description: "Constant (float) to be used to regularize genotype likelihoods read from input VCF (default None).", category: "advanced"}
        changed_genotype_list: {description: "Write list of changed genotypes to FILE.", category: "advanced"}
        ped: {description: "Use pedigree information in PED file to improve phasing (switches to PedMEC algorithm). Columns 2, 3, 4 must refer to child, mother, and father sample names as used in the VCF and BAM/CRAM. Other columns are ignored.", category: "advanced"}
        recombination_list: {description: "Write putative recombination events to FILE.", category: "advanced"}
        recomb_rate: {description: "Recombination rate in cM/Mb (used with --ped). If given, a constant recombination rate is assumed (default: {description: 1.26cM/Mb).", category: "advanced"}
        gen_map: {description: "File with genetic map (used with --ped) to be used instead of constant recombination rate, i.e. overrides option --recombrate.", category: "advanced"}
        no_genetic_haplo_typing: {description: "Do not merge blocks that are not connected by reads (i.e. solely based on genotype status). Default: {description: when in --ped mode, merge all blocks that contain at least one homozygous genotype in at least one individual into one block.", category: "advanced"}
        use_ped_samples: {description: "Only work on samples mentioned in the provided PED file.", category: "advanced"}
        vcf: {description: "VCF or BCF file with variants to be phased (can be gzip-compressed)", category: "required"}
        vcfIndex: {description: "Index for the VCF or BCF file with variants to be phased", category: "required"}
        phaseInput: {description: "BAM, CRAM, VCF or BCF file(s) with phase information, either through sequencing reads (BAM, CRAM) or through phased blocks (VCF, BCF)", category: "required"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}

task Stats {
    input {
        String? gtf
        String? sample
        String? chr_lengths
        String? tsv
        Boolean? only_sn_vs
        String? block_list
        String? chromosome
        File vcf

        String memory = "4G"
        Int timeMinutes = 120
        # Whatshap 1.0, tabix 0.2.5
        String dockerImage = "quay.io/biocontainers/mulled-v2-5c61fe1d8c284dd05d26238ce877aa323205bf82:89b4005d04552bdd268e8af323df83357e968d83-0"
      }

    command {
      whatshap stats \
        ~{vcf} \
        ~{if defined(gtf) then ("--gtf " +  '"' + gtf + '"') else ""} \
        ~{if defined(sample) then ("--sample " +  '"' + sample + '"') else ""} \
        ~{if defined(chr_lengths) then ("--chr-lengths " +  '"' + chr_lengths + '"') else ""} \
        ~{if defined(tsv) then ("--tsv " +  '"' + tsv + '"') else ""} \
        ~{true="--only-snvs" false="" only_sn_vs} \
        ~{if defined(block_list) then ("--block-list " +  '"' + block_list + '"') else ""} \
        ~{if defined(chromosome) then ("--chromosome " +  '"' + chromosome + '"') else ""}
    }

    output {
      File? phasedGTF = gtf
      File? phasedTSV = tsv
      File? phasedBlockList = block_list
    }

    runtime {
        docker: dockerImage
        time_minutes: timeMinutes
        memory: memory
    }

    parameter_meta {
        gtf: "Write phased blocks to GTF file."
        sample: "Name of the sample to process. If not given, use first sample found in VCF."
        chr_lengths: "File with chromosome lengths (one line per chromosome, tab separated '<chr> <length>') needed to compute N50 values."
        tsv: "Filename to write statistics to (tab-separated)."
        only_sn_vs: "Only process SNVs and ignore all other variants."
        block_list: "Filename to write list of all blocks to (one block per line)."
        chromosome: "Name of chromosome to process. If not given, all chromosomes in the input VCF are considered."
        vcf: "Phased VCF file"
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}

task Haplotag {
    input {
        String outputFile
        File? reference
        File? referenceFastaIndex
        String? regions
        Boolean? ignore_linked_read
        String? linked_read_distance_cut_off
        Boolean? ignore_read_groups
        String? sample
        String? output_haplo_tag_list
        Boolean? tag_supplementary
        File vcf
        File vcfIndex
        File alignments
        File alignmentsIndex

        String memory = "4G"
        Int timeMinutes = 120
        # Whatshap 1.0, tabix 0.2.5
        String dockerImage = "quay.io/biocontainers/mulled-v2-5c61fe1d8c284dd05d26238ce877aa323205bf82:89b4005d04552bdd268e8af323df83357e968d83-0"
    }

    command {
      whatshap haplotag \
          ~{vcf} \
          ~{alignments} \
          ~{if defined(outputFile) then ("--output " +  '"' + outputFile+ '"') else ""} \
          ~{if defined(reference) then ("--reference " +  '"' + reference + '"') else ""} \
          ~{if defined(regions) then ("--regions " +  '"' + regions + '"') else ""} \
          ~{true="--ignore-linked-read" false="" ignore_linked_read} \
          ~{if defined(linked_read_distance_cut_off) then ("--linked-read-distance-cutoff " +  '"' + linked_read_distance_cut_off + '"') else ""} \
          ~{true="--ignore-read-groups" false="" ignore_read_groups} \
          ~{if defined(sample) then ("--sample " +  '"' + sample + '"') else ""} \
          ~{if defined(output_haplo_tag_list) then ("--output-haplotag-list " +  '"' + output_haplo_tag_list + '"') else ""} \
          ~{true="--tag-supplementary" false="" tag_supplementary} && \
          python3 -c "import pysam; pysam.index('~{outputFile}')"
    }

    output {
      File bam = outputFile
      File bamIndex = outputFile + ".bai"
    }

    runtime {
        docker: dockerImage
        time_minutes: timeMinutes
        memory: memory
    }

    parameter_meta {
        outputFile: "Output file. If omitted, use standard output."
        reference: "Reference file. Provide this to detect alleles through re-alignment. If no index (.fai) exists, it will be created"
        regions: "Specify region(s) of interest to limit the tagging to reads/variants overlapping those regions. You can specify a space-separated list of regions in the form of chrom:start-end, chrom (consider entire chromosome), or chrom:start (consider region from this start to end of chromosome)."
        ignore_linked_read: "Ignore linkage information stored in BX tags of the reads."
        linked_read_distance_cut_off: "Assume reads with identical BX tags belong to different read clouds if their distance is larger than LINKEDREADDISTANCE (default: 50000)."
        ignore_read_groups: "Ignore read groups in BAM/CRAM header and assume all reads come from the same sample."
        sample: "Name of a sample to phase. If not given, all samples in the input VCF are phased. Can be used multiple times."
        output_haplo_tag_list: "Write assignments of read names to haplotypes (tab separated) to given output file. If filename ends in .gz, then output is gzipped."
        tag_supplementary: "Also tag supplementary alignments. Supplementary alignments are assigned to the same haplotype the primary alignment has been assigned to (default: only tag primary alignments)."
        vcf: "VCF file with phased variants (must be gzip-compressed and indexed)"
        alignments: "File (BAM/CRAM) with read alignments to be tagged by haplotype"
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}
