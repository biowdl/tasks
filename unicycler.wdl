task unicycler {
    String? preCommand
    File? short1
    File? short2
    File? unpaired
    File? long
    String out
    Int? verbosity
    Int? minFastaLength
    Int? keep
    Boolean? vcf
    Int? threads
    Int? memory
    Int finalThreads = select_first(threads, 1)
    Int finalMemory = select_first(memory, 4)
    String? mode
    Float? minBridgeQual
    Int? linearSeqs
    File? spadesPath
    Boolean? noCorrect
    Float? minKmerFrac
    Float? maxKmerFrac
    Int? kmerCount
    Float? depthFilter
    Boolean? noMiniasm
    File? raconPath
    File? existingLongReadAssembly
    Boolean? noRotate
    File? startGenes
    Float? startGeneId
    Float? startGeneCov
    String? makeblastdbPath
    File? tblastnPath
    Boolean? noPilon
    File? bowtie2Path
    File? bowtie2buildPath
    File? samtoolsPath
    File? pilonPath
    File? javaPath
    Int? minPolishSize
    File? bcftoolsPath
    Int? minComponentSize
    Int? minDeadEndSize
    File? contamination
    String? scores
    String? lowScore
    command {
        set -e -o pipefail
        mkdir -p ${out}
        ${preCommand}
        unicycler \
        ${"--short1 " + short1} \
        ${"--short2 " + short2} \
        ${"--unpaired " + unpaired} \
        ${"--long " + long} \
        --out ${out} \
        ${"--min_fasta_length " + minFastaLength} \
        ${"--keep " + keep } \
        ${true="--vcf" false="" vcf } \
        ${"--threads " + finalThreads } \
        ${"--mode " + mode } \
        ${"--min_bridge_qual " + minBridgeQual } \
        ${"--linear_seqs " + linearSeqs } \
        ${"--spades_path " + spadesPath } \
        ${true="--no_correct" false="" noCorrect } \
        ${"--min_kmer_frac " + minKmerFrac } \
        ${"--max_kmer_frac " + maxKmerFrac } \
        ${"--kmer_count " + kmerCount } \
        ${"--depth_filter " + depthFilter } \
        ${true="--no_miniasm" false="" noMiniasm } \
        ${"--racon_path " + raconPath } \
        ${"--existing_long_read_assembly " + existingLongReadAssembly } \
        ${true="--no_rotate" false="" noRotate } \
        ${"--start_genes " + startGenes } \
        ${"--start_gene_id " + startGeneId } \
        ${"--start_gene_cov " + startGeneCov } \
        ${"--makeblastdb_path " + makeblastdbPath } \
        ${"--tblastn_path " + tblastnPath } \
        ${true="--no_pilon" false="" noPilon } \
        ${"--bowtie2_path " + bowtie2Path } \
        ${"--bowtie2_build_path " + bowtie2buildPath } \
        ${"--samtools_path " + samtoolsPath } \
        ${"--pilon_path " + pilonPath } \
        ${"--java_path " + javaPath } \
        ${"--min_polish_size " + minPolishSize } \
        ${"--bcftools_path " + bcftoolsPath } \
        ${"--min_component_size " + minComponentSize } \
        ${"--min_dead_end_size " + minDeadEndSize } \
        ${"--contamination " + contamination } \
        ${"--scores " + scores } \
        ${"--low_score " + lowScore }
    }
    output {
        File assemblyFasta = out + "/assembly.fasta"
        File assemblyGfa = out + "/assembly.gfa"
        File log = out + "/unicycler.log"
    }
    runtime {
        cpu: finalThreads
        memory: finalMemory
    }
}