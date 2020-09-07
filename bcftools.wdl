version 1.0

# MIT License
#
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

task Annotate {
    input {
        File? annsFile
        String? collapse
        Array[String] columns = []
        String? exclude
        Boolean force = false
        File? headerLines
        String? newId
        String? include
        Boolean keepSites = false
        String? markSites
        Boolean noVersion = false
        String outputType = "z"
        String? regions
        File? regionsFile
        File? renameChrs
        Array[String] samples = []
        File? samplesFile
        Boolean singleOverlaps = false
        Array[String] removeAnns = []
        File inputFile
        String outputPath = "output.vcf.gz"
        
        Int threads = 0
        String memory = "256M"
        Int timeMinutes = 1 + ceil(size(inputFile, "G"))
        String dockerImage = "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
    }

    Boolean indexing = if outputType == "z" then true else false

    command {
        set -e 
        mkdir -p "$(dirname ~{outputPath})"
        bcftools annotate \
        -o ~{outputPath} \
        -O ~{outputType} \
        ~{"--annotations " + annsFile} \
        ~{"--collapse " + collapse} \
        ~{true="--columns" false="" length(columns) > 0} ~{sep="," columns} \
        ~{"--exclude " + exclude} \
        ~{true="--force" false="" force} \
        ~{"--header-lines " + headerLines} \
        ~{"--set-id " + newId} \
        ~{"--include " + include} \
        ~{true="--keep-sites" false="" keepSites} \
        ~{"--mark-sites " + markSites} \
        ~{true="--no-version" false="" noVersion} \
        ~{"--regions " + regions} \
        ~{"--regions-file " + regionsFile} \
        ~{"--rename-chrs " + renameChrs} \
        ~{true="--samples" false="" length(samples) > 0} ~{sep="," samples} \
        ~{"--samples-file " + samplesFile} \
        ~{true="--single-overlaps" false="" singleOverlaps} \
        ~{true="--remove" false="" length(removeAnns) > 0} ~{sep="," removeAnns} \
        ~{inputFile}

        ~{if indexing then 'bcftools index --tbi ~{outputPath}' else ''}

    }
    
    output {
        File outputVcf = outputPath
        File? outputVcfIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        outputType: {description: "Output type: v=vcf, z=vcf.gz, b=bcf, u=uncompressed bcf", category: "advanced"}
        annsFile: {description: "Bgzip-compressed and tabix-indexed file with annotations (see man page for details).", category: "advanced"}
        collapse: {description: "Treat as identical records with <snps|indels|both|all|some|none>, see man page for details.", category: "advanced"}
        columns: {description: "Comma-separated list of columns or tags to carry over from the annotation file (see man page for details).", category: "advanced"}
        exclude: {description: "Exclude sites for which the expression is true (see man page for details).", category: "advanced"}
        force: {description: "Continue even when parsing errors, such as undefined tags, are encountered.", category: "advanced"}
        headerLines: {description: "Lines to append to the VCF header (see man page for details).", category: "advanced"}
        newId: {description: "Assign ID on the fly (e.g. --set-id +'%CHROM\_%POS').", category: "advanced"}
        include: {description: "Select sites for which the expression is true (see man page for details).", category: "advanced"}
        keepSites: {description: "Keep sites which do not pass -i and -e expressions instead of discarding them.", category: "advanced"}
        markSites: {description: "Annotate sites which are present ('+') or absent ('-') in the -a file with a new INFO/TAG flag.", category: "advanced"}
        noVersion: {description: "Do not append version and command line information to the output VCF header.", category: "advanced"}
        regions: {description: "Restrict to comma-separated list of regions.", category: "advanced"}
        regionsFile: {description: "Restrict to regions listed in a file.", category: "advanced"}
        renameChrs: {description: "rename chromosomes according to the map in file (see man page for details).", category: "advanced"}
        samples: {description: "List of samples for sample stats, \"-\" to include all samples.", category: "advanced"}
        samplesFile: {description: "File of samples to include.", category: "advanced"}
        singleOverlaps: {description: "keep memory requirements low with very large annotation files.", category: "advanced"}
        removeAnns: {description: "List of annotations to remove (see man page for details).", category: "advanced"}
        inputFile: {description: "A vcf or bcf file.", category: "required"}

        threads: {description: "Number of extra decompression threads [0].", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
    }
}

task Sort {
    input {
        File inputFile
        String outputPath = "output.vcf.gz"
        String memory = "256M"
        Int timeMinutes = 1 + ceil(size(inputFile, "G"))
        String dockerImage = "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
        String outputType = "z"
    }

    Boolean indexing = if outputType == "z" then true else false

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        bcftools sort \
        -o ~{outputPath} \
        -O ~{outputType} \
        ~{inputFile}

        ~{if indexing then 'bcftools index --tbi ~{outputPath}' else ''}
    }

    output {
        File outputVcf = outputPath
        File? outputVcfIndex = outputPath + ".tbi"
    }
    
    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        inputFile: {description: "A vcf or bcf file.", category: "required"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        outputType: {description: "Output type: v=vcf, z=vcf.gz, b=bcf, u=uncompressed bcf", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }

    
}

task Stats {
    input {
        File inputVcf
        File inputVcfIndex
        File? compareVcf
        File? compareVcfIndex
        String outputPath = basename(inputVcf) + ".stats"
        String? afBins
        String? afTag
        Boolean firstAlleleOnly = false 
        String? collapse
        String? depth
        String? exclude
        File? exons 
        String? applyFilters
        File? fastaRef
        File? fastaRefIndex
        String? include 
        Boolean splitByID = false 
        String? regions
        File? regionsFile
        Array[String] samples = []
        File? samplesFile 
        String? targets 
        File? targetsFile
        String? userTsTv
        Boolean verbose = false

        Int threads = 0
        Int timeMinutes = 1 + 2* ceil(size(select_all([inputVcf, compareVcf]), "G"))  # TODO: Estimate, 2 minutes per GB, refine later.
        String memory = "256M" 
        String dockerImage = "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
    }
    
    command {
        set -e 
        mkdir -p $(dirname ~{outputPath})
        bcftools stats \
        ~{"--af-bins " + afBins} \
        ~{"--af-tag " + afTag} \
        ~{true="--1st-allele-only" false="" firstAlleleOnly} \
        ~{"--collapse " + collapse} \
        ~{"--depth " + depth} \
        ~{"--exclude " + exclude} \
        ~{"--exons " + exons} \
        ~{"--apply-filters " + applyFilters} \
        ~{"--fasta-ref " + fastaRef} \
        ~{"--include " + include} \
        ~{true="--split-by-ID" false="" splitByID} \
        ~{"--regions " + regions} \
        ~{"--regions-file " + regionsFile} \
        ~{true="--samples" false="" length(samples) > 0} ~{sep="," samples} \
        ~{"--samples-file " + samplesFile} \
        ~{"--targets " + targets} \
        ~{"--targets-file " + targetsFile} \
        ~{"--user-tstv " + userTsTv} \
        --threads ~{threads} \
        ~{true="--verbose" false="" verbose} \
        ~{inputVcf} ~{compareVcf} > ~{outputPath}
    }

    output {
        File stats = outputPath
    }

    runtime {
        cpu: threads + 1
        time_minutes: timeMinutes
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        inputVcf: {description: "The VCF to be analysed.", category: "required"}
        inputVcfIndex: {description: "The index for the input VCF.", category: "required"}
        compareVcf: {description: "When inputVcf and compareVCF are given, the program generates separate stats for intersection and the complements. By default only sites are compared, samples must be given to include also sample columns.", category: "common"}
        compareVcfIndex: {description: "Index for the compareVcf.", category: "common"}
        afBins: {description: "Allele frequency bins, a list (0.1,0.5,1) or a file (0.1\n0.5\n1).", category: "advanced"}
        afTag: {description: "Allele frequency tag to use, by default estimated from AN,AC or GT.", category: "advanded"}
        firstAlleleOnly: {description: "Include only 1st allele at multiallelic sites.", category: "advanced"}
        collapse: {description: "Treat as identical records with <snps|indels|both|all|some|none>, see man page for details.", category: "advanced"}
        depth: {description: "Depth distribution: min,max,bin size [0,500,1].", category: "advanced"}
        exclude: {description: "Exclude sites for which the expression is true (see man page for details).", category: "advanced"}
        exons: {description: "Tab-delimited file with exons for indel frameshifts (chr,from,to; 1-based, inclusive, bgzip compressed).", category: "advanced"}
        applyFilters: {description: "Require at least one of the listed FILTER strings (e.g. \"PASS,.\").", category: "advanced"}
        fastaRef: {description: "Faidx indexed reference sequence file to determine INDEL context.", category: "advanced"}
        fastaRefIndex: {description: "Index file (.fai) for fastaRef. Must be supplied if fastaRef is supplied.", category: "advanced"}
        include: {description: "Select sites for which the expression is true (see man page for details).", category: "advanced"}
        splitByID: {description: "Collect stats for sites with ID separately (known vs novel).", category: "advanced"}
        regions: {description: "Restrict to comma-separated list of regions.", category: "advanced"}
        regionsFile: {description: "Restrict to regions listed in a file.", category: "advanced"}
        samples: {description: "List of samples for sample stats, \"-\" to include all samples.", category: "advanced"}
        samplesFile: {description: "File of samples to include.", category: "advanced"}
        targets: {description: "Similar to regions but streams rather than index-jumps.", category: "advanced"}
        targetsFile: {description: "Similar to regionsFile but streams rather than index-jumps.", category: "advanced"}
        userTsTv: {description: "<TAG[:min:max:n]>. Collect Ts/Tv stats for any tag using the given binning [0:1:100].", category: "advanced"}
        threads: {description: "Number of extra decompression threads [0].", category: "advanced"}
        verbose: {description: "Produce verbose per-site and per-sample output.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
    }
}

task View {
    input {
        File inputFile
        String outputPath = "output.vcf"
        Int compressionLevel = 0
        String memory = "256M"
        Int timeMinutes = 1 + ceil(size(inputFile, "G"))
        String dockerImage = "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
    }

    String outputType = if compressionLevel > 0 then "z" else "v"
    Boolean indexing = if compressionLevel > 0 then true else false
    String outputFilePath = if compressionLevel > 0 then outputPath + ".gz" else outputPath

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        bcftools view \
        -o ~{outputPath} \
        -l ~{compressionLevel} \
        -O ~{outputType} \
        ~{inputFile}

        ~{if indexing then 'bcftools index --tbi ~{outputPath}' else ''}
    }
    output {
        File outputVcf = outputPath
        File? outputVcfIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        inputFile: {description: "A vcf or bcf file.", category: "required"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}
