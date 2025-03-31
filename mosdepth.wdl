version 1.0

# Copyright (c) 2025 Leiden University Medical Center
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

task Mosdepth {
    input {
        File bam
        File bamIndex
        String prefix = "./out"
        
        String? chrom
        # --by flag takes a BED file or an integer. So there need to be two inputs in WDL's typed system.
        File? byBed
        Int? byWindow
        File? fasta
        Int? flag
        Int? includeFlag

        Boolean noPerBase = false
        Boolean d4 = false
        Boolean fastMode = false
 
        Int threads = 1
        String memory = "1GiB"
        Int timeMinutes = 10 + ceil(size(bam, "G")) * 4
        String dockerImage = "quay.io/biocontainers/mosdepth:0.3.10--h4e814b3_1"
    }

    command <<<
        set -e 
        mkdir -p $(dirname ~{prefix})
        mosdepth \
        --threads ~{threads} \
        ~{"--chrom " + chrom} \
        ~{"--by " + byBed} \
        ~{"--by " + byWindow} \
        ~{"--fasta " + fasta} \
        ~{true="--no-per-base" false="" noPerBase} \
        ~{true="--d4" false="" d4} \
        ~{"--flag " + flag} \
        ~{"--include-flag " + includeFlag} \
        ~{true="--fast-mode" false="" fastMode} \
        ~{prefix} ~{bam} 
    >>>

    output {
        File globalDist = "~{prefix}.mosdepth.global.dist.txt"
        File summary = "~{prefix}.mosdepth.summary.txt"
        File? perBaseBed = "~{prefix}.per-base.bed.gz"
        File? regionsBed = "~{prefix}.regions.bed.gz"
    }

    runtime {
        cpu: threads
        memory: memory 
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bam: {description: "Input BAM or CRAM file.", category: "required"}
        bamIndex: {description: "Index for the input BAM or CRAM file.", category: "required"}
        prefix: {description: "Output prefix.", category: "common"}

        chrom: {description: "Chromosome to restrict depth calculation.", category: "advanced"}
        byBed: {description: "Bed file with windows to include for the --by flag. Should not be used together with byWindow.", category: "common"}
        byWindow: {description: "Integer window size for the --by flag. Should not be used together with byBed.", category: "advanced"}
        fasta: {description: "FASTA file, only necessary when CRAM input is used.", category: "advanced"}
        flag: {description: "Exclude reads with any of the bits in FLAG set.", category: "advanced"}
        includeFlag: {description: "Only include reads with any of the bits in FLAG set.", category: "advanced"}

        noPerBase: {description: "Don't output per-base depth. Skipping this output will speed execution.", category: "common"}
        d4: {description: "output per-base depth in d4 format.", category: "advanced"}
        fastMode: {description: "Don't look at internal cigar operations or correct mate overlaps (recommended for most use-cases).", category: "common"}

        threads: {description: "How many threads to use.", category: "common"}
        memory: {description: "How much memory to allocate.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
        
        # outputs
        globalDist: {description: "Global distribution table file."}
        summary: {description: "Summary table file."}
        perBaseBed: {description: "Per base coverage BED file."}
        regionsBed: {description: "Per region BED file, if byBed or byWindow is used."}
    }
}