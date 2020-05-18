version 1.0

# Copyright (c) 2017 Leiden University Medical Center
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

task GenomeGenerate {
    input {
        String genomeDir = "STAR_index"
        File referenceFasta
        File? referenceGtf
        Int? sjdbOverhang

        Int threads = 4
        String memory = "32G"
        Int timeMinutes = ceil(size(referenceFasta, "G") * 240 / threads)
        String dockerImage = "quay.io/biocontainers/star:2.7.3a--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{genomeDir})"
        STAR \
        --runMode genomeGenerate \
        --runThreadN ~{threads} \
        --genomeDir ~{genomeDir} \
        --genomeFastaFiles ~{referenceFasta} \
        ~{"--sjdbGTFfile " + referenceGtf} \
        ~{"--sjdbOverhang " + sjdbOverhang}
    }

    output {
        File chrLength = "~{genomeDir}/chrLength.txt"
        File chrNameLength = "~{genomeDir}/chrNameLength.txt"
        File chrName = "~{genomeDir}/chrName.txt"
        File chrStart = "~{genomeDir}/chrStart.txt"
        File genome = "~{genomeDir}/genome.txt"
        File genomeParameters = "~{genomeDir}/genomeParameters.txt"
        File sa = "~{genomeDir}/SA"
        File saIndex = "~{genomeDir}/SAindex"
        File? exonGeTrInfo = "~{genomeDir}/exonGeTrInfo.tab"
        File? exonInfo = "~{genomeDir}/exonInfo.tab"
        File? geneInfo = "~{genomeDir}/geneInfo.tab"
        File? sjdbInfo = "~{genomeDir}/sjdbInfo.txt"
        File? sjdbListFromGtfOut = "~{genomeDir}/sjdbList.fromGTF.out.tab"
        File? sjdbListOut = "~{genomeDir}/sjdbList.out.tab"
        File? transcriptInfo = "~{genomeDir}/transcriptInfo.tab"
        Array[File] starIndex = select_all([chrLength, chrNameLength, chrName, chrStart, genome, genomeParameters,
                                            sa, saIndex, exonGeTrInfo, exonInfo, geneInfo, sjdbInfo, sjdbListFromGtfOut,
                                            sjdbListOut, transcriptInfo])
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        genomeDir: {description:"The directory the STAR index should be written to.", categroy: "common"}
        referenceFasta: {description: "The reference Fasta file.", category: "required"}
        referenceGtf: {description: "The reference GTF file.", category: "common"}
        sjdbOverhang: {description: "Equivalent to STAR's `--sjdbOverhang` option.", category: "advanced"}

        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task Star {
    input {
        Array[File]+ inputR1
        Array[File]? inputR2
        Array[File]+ indexFiles
        String outFileNamePrefix
        String outSAMtype = "BAM SortedByCoordinate"
        String readFilesCommand = "zcat"
        String? outStd
        String? twopassMode = "Basic"
        Array[String]? outSAMattrRGline
        String? outSAMunmapped = "Within KeepPairs"
        Int? limitBAMsortRAM

        Int runThreadN = 4
        String memory = "48G"
        String dockerImage = "quay.io/biocontainers/star:2.7.3a--0"
    }

    #TODO Could be extended for all possible output extensions
    Map[String, String] samOutputNames = {"BAM SortedByCoordinate": "sortedByCoord.out.bam"}

    command {
        set -e
        mkdir -p "$(dirname ~{outFileNamePrefix})"
        STAR \
        --readFilesIn ~{sep=',' inputR1} ~{sep="," inputR2} \
        --outFileNamePrefix ~{outFileNamePrefix} \
        --genomeDir ~{sub(indexFiles[0], basename(indexFiles[0]), "")} \
        --outSAMtype ~{outSAMtype} \
        --readFilesCommand ~{readFilesCommand} \
        ~{"--outSAMunmapped " + outSAMunmapped} \
        ~{"--runThreadN " + runThreadN} \
        ~{"--outStd " + outStd} \
        ~{"--twopassMode " + twopassMode} \
        ~{"--limitBAMsortRAM " + limitBAMsortRAM} \
        ~{true="--outSAMattrRGline " false="" defined(outSAMattrRGline)} ~{sep=" , " outSAMattrRGline}
    }

    output {
        File bamFile = outFileNamePrefix + "Aligned." +  samOutputNames[outSAMtype]
        File logFinalOut = outFileNamePrefix + "Log.final.out"
    }

    runtime {
        cpu: runThreadN
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        inputR1: {description: "The first-/single-end FastQ files.", category: "required"}
        inputR2: {description: "The second-end FastQ files (in the same order as the first-end files).", category: "common"}
        indexFiles: {description: "The star index files.", category: "required"}
        outFileNamePrefix: {description: "The prefix for the output files. May include directories.", category: "required"}
        outSAMtype: {description: "The type of alignment file to be produced. Currently only `BAM SortedByCoordinate` is supported.", category: "advanced"}
        readFilesCommand: {description: "Equivalent to star's `--readFilesCommand` option.", category: "advanced"}
        outStd: {description: "Equivalent to star's `--outStd` option.", category: "advanced"}
        twopassMode: {description: "Equivalent to star's `--twopassMode` option.", category: "advanced"}
        outSAMattrRGline: {description: "The readgroup lines for the fastq pairs given (in the same order as the fastq files).", category: "common"}
        outSAMunmapped: {description: "Equivalent to star's `--outSAMunmapped` option.", category: "advanced"}
        limitBAMsortRAM: {description: "Equivalent to star's `--limitBAMsortRAM` option.", category: "advanced"}
        runThreadN: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task MakeStarRGline {
    input {
        String sample
        String library
        String platform = "ILLUMINA"
        String readgroup
    }

    command {
        printf '"ID:~{readgroup}" "LB:~{library}" "PL:~{platform}" "SM:~{sample}"'
    }

    output {
        String rgLine = read_string(stdout())
    }
}
