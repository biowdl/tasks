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

task Complement {
    input {
        File faidx
        File inputBed
        String outputBed = basename(inputBed, "\.bed") + ".complement.bed"

        String memory = "~{512 + ceil(size([inputBed, faidx], "M"))}M"
        Int timeMinutes = 1 + ceil(size([inputBed, faidx], "G"))
        String dockerImage = "quay.io/biocontainers/bedtools:2.23.0--hdbcaa40_3"
    }

    # Use a fasta index file to get the genome sizes. And convert that to the
    # bedtools specific "genome" format.
    command {
        set -e
        cut -f1,2 ~{faidx} > sizes.genome
        bedtools complement \
        -g sizes.genome \
        -i ~{inputBed} \
        > ~{outputBed}
    }

    output {
        File complementBed = outputBed
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        faidx: {description: "The fasta index (.fai) file from which to extract the genome sizes.", category: "required"}
        inputBed: {description: "The inputBed to complement.", category: "required"}
        outputBed: {description: "The path to write the output to.", category: "advanced"}
        memory: {description: "The amount of memory needed for the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        complementBed: {description: "All intervals in a genome that are not covered by at least one interval in the input file."}
    }
}

task Coverage {
    input {
        File genomeFile
        File a
        File? aIndex
        File b
        File? bIndex
        String outputPath = "./coverage.tsv"

        String memory = "8G"
        Int timeMinutes = 320
        String dockerImage = "quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2"
    }

    command {
        bedtools coverage \
        -sorted \
        -g ~{genomeFile} \
        -a ~{a} \
        -b ~{b} \
        -d \
        > ~{outputPath}
    }

    output {
        File coverageTsv = outputPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        genomeFile: {description: "A file listing the chromosomes and their lengths.", category: "required"}
        a: {description: "The file containing the regions for which the coverage will be counted.", category: "required"}
        aIndex: {description: "An index for the file given as `a`.", category: "common"}
        b: {description: "The file in which the coverage will be counted. Likely a BAM file.", category: "required"}
        bIndex: {description: "An index for the file given as `b`.", category: "common"}
        outputPath: {description: "The path the ouptu will be written to.", category: "common"}

        memory: {description: "The amount of memory needed for the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

    }
}

task Merge {
    input {
        File inputBed
        String outputBed = "merged.bed"

        String memory = "~{512 + ceil(size(inputBed, "M"))}M"
        Int timeMinutes = 1 + ceil(size(inputBed, "G"))
        String dockerImage = "quay.io/biocontainers/bedtools:2.23.0--hdbcaa40_3"
    }

    command {
        set -e
        bedtools merge -i ~{inputBed} > ~{outputBed}
    }

    output {
        File mergedBed = outputBed
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBed: {description: "The bed to merge.", category: "required"}
        outputBed: {description: "The path to write the output to.", category: "advanced"}
        memory: {description: "The amount of memory needed for the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        mergedBed: {description: "Merged bed file."}
    }
}

# Use cat, bedtools sort and bedtools merge to merge bedfiles in a single task.
task MergeBedFiles {
    input {
        Array[File]+ bedFiles
        String outputBed = "merged.bed"

        String memory = "~{512 + ceil(size(bedFiles, "M"))}M"
        Int timeMinutes = 1 + ceil(size(bedFiles, "G"))
        String dockerImage = "quay.io/biocontainers/bedtools:2.23.0--hdbcaa40_3"
    }

    # A sorted bed is needed for bedtools merge
    command {
        set -e -o pipefail
        cat ~{sep=" " bedFiles} | bedtools sort | bedtools merge > ~{outputBed}
    }

    output {
        File mergedBed = outputBed
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bedFiles: {description: "The bed files to merge.", category: "required"}
        outputBed: {description: "The path to write the output to.", category: "advanced"}
        memory: {description: "The amount of memory needed for the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        mergedBed: {description: "Merged bed file."}
    }
}

task Sort {
    input {
        File inputBed
        Boolean sizeA = false
        Boolean sizeD = false
        Boolean chrThenSizeA = false
        Boolean chrThenSizeD = false
        Boolean chrThenScoreA = false
        Boolean chrThenScoreD = false
        String outputBed = "output.sorted.bed"

        File? genome
        File? faidx

        String memory = "~{512 + ceil(size(inputBed, "M"))}M"
        Int timeMinutes = 1 + ceil(size(inputBed, "G"))
        String dockerImage = "quay.io/biocontainers/bedtools:2.23.0--hdbcaa40_3"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputBed})"
        bedtools sort \
        -i ~{inputBed} \
        ~{true="-sizeA" false="" sizeA} \
        ~{true="-sizeD" false="" sizeD} \
        ~{true="-chrThenSizeA" false="" chrThenSizeA} \
        ~{true="-chrThenSizeD" false="" chrThenSizeD} \
        ~{true="-chrThenScoreA" false="" chrThenScoreA} \
        ~{true="-chrThenScoreD" false="" chrThenScoreD} \
        ~{"-g " + genome} \
        ~{"-faidx" + faidx} \
        > ~{outputBed}
    }

    output {
        File sortedBed = outputBed
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBed: {description: "The bed to sort.", category: "required"}
        sizeA: {description: "Sort by feature size in ascending order.", category: "common"}
        sizeD: {description: "Sort by feature size in descending order.", category: "common"}
        chrThenSizeA: {description: "Sort by chromosome (asc), then by feature size (asc).", category: "common"}
        chrThenSizeD: {description: "Sort by chromosome (asc), then by feature size (desc).", category: "common"}
        chrThenScoreA: {description: "Sort by chromosome (asc), then by score (asc).", category: "common"}
        chrThenScoreD: {description: "Sort by chromosome (asc), then by score (desc).", category: "common"}
        outputBed: {description: "The path to write the output to.", category: "advanced"}
        genome: {description: "Define sort order by order of tab-delimited file with chromosome names in the first column.", category: "advanced"}
        faidx: {description: "Define sort order by order of tab-delimited file with chromosome names in the first column. Sort by specified chromosome order.", category: "advanced"}
        memory: {description: "The amount of memory needed for the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        sortedBed: {description: "The sorted bed file."}
    }
}

task Intersect {
    input {
        File regionsA
        File regionsB
        String outputBed = "intersect.bed"

        File? faidx # Giving a faidx file will set the sorted option.

        String memory = "~{512 + ceil(size([regionsA, regionsB], "M"))}M"
        Int timeMinutes = 1 + ceil(size([regionsA, regionsB], "G"))
        String dockerImage = "quay.io/biocontainers/bedtools:2.23.0--hdbcaa40_3"
    }

    Boolean sorted = defined(faidx)

    command {
        set -e
        ~{"cut -f1,2 " + faidx} ~{true="> sorted.genome" false ="" sorted}
        bedtools intersect \
        -a ~{regionsA} \
        -b ~{regionsB} \
        ~{true="-sorted" false="" sorted} \
        ~{true="-g sorted.genome" false="" sorted} \
        > ~{outputBed}
    }

    output {
        File intersectedBed = outputBed
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        regionsA: {description: "Region file a to intersect.", category: "required"}
        regionsB: {description: "Region file b to intersect.", category: "required"}
        outputBed: {description: "The path to write the output to.", category: "advanced"}
        faidx: {description: "The fasta index (.fai) file that is used to create the genome file required for sorted output. Implies sorted option.", category: "common"}
        memory: {description: "The amount of memory needed for the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        intersectedBed: {description: "The intersected bed file."}
    }
}
