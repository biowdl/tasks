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
        File genome
        File bedFile
        String dockerImage = "quay.io/biocontainers/bedtools:2.23.0--hdbcaa40_3"
        String outputFile = basename(bedFile, "\.bed") + ".complement.bed"
    }

    command {
        bedtools complement \
        -g ~{genome} \
        -i ~{bedFile} \
        > ~{outputFile}
    }

    output {
        File complementBed = outputFile
    }

    runtime {
        docker: dockerImage
    }

    parameter_meta {
        genome: {description: "Genome file with names and sizes",
                category: "required"}
        bedFile: {description: "The bedfile to complement",
                category: "required"}
        outputFile: {description: "The path to write the output to",
                     catgory: "advanced"}
        dockerImage: {
            description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
            category: "advanced"
        }
    }
}

# Technically not a bedtools task, but needed for bedtools complement.
task GetChromSizes {
    input {
        File faidx
        # Debian for proper GNU Coreutils. Busybox sucks!
        String dockerImage = "debian@sha256:f05c05a218b7a4a5fe979045b1c8e2a9ec3524e5611ebfdd0ef5b8040f9008fa"
        String outputFile = basename(faidx, "\.fai") + ".genome"
    }

    # Get first two columns from the fasta index which note name and size.
    command {
        cut -f1,2 ~{faidx} \
        > ~{outputFile}
    }

    output {
        File chromSizes = outputFile
    }

    runtime {
        docker: dockerImage
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
        File? g
        File? faidx
        String outputBed = "output.sorted.bed"
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
        ~{"-g " + g} \
        ~{"-faidx" + faidx} \
        > ~{outputBed}
    }

    output {
        File bedFile = outputBed
    }

    runtime {
        docker: dockerImage
    }
}
