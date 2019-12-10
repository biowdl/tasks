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
        mkdir -p $(dirname ~{outputBed})
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
