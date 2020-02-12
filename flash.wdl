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

import "common.wdl" as common

task Flash {
    input {
        String? preCommand
        FastqPair inputFastq
        String outdirPath
        String outPrefix = "flash"
        Int? minOverlap
        Int? maxOverlap
        Boolean compress = true

        Int threads = 2
        String memory = "2G"
    }

    command {
        set -e -o pipefail
        mkdir -p ~{outdirPath}
        ~{preCommand}
        flash \
        ~{"--threads=" + threads} \
        ~{"--output-directory=" + outdirPath} \
        ~{"--output-prefix=" + outPrefix} \
        ~{true="--compress " false="" compress} \
        ~{"--min-overlap=" + minOverlap} \
        ~{"--max-overlap=" + maxOverlap} \
        ~{inputFastq.R1} ~{inputFastq.R2}
    }

    output {
        File extendedFrags = outdirPath + "/" + outPrefix + ".extendedFrags.fastq.gz"
        File notCombined1 = outdirPath + "/" + outPrefix + ".notCombined_1.fastq.gz"
        File notCombined2 = outdirPath + "/" + outPrefix + ".notCombined_2.fastq.gz"
        FastqPair notCombined = object {
          R1: notCombined1,
          R2: notCombined2
        }
        File hist = outdirPath + "/" + outPrefix + ".hist"
        File histogram = outdirPath + "/" + outPrefix + ".histogram"
    }

    runtime {
        cpu: threads
        memory: memory
    }

}