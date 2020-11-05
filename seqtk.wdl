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

task Sample {
    input {
        File sequenceFile
        String outFilePath = "subsampledReads.fq.gz"
        Boolean twoPassMode = false
        Float fractionOrNumber # When above 1.0 is the number of reads, otherwise it's a fraction.
        Boolean zip = true

        String? preCommand
        Int? seed
    }

    command {
        set -e -o pipefail
        mkdir -p "$(dirname ~{outFilePath})"
        ~{preCommand}
        seqtk sample \
        ~{"-s " + seed} \
        ~{true="-2 " false="" twoPassMode} \
        ~{sequenceFile} \
        ~{fractionOrNumber} \
        ~{true="| gzip" false="" zip} \
        >  ~{outFilePath}
    }

    output {
        File subsampledReads = outFilePath
    }
}
