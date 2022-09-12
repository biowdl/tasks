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

task Spades {
    input {
        File read1
        File? read2
        String outputDir

        String? preCommand
        File? interlacedReads
        File? sangerReads
        File? pacbioReads
        File? nanoporeReads
        File? tslrContigs
        File? trustedContigs
        File? untrustedContigs
        Boolean singleCell = false
        Boolean metagenomic = false
        Boolean rna = false
        Boolean plasmid = false
        Boolean ionTorrent = false
        Boolean onlyErrorCorrection = false
        Boolean onlyAssembler = false
        Boolean careful = false
        Boolean disableGzipOutput = false
        Boolean disableRepeatResolution = false
        File? dataset
        File? tmpDir
        String? k
        Float? covCutoff
        Int? phredOffset

        Int threads = 1
        Int memoryGb = 16
    }

    command {
        set -e -o pipefail
        ~{preCommand}
        spades.py \
        ~{"-o " + outputDir} \
        ~{true="--sc" false="" singleCell} \
        ~{true="--meta" false="" metagenomic} \
        ~{true="--rna" false="" rna} \
        ~{true="--plasmid" false="" plasmid} \
        ~{true="--iontorrent" false="" ionTorrent} \
        ~{"--12 " + interlacedReads} \
        ~{true="-1" false="-s" defined(read2)} ~{read1}  \
        ~{"-2 " + read2} \
        ~{"--sanger " + sangerReads} \
        ~{"--pacbio " + pacbioReads} \
        ~{"--nanopore " + nanoporeReads} \
        ~{"--tslr " + tslrContigs} \
        ~{"--trusted-contigs " + trustedContigs} \
        ~{"--untrusted-contigs " + untrustedContigs} \
        ~{true="--only-error-correction" false="" onlyErrorCorrection} \
        ~{true="--only-assembler" false="" onlyAssembler} \
        ~{true="--careful" false="" careful} \
        ~{true="--disable-gzip-output" false="" disableGzipOutput} \
        ~{true="--disable-rr" false="" disableRepeatResolution} \
        ~{"--dataset " + dataset} \
        ~{"--threads " + threads} \
        ~{"--memory " + memoryGb} \
        ~{"-k " + k} \
        ~{"--cov-cutoff " + covCutoff} \
        ~{"--phred-offset " + phredOffset}
    }

    output {
        Array[File] correctedReads = glob(outputDir + "/corrected/*.fastq*")
        File scaffolds = outputDir + "/scaffolds.fasta"
        File contigs = outputDir + "/contigs.fasta"
        File assemblyGraphWithScaffoldsGfa = outputDir + "/assembly_graph_with_scaffolds.gfa"
        File assemblyGraphFastg = outputDir + "/assembly_graph.fastg"
        File contigsPaths = outputDir + "/contigs.paths"
        File scaffoldsPaths = outputDir + "/scaffolds.paths"
        File params = outputDir + "/params.txt"
        File log = outputDir + "/spades.log"
    }

    runtime {
        cpu: threads
        memory: "~{memoryGb}GiB"
    }
}
