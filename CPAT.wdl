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

task CPAT {
    input {
        File gene
        String outputPrefix
        File hex
        File logitModel

        File? referenceGenome
        # Should be added as input if CPAT should not index the
        # reference genome.
        File? referenceGenomeIndex
        Array[String]? startCodons
        Array[String]? stopCodons

        String memory = "4G"
        Int timeMinutes = 10 + ceil(size(gene, "G") * 30)
        String dockerImage = "biocontainers/cpat:3.0.4--py39hcbe4a3b_0"
    }

    # Some WDL magic in the command section to properly output the start and
    # stopcodons to the command.
    # select_first is needed in order to convert the optional arrays
    # to non-optionals.
    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        cpat.py \
        --gene ~{gene} \
        --outfile ~{outputPrefix} \
        --hex ~{hex} \
        --logitModel ~{logitModel} \
        ~{"--ref " + referenceGenome} \
        ~{true="--start" false="" defined(startCodons)} ~{sep="," select_first([startCodons, [""]])} \
        ~{true="--stop" false="" defined(stopCodons)} ~{sep="," select_first([stopCodons, [""]])}
    }

    output {
        File orfSeqs = "~{outputPrefix}.ORF_seqs.fa"
        File orfProb = "~{outputPrefix}.ORF_prob.tsv"
        File orfProbBest = "~{outputPrefix}.ORF_prob.best.tsv"
        File noOrf = "~{outputPrefix}.no_ORF.txt"
        File rScript = "~{outputPrefix}.r"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        gene: {description: "Equivalent to CPAT's `--gene` option.", category: "required"}
        outputPrefix: {description: "Equivalent to CPAT's `--outfile` option.", category: "required"}
        hex: {description: "Equivalent to CPAT's `--hex` option.", category: "required"}
        logitModel: {description: "Equivalent to CPAT's `--logitModel` option.", category: "required"}
        referenceGenome: {description: "Equivalent to CPAT's `--ref` option.", category: "advanced"}
        referenceGenomeIndex: {description: "The index of the reference. Should be added as input if CPAT should not index the reference genome.", category: "advanced"}
        startCodons: {description: "Equivalent to CPAT's `--start` option.", category: "advanced"}
        stopCodons: {description: "Equivalent to CPAT's `--stop` option.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}

# There is also make_hexamer_tab.py and make_logitModel.py
# that can be added as tasks here.
