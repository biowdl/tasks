version 1.0

# MIT License
#
# Copyright (c) 2020 Leiden University Medical Center
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

task Annotation {
    input {
        File vcfFile
        File customs
        String cacheDir
        String cacheVersion
        String outputPath = "./annotated.vcf.gz"
        String? chromosome
        String inFormat = "vcf"
        String outFormat = "vcf"
        Boolean protein = false
        Boolean symbol = false
        Boolean biotype = false
        Boolean regulatory = false
        Boolean phenotype = false
        Boolean stats = false
        Boolean canonical = false
        Boolean coding = false
        Boolean intergenic = false
        Boolean uniprot = false
        String memory = "15G"
        String dockerImage = "quay.io/biocontainers/ensembl-vep:100.1--pl526hecc5488_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        vep \
        --offline \
        --dir_cache ${cacheDir} \
        --cache_version ~{cacheVersion} \
        --input_file ~{vcfFile} \
        --format ~{inFormat} \
        --output_file ~{outputPath} \
        --~{outFormat} \
        --force_overwrite \
        --compress_output bgzip \
        ~{"--chr " + chromosome} \
        ~{true="" false="--no_stats" stats} \
        ~{true="--protein" false="" protein} \
        ~{true="--symbol" false="" symbol} \
        ~{true="--uniprot" false="" uniprot} \
        ~{true="--gene_phenotype" false="" phenotype} \
        ~{true="--biotype" false="" biotype} \
        ~{true="--regulatory" false="" regulatory} \
        ~{true="--canonical" false="" canonical} \
        ~{true="--coding_only" false="" coding} \
        ~{true="" false="--no_intergenic" intergenic} \
        ~{if defined(customs) then prefix("--customs ", read_lines(customs)) else ""}

    }
    
    output {
        File outputVcf = outputPath
    }
 
    runtime {
        memory: memory
        docker: dockerImage
    }

}