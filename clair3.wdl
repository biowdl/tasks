version 1.0

# Copyright (c) 2024 Leiden University Medical Center
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

task Clair3 {
    input {
        File bam 
        File bamIndex 
        File referenceFasta 
        File referenceFastaFai
        String outputPrefix 
        String? sampleName
        File? modelTar
        String? builtinModel
        String platform
        Int threads = 8
        Boolean includeAllCtgs = false
        String memory = "~{threads + 16}GiB"
        Int timeMinutes = 10 + ceil(size(bam, "G") * 200 / threads)
        String dockerImage = "quay.io/biocontainers/clair3:1.0.10--py39h46983ab_0"   
    }

    String modelArg = "~{if defined(modelTar) then basename(select_first([modelTar]), '.tar.gz') else builtinModel}"

    command <<<
    set -e
    ~{if defined(modelTar) then "tar -xvf " + modelTar else "" }
    mkdir -p $(dirname ~{outputPrefix})
    run_clair3.sh \
    --model=~{modelArg} \
    --ref_fn=~{referenceFasta} \
    --bam_fn=~{bam} \
    --output=out \
    --threads=~{threads} \
    --platform=~{platform} \
    ~{"--sample_name=" + sampleName} \
    ~{true="--include_all_ctgs" false ="" includeAllCtgs}  
    mv out/merge_output.vcf.gz ~{outputPrefix}.vcf.gz
    mv out/merge_output.vcf.gz.tbi ~{outputPrefix}.vcf.gz.tbi
    >>>

    output {
        File vcf = "~{outputPrefix}.vcf.gz"
        File vcfIndex = "~{outputPrefix}.vcf.gz.tbi"  
    }

    runtime {
        cpu: threads 
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    } 

}