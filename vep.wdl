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

task Vep {
    input {
        File inputFile
        String outputPath = "vep.annotated.vcf.gz"
        File cacheTar
        File? pluginsTar
        String? species 
        Array[String] plugins = []
        Boolean refseq = false 
        Boolean merged = false

        Boolean everything = false
        Boolean symbol = false

        String memory = "8GiB"
        # Account time for unpacking the cache.
        Int timeMinutes = 1 + ceil(size(cacheTar, "GiB")) + ceil(size(inputFile, "MiB") * 3)
        String dockerImage = "quay.io/biocontainers/ensembl-vep:113.3--pl5321h2a3209d_0"
    }

    command <<< 
        set -e
        mkdir vep_cache
        tar -x --directory vep_cache -f ~{cacheTar}
        ~{"tar -x --directory vep_cache -f " + pluginsTar}

        # Output all stats files by default for MultiQC integration
        vep \
        --input_file ~{inputFile} \
        --output_file ~{outputPath} \
        ~{"--species " + species} \
        --stats_html --stats_text \
        --dir vep_cache \
        --offline \
        ~{true="--plugin" false="" length(plugins) > 0} ~{sep=" --plugin " plugins} \
        --vcf \
        --compress_output bgzip \
        ~{true="--refseq" false="" refseq} \
        ~{true="--merged" false="" merged} \
        ~{true="--everything" false="" everything} \
        ~{true="--symbol" false="" symbol} \


        # Cleanup the tar extract to save filesystem space
        rm -rf vep_cache
        

    >>>

    output {
        File outputFile = outputPath
        File statsHtml = outputPath + "_summary.html"
        File statsTxt = outputPath + "_summary.txt"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}
