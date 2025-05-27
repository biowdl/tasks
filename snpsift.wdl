version 1.0

# MIT License
#
# Copyright (c) 2025 Leiden University Medical Center
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

task SnpSiftFilter {
    input {
        File vcf
        File? vcfIndex
        String filterExpression
        String outputPath = "./snpsift_filter.vcf"

        String memory = "9GiB"
        String javaXmx = "8G"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/snpsift:5.2--hdfd78af_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        SnpSift -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        filter \
        "~{filterExpression}" \
        ~{vcf} \
        > ~{outputPath}
    }

    output {
        File outputVcf = outputPath
    }

    runtime {
        docker: dockerImage
        time_minutes: timeMinutes # !UnknownRuntimeKey
        memory: memory
    }

    parameter_meta {
        vcf: {description: "A VCF file to filter.", category: "required"}
        vcfIndex: {description: "The index for the VCF file.", category: "common"}
        filterExpression: {description: "The SnpSift filtering expression.", category: "required"}
        outputPath: {description: "The path to write the output to.", category: "common"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
