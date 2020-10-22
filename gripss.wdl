version 1.0

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

task ApplicationKt {
    input {
        File inputVcf
        String outputPath = "gripss.vcf.gz"
        File referenceFasta
        File breakpointHotspot
        File breakendPon
        File breakpointPon

        String memory = "25G"
        String javaXmx = "24G"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/hmftools-gripss:1.8--0"
    }

    command {
        java -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -cp /usr/local/share/hmftools-gripss-1.8-0/gripss.jar \
        com.hartwig.hmftools.gripss.GripssApplicationKt \
        -ref_genome	~{referenceFasta} \
        -breakpoint_hotspot ~{breakpointHotspot} \
        -breakend_pon ~{breakendPon} \
        -breakpoint_pon ~{breakpointPon} \
        -input_vcf ~{inputVcf} \
        -output_vcf ~{outputPath} 
    }

    output {
        File outputVcf = outputPath
        File outputVcfIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        inputVcf: {description: "The input VCF.", category: "required"}
        outputPath: {description: "The path where th eoutput VCF will be written.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "advanced"}
        breakpointHotspot: {description: "Equivalent to the `-breakpoint_hotspot` option.", category: "required"}
        breakendPon: {description: "Equivalent to the `-breakend_pon` option.", category: "required"}
        breakpointPon: {description: "Equivalent to the `breakpoint_pon` option.", category: "required"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}

task HardFilterApplicationKt {
    input {
        File inputVcf
        String outputPath = "gripss_hard_filter.vcf.gz"

        String memory = "25G"
        String javaXmx = "24G"
        Int timeMinutes = 60
        String dockerImage = "quay.io/biocontainers/hmftools-gripss:1.8--0"
    }

    command {
        java -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -cp /usr/local/share/hmftools-gripss-1.8-0/gripss.jar \
        com.hartwig.hmftools.gripss.GripssHardFilterApplicationKt \
        -input_vcf ~{inputVcf} \
        -output_vcf ~{outputPath} 
    }

    output {
        File outputVcf = outputPath
        File outputVcfIndex = outputPath + ".tbi"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        inputVcf: {description: "The input VCF.", category: "required"}
        outputPath: {description: "The path where th eoutput VCF will be written.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}