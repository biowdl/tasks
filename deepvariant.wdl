version 1.0

# Copyright (c) 2018 Leiden University Medical Center
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

task RunDeepVariant {
    input {
        File referenceFasta
        File referenceFastaIndex
        File inputBam
        File inputBamIndex
        String modelType
        String outputVcf
        String? postprocessVariantsExtraArgs
        File? customizedModel
        Int? numShards
        String? outputGVcf
        File? regions
        String? sampleName
        Boolean? VCFStatsReport = true

        String memory = "3G"
        Int timeMinutes = 5000
        String dockerImage = "google/deepvariant:1.0.0"
    }

    command {
        set -e

        /opt/deepvariant/bin/run_deepvariant \
        --ref ~{referenceFasta} \
        --reads ~{inputBam} \
        --model_type ~{modelType} \
        --output_vcf ~{outputVcf} \
        ~{"--output_gvcf " + outputGVcf} \
        ~{"--customized_model " + customizedModel} \
        ~{"--num_shards " + numShards} \
        ~{"--regions "  + regions} \
        ~{"--sample_name " + sampleName} \
        ~{"--postprocess_variants_extra_args " + postprocessVariantsExtraArgs} \
        ~{true="--vcf_stats_report" false="--novcf_stats_report" VCFStatsReport}
    }

    runtime {
        docker: dockerImage
        time_minutes: timeMinutes
        memory: memory
    }

    output {
        File outputVCF = outputVcf
        File outputVCFIndex = outputVCF + ".tbi"
        File? outputGVCF = outputGVcf
        File? outputGVCFIndex = outputGVcf + ".tbi"
        Array[File] outputVCFStatsReport = glob("*.visual_report.html")
    }
    
    parameter_meta {
        referenceFasta: {description: "Genome reference to use", category: "required"}
        referenceFastaIndex: {description: "Index for the genome reference file.", category: "required"}
        inputBam: {description: "Aligned, sorted, indexed BAM file containing the reads we want to call.", category: "required"}
        inputBamIndex: {description: "Index for the input bam file.", category: "required"}
        modelType: {description: "<WGS|WES|PACBIO>. Type of model to use for variant calling. Each model_type has an associated default model, which can be overridden by the --customized_model flag", category: "required"}
        outputVcf: {description: "Path where we should write VCF file.", category: "required"}
        customizedModel: {description: "A path to a model checkpoint to load for the `call_variants` step. If not set, the default for each --model_type will be used", category: "advanced"}
        numShards: {description: "Number of shards for make_examples step.", category: "common"}
        outputGVcf: {description: "Path where we should write gVCF file.", category: "common"}
        regions: {description: "List of regions we want to process, in BED/BEDPE format.", category: "advanced"}
        sampleName: {description: "Sample name to use instead of the sample name from the input reads BAM (SM tag in the header).", category: "common"}
        VCFStatsReport: {description: "Output a visual report (HTML) of statistics about the output VCF.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
