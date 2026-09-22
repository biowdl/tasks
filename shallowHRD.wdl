version 1.0

# Copyright (c) 2026 Leiden University Medical Center
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

task ShallowHRD_hg19_controlfreec_chrX {
    input {
        File ratioTxt

        String outputDir = "./"

        String memory = "8GiB"
        Int timeMinutes = 240
        String dockerImage = "quay.io/biowdl/shallowhrd@sha256:0897bc9d840229da4e9a5712e7b6076ef776c03a8286359642ffa277046d2344"
    }

    String sampleName = basename(ratioTxt, ".bam_ratio.txt")

    command {
        set -e
        mkdir -p ~{outputDir}

        Rscript /usr/share/shallowHRD/shallowHRD_hg19_1.13_controlfreec_chrX.R \
        ~{ratioTxt} \
        ~{outputDir} \
        /usr/share/shallowHRD/cytoband_adapted_hg19.csv
        mv Rplots.pdf -t ~{outputDir}
    }

    output {
        File II = "~{outputDir}/~{sampleName}_II.jpeg"
        File III = "~{outputDir}/~{sampleName}_III.jpeg"
        File IV = "~{outputDir}/~{sampleName}_IV.jpeg"
        File IV_txt = "~{outputDir}/~{sampleName}_IV.txt"
        File LGAs = "~{outputDir}/~{sampleName}_LGAs.jpeg"
        File LGAs_txt = "~{outputDir}/~{sampleName}_LGAs.txt"
        File LGAs_intermediary = "~{outputDir}/~{sampleName}_LGAs_intermediary.jpeg"
        File THR = "~{outputDir}/~{sampleName}_THR.jpeg"
        File THR_intermediary = "~{outputDir}/~{sampleName}_THR_intermediary.jpeg"
        File amplification_deletion_table = "~{outputDir}/~{sampleName}_amplification_deletion_table.txt"
        File amplifications_deletions = "~{outputDir}/~{sampleName}_amplifications_deletions.jpeg"
        File beginning_segmentation = "~{outputDir}/~{sampleName}_beginning_segmentation.jpeg"
        File final_segmentation = "~{outputDir}/~{sampleName}_final_segmentation.jpeg"
        File final_segmentation_txt = "~{outputDir}/~{sampleName}_final_segmentation.txt"
        File final_segmentation_intermediary = "~{outputDir}/~{sampleName}_final_segmentation_intermediary.jpeg"
        File final_segmentation_visual = "~{outputDir}/~{sampleName}_final_segmentation_visual.jpeg"
        File final_segmentation_zoomed = "~{outputDir}/~{sampleName}_final_segmentation_zoomed.jpeg"
        File normalised_read_count = "~{outputDir}/~{sampleName}_normalised_read_count.jpeg"
        File number_LGAs = "~{outputDir}/~{sampleName}_number_LGAs.txt"
        File ratio_median_gathered = "~{outputDir}/~{sampleName}_ratio_median_gathered.txt"
        File summary_plot = "~{outputDir}/~{sampleName}_summary_plot.jpeg"
        File rplots = "~{outputDir}/Rplots.pdf"
        Array[File] all = [II, III, IV, IV, LGAs, LGAs,
                           LGAs_intermediary, THR, THR_intermediary,
                           amplification_deletion_table, amplifications_deletions,
                           beginning_segmentation, final_segmentation,
                           final_segmentation, final_segmentation_intermediary,
                           final_segmentation_visual, final_segmentation_zoomed,
                           normalised_read_count, number_LGAs, ratio_median_gathered,
                           summary_plot, rplots]
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        ratioTxt: {description: "The _ratio.txt as produced by control-FREEC.", category: "required"}
        outputDir: {description: "The directory to write the output to.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}