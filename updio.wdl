version 1.0

# Copyright (c) 2022 Leiden University Medical Center
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

task Updio {
    input {
        File childVcf
        File momVcf 
        File dadVcf
        File? commonCnvFile
        String outputPath = "output_dir"
        Boolean includeX = false
    }

    # "catch .vcf and .gz"
    String outputPrefix = outputPath + "/" + sub(basename(childVcf), "\.vcf(\.gz)?$", "")

    command <<<
        set -e 
        mkdir -p ~{outputPath} 
        updio \
        --output_path ~{outputPath} \
        --child_vcf ~{childVcf} \
        --mom_vcf ~{momVcf} \
        --dad_vcf ~{dadVcf} \
        ~{"--common_cnv_file " + commonCnvFile} \
        ~{true="--include_X=1" false="" includeX}
    >>>

    output {
        File eventsList = outputPrefix + ".events_list"
        File eventsPlot = outputPrefix + ".events_plot.png"
        File log = outputPrefix + ".log"
        File table = outputPrefix + ".table"
        File upd = outputPrefix + ".upd"
        Array[File] files = [eventsList, eventsPlot, log, table, upd]
    }

    runtime {
        # Should be replaced with a tagged version after shake-out
        docker: "quay.io/biowdl/updio:1.0"
    }
}

task UpdioMultisample {
    input {
        File multisampleVcf
        File multisampleVcfIndex
        String childId
        String momId
        String dadId
        File? commonCnvFile
        String outputPath = "output_dir"
        Boolean includeX = false
    }

    # "catch .vcf and .gz"
    String outputPrefix = outputPath + "/" + childId

    command <<<
        set -e
        mkdir -p ~{outputPath}
        updio \
        --output_path ~{outputPath} \
        --multisample_vcf ~{multisampleVcf} \
        --childID ~{childId} \
        --momID ~{momId} \
        --dadID ~{dadId} \
        ~{"--common_cnv_file " + commonCnvFile} \
        ~{true="--include_X=1" false="" includeX}
    >>>

    output {
        File eventsList = outputPrefix + ".events_list"
        File eventsPlot = outputPrefix + ".events_plot.png"
        File log = outputPrefix + ".log"
        File table = outputPrefix + ".table"
        File upd = outputPrefix + ".upd"
        Array[File] files = [eventsList, eventsPlot, log, table, upd]
    }

    runtime {
        # Should be replaced with a tagged version after shake-out
        docker: "quay.io/biowdl/updio:1.0"
    }
}