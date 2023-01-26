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

task PeakCalling {
    input {
        Array[File]+ inputBams
        Array[File]+ inputBamsIndex
        Array[File] controlBams
        Array[File] controlBamsIndex
        String outDir = "macs2"
        String sampleName
        String format = "AUTO"
        Boolean nomodel = false
        String? gensz
        Int? extsize
        Int? shiftsize
        Float? pval_thres
        Boolean? bdg
        String? keepdup
        Boolean? callsummits
        Int timeMinutes = 600  # Default to 10 hours
        String memory = "8GiB"
        String dockerImage = "quay.io/biocontainers/macs2:2.1.2--py27r351_0"
    }

    command {
        set -e
        macs2 callpeak \
        --treatment ~{sep = ' ' inputBams} \
        ~{true="--control" false="" length(controlBams) > 0} ~{sep = ' ' controlBams} \
        --outdir ~{outDir} \
        --name ~{sampleName} \
        ~{"-f" + format} \
        ~{"-g" + gensz} \
        ~{"-p" + pval_thres} \
        ~{"--shift" + shiftsize} \
        ~{"--extsize" + extsize} \
        ~{true='--nomodel' false='' nomodel} \
        ~{true='-B' false='' bdg} \
        ~{"--keep-dup" + keepdup} \
        ~{true='--call-summits' false='' callsummits}
    }

    output {
        File peakFile = outDir + "/" + sampleName + "_peaks.narrowPeak"
    }

    runtime {
        cpu: 1
        memory: memory
        docker: dockerImage
        time_minutes: timeMinutes
    }
    parameter_meta {
        inputBams: {description: "The BAM files on which to perform peak calling.", category: "required"}
        inputBamsIndex: {description: "The indexes for the input BAM files.", category: "required"}
        controlBams: {description: "Control BAM files for the input bam files.", category: "common"}
        controlBamsIndex: {description: "The indexes for the control BAM files.", category: "common"}
        sampleName: {description: "Name of the sample to be analysed", category: "required"}
        outDir: {description: "All output files will be written in this directory.", category: "advanced"}
        nomodel: {description: "Whether or not to build the shifting model.", category: "advanced"}
        gensz: {description: "macs2 argument for setting the mappable genome size or effective genome size which is defined as the genome size which can be sequenced."}
        pval_thres: {description: "macs2 argument for setting the p-value cutoff. If -p is specified, MACS2 will use p-value instead of q-value."}
        shiftsize: {description: "macs2 argument to set an arbitrary shift in bp. Can be negative to indicate direction"}
        extsize: {description: "macs2 argument to extend reads in 5'->3' direction to fix-sized fragments."}
        bdg: {description: "macs2 argument that ebanbles the storage of the fragment pileup, control lambda in bedGraph files."}
        keepdup: {description: "macs2 argument that controls the behavior towards duplicate tags at the exact same location."}
        callsummits: {description: "macs2 argument to reanalyze the shape of signal profile to deconvolve subpeaks within each peak called from the general procedure."}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
        format: {description: "Which format to use. Use BAMPE for paired-end reads.", category: "common"}
    }
}
