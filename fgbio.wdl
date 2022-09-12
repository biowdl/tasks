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

task AnnotateBamWithUmis {
    input {
        File inputBam
        File inputUmi
        String outputPath

        String memory = "120GiB"
        Int timeMinutes = 360
        String javaXmx="100G"
        String dockerImage = "quay.io/biocontainers/fgbio:1.4.0--hdfd78af_0"
    }

    command {
        set -e 
        mkdir -p "$(dirname ~{outputPath})"
        fgbio -Xmx~{javaXmx} \
        AnnotateBamWithUmis \
        -i ~{inputBam} \
        -f ~{inputUmi} \
        -o ~{outputPath} 
    }

    output {
        File outputBam = outputPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The input BAM file.", category: "required"}
        inputUmi: {description: "The input fastq file with UMIs.", category: "required"}
        outputPath: {description: "Output directory path + output file.", category: "required"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputBam: {description: "UMI-annotated output BAM file."}
    }
}
