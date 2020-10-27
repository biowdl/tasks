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

task ChunkedScatter {
    input {
        File inputFile
        String prefix = "./scatter"
        Boolean splitContigs = false
        Int? chunkSize
        Int? overlap
        Int? minimumBasesPerFile

        String memory = "256M"
        Int timeMinutes = 2
        String dockerImage = "quay.io/biocontainers/chunked-scatter:1.0.0--py_0"
    }

    command {
        chunked-scatter \
        --print-paths \
        -p ~{prefix} \
        ~{"-c " + chunkSize} \
        ~{"-o " + overlap} \
        ~{"-m " + minimumBasesPerFile} \
        ~{true="--split-contigs " false="" splitContigs} \
        ~{inputFile}
    }

    output {
        Array[File] scatters = read_lines(stdout())
    }

    runtime {
        cpu: 1
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        inputFile: {description: "Either a bed file describing regiosn of intrest or a sequence dictionary.", category: "required"}
        prefix: {description: "The prefix for the output files.", category: "advanced"}
        chunkSize: {description: "Equivalent to chunked-scatter's `-c` option.", category: "advanced"}
        overlap: {description: "Equivalent to chunked-scatter's `-o` option.", category: "advanced"}
        minimumBasesPerFile: {description: "Equivalent to chunked-scatter's `-m` option.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}


task ScatterRegions {
    input {
        File inputFile
        String prefix = "scatters/scatter-" 
        Boolean splitContigs = false
        Int scatterSizeMillions = 1000
        Int? scatterSize
        Int timeMinutes = 2
        String memory = "256M"
        String dockerImage = "quay.io/biocontainers/chunked-scatter:0.2.0--py_0"
    }

    String finalSize = if defined(scatterSize) then "~{scatterSize}" else "~{scatterSizeMillions}000000"
    
    command {
        scatter-regions \
        --print-paths \
        --scatter-size ~{finalSize} \
        ~{true="--split-contigs" false="" splitContigs} \
        ~{"--prefix " + prefix} \
        ~{inputFile} 
    }

    output {
        Array[File] scatters = read_lines(stdout())
    }
    
    runtime {
        cpu: 1
        memory: memory
        docker: dockerImage
        time_minutes: timeMinutes
    }

    parameter_meta {
        inputFile: {description: "The input file, either a bed file or a sequence dict. Which format is used is detected by the extension: '.bed', '.fai' or '.dict'.", category: "required"}
        prefix: {description: "The prefix of the ouput files. Output will be named like: <PREFIX><N>.bed, in which N is an incrementing number. Default 'scatter-'.", category: "advanced"}
        splitContigs: {description: "If set, contigs are allowed to be split up over multiple files.", category: "advanced"}
        scatterSizeMillions: {description: "Over how many million base pairs should be scattered.", category: "common"}
        scatterSize: {description: "Overrides scatterSizeMillions with a smaller value if set.", category: "advanced"}

        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                category: "advanced"}
    }
}
