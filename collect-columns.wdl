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

task CollectColumns {
    input {
        Array[File]+ inputTables
        String outputPath
        Int? featureColumn
        Int? valueColumn
        Int? separator
        Array[String]? sampleNames
        Boolean header = false
        Boolean sumOnDuplicateId = false
        Array[String]? additionalAttributes
        File? referenceGtf
        String? featureAttribute

        Int memoryGb = 4 + ceil(0.5 * length(inputTables))
        Int timeMinutes = 10
        String dockerImage = "quay.io/biocontainers/collect-columns:1.0.0--py_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        collect-columns \
        ~{outputPath} \
        ~{sep=" " inputTables} \
        ~{"-f "  + featureColumn} \
        ~{"-c " + valueColumn} \
        ~{"-s " + separator} \
        ~{true="-n" false="" defined(sampleNames)} ~{sep=" " sampleNames} \
        ~{true="-H" false="" header} \
        ~{true="-S" false="" sumOnDuplicateId} \
        ~{true="-a" false="" defined(additionalAttributes)} ~{sep=" " additionalAttributes} \
        ~{"-g " + referenceGtf} \
        ~{"-F " + featureAttribute}
    }

    output {
        File outputTable = outputPath
    }

    runtime {
        memory: "~{memoryGb}G"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        inputTables: {description: "The tables from which columns should be taken.", category: "required"}
        outputPath: {description: "The path to which the output should be written.", category: "required"}
        featureColumn: {description: "Equivalent to the -f option of collect-columns.", category: "advanced"}
        valueColumn: {description: "Equivalent to the -c option of collect-columns.", category: "advanced"}
        separator: {description: "Equivalent to the -s option of collect-columns.", category: "advanced"}
        sampleNames: {description: "Equivalent to the -n option of collect-columns.", category: "advanced"}
        header: {description: "Equivalent to the -H flag of collect-columns.", category: "advanced"}
        sumOnDuplicateId: {description: "Equivalent to the -S flag of collect-columns.", category: "advanced"}
        additionalAttributes: {description: "Equivalent to the -a option of collect-columns.", category: "advanced"}
        referenceGtf: {description: "Equivalent to the -g option of collect-columns.", category: "advanced"}
        featureAttribute: {description: "Equivalent to the -F option of collect-columns.", category: "advanced"}
        memoryGb: {description: "The maximum amount of memory the job will need in GB", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}