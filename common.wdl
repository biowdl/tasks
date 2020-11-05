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

task AppendToStringArray {
    input {
        Array[String] array
        String string
    }

    command {
        echo "~{sep='\n' array}
        ~{string}"
    }

    output {
        Array[String] outArray = read_lines(stdout())
    }

    runtime {
        memory: "1G"
    }
}

# This task will fail if the MD5sum doesn't match the file.
task CheckFileMD5 {
    input {
        File file
        String md5
        # By default cromwell expects /bin/bash to be present in the container.
        # The 'bash' container does not fill this requirement. (It is in /usr/local/bin/bash)
        # Use a stable version of debian:stretch-slim for this. (Smaller than ubuntu)
        String dockerImage = "debian@sha256:f05c05a218b7a4a5fe979045b1c8e2a9ec3524e5611ebfdd0ef5b8040f9008fa"
    }

    command {
        bash -c '
        set -e -o pipefail
        echo "~{md5}  ~{file}" | md5sum -c
        '
    }

    runtime {
        docker: dockerImage
    }
}

task ConcatenateTextFiles {
    input {
        Array[File] fileList
        String combinedFilePath
        Boolean unzip = false
        Boolean zip = false
    }

    # When input and output is both compressed decompression is not needed.
    String cmdPrefix = if (unzip && !zip) then "zcat " else "cat "
    String cmdSuffix = if (!unzip && zip) then " | gzip -c " else ""

    command {
        set -e -o pipefail
        mkdir -p "$(dirname ~{combinedFilePath})"
        ~{cmdPrefix} ~{sep=" " fileList} ~{cmdSuffix} > ~{combinedFilePath}
    }

    output {
        File combinedFile = combinedFilePath
    }

    runtime {
        memory: "1G"
    }
}

task Copy {
    input {
        File inputFile
        String outputPath
        Boolean recursive = false

        # Version not that important as long as it is stable.
        String dockerImage = "debian@sha256:f05c05a218b7a4a5fe979045b1c8e2a9ec3524e5611ebfdd0ef5b8040f9008fa"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        cp ~{true="-r" false="" recursive} ~{inputFile} ~{outputPath}
    }

    output {
        File outputFile = outputPath
    }

    runtime {
        docker: dockerImage
    }
}

task CreateLink {
    # Making this of type File will create a link to the copy of the file in
    # the execution folder, instead of the actual file.
    # This cannot be propperly call-cached or used within a container.
    input {
        String inputFile
        String outputPath
    }

    command {
        ln -sf ~{inputFile} ~{outputPath}
    }

    output {
        File link = outputPath
    }
}

task MapMd5 {
    input {
        Map[String,String] map

        String dockerImage = "debian@sha256:f05c05a218b7a4a5fe979045b1c8e2a9ec3524e5611ebfdd0ef5b8040f9008fa"
    }

    command {
        set -e -o pipefail
        md5sum "~{write_map(map)}" | cut -f 1 -d ' '
    }

    output {
        String md5sum = read_string(stdout())
    }

    runtime {
        memory: "1G"
        docker: dockerImage
    }
}


task StringArrayMd5 {
    input {
        Array[String] stringArray

        String dockerImage = "debian@sha256:f05c05a218b7a4a5fe979045b1c8e2a9ec3524e5611ebfdd0ef5b8040f9008fa"
    }

    command {
        set -eu -o pipefail
        echo ~{sep=',' stringArray} | md5sum - | sed -e 's/  -//'
    }

    output {
        String md5sum = read_string(stdout())
    }

    runtime {
        memory: "1G"
        docker: dockerImage
    }
}

task TextToFile {
    input {
        String text
        String outputFile = "out.txt"

        Int timeMinutes = 1
        String dockerImage = "debian@sha256:f05c05a218b7a4a5fe979045b1c8e2a9ec3524e5611ebfdd0ef5b8040f9008fa"
    }

    command <<<
        echo ~{text} > ~{outputFile}
    >>>

    output {
        File out = outputFile
    }

    runtime {
        memory: "1G"
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        text: {description: "The text to print.", category: "required"}
        outputFile: {description: "The name of the output file.", category: "common"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}

task YamlToJson {
    input {
        File yaml
        String outputJson = basename(yaml, "\.ya?ml$") + ".json"

        String  memory = "128M"
        Int timeMinutes = 1
        # biowdl-input-converter has python and pyyaml.
        String dockerImage = "quay.io/biocontainers/biowdl-input-converter:0.2.1--py_0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputJson})"
        python <<CODE
        import json
        import yaml
        with open("~{yaml}", "r") as input_yaml:
            content = yaml.load(input_yaml)
        with open("~{outputJson}", "w") as output_json:
            json.dump(content, output_json)
        CODE
    }

    output {
        File json = outputJson
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        yaml: {description: "The YAML file to convert.", category: "required"}
        outputJson: {description: "The location the output JSON file should be written to.", category: "advanced"}
        memory: {description: "The maximum amount of memory the job will need.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}

struct Reference {
    File fasta
    File fai
    File dict
}

struct IndexedVcfFile {
    File file
    File index
    String? md5sum
}

struct IndexedBamFile {
    File file
    File index
    String? md5sum
}

struct FastqPair {
    File R1
    String? R1_md5
    File? R2
    String? R2_md5
}

struct CaseControl {
    String inputName
    IndexedBamFile inputFile
    String controlName
    IndexedBamFile controlFile
}

struct CaseControls {
    Array[CaseControl] caseControls
}
