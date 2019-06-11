version 1.0

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
        memory: 1
    }
}

# This task will fail if the MD5sum doesn't match the file.
task CheckFileMD5 {
    input {
        File file
        String md5
        # By default cromwell expects /bin/bash to be present in the container
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

    # When input and output is both compressed decompression is not needed
    String cmdPrefix = if (unzip && !zip) then "zcat " else "cat "
    String cmdSuffix = if (!unzip && zip) then " | gzip -c " else ""

    command {
        set -e -o pipefail
        ~{"mkdir -p $(dirname " + combinedFilePath + ")"}
        ~{cmdPrefix} ~{sep=' ' fileList} ~{cmdSuffix} > ~{combinedFilePath}
    }

    output {
        File combinedFile = combinedFilePath
    }

    runtime {
        memory: 1
    }
}

task Copy {
    input {
        File inputFile
        String outputPath
        Boolean recursive = false

        # Version not that important as long as it is stable.
        String dockerTag = "5.0.2"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPath})
        cp ~{true="-r" false="" recursive} ~{inputFile} ~{outputPath}
    }

    output {
        File outputFile = outputPath
    }

    runtime {
        docker: "bash:" + dockerTag
    }
}

task CreateLink {
    # Making this of type File will create a link to the copy of the file in the execution
    # folder, instead of the actual file.
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
    }

    command {
        cat ~{write_map(map)} | md5sum - | sed -e 's/  -//'
    }

    output {
        String md5sum = read_string(stdout())
    }

    runtime {
        memory: 1
    }
}

task ObjectMd5 {
    input {
        Object the_object
    }

    command {
        cat ~{write_object(the_object)} |  md5sum - | sed -e 's/  -//'
    }

    output {
        String md5sum = read_string(stdout())
    }

    runtime {
        memory: 1
    }
}

task StringArrayMd5 {
    input {
        Array[String] stringArray
    }

    command {
        set -eu -o pipefail
        echo ~{sep=',' stringArray} | md5sum - | sed -e 's/  -//'
    }

    output {
        String md5sum = read_string(stdout())
    }

    runtime {
        memory: 1
    }
}

task YamlToJson {
    input {
        File yaml
        String outputJson = basename(yaml, "\.ya?ml$") + ".json"
        String dockerTag = "3.13-py37-slim"
    }
    command {
        set -e
        mkdir -p $(dirname ~{outputJson})
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
        docker: "biowdl/pyyaml:" + dockerTag
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
