task objectMd5 {
    Object the_object

    command {
        cat ${write_object(the_object)} |  md5sum - | sed -e 's/  -//'
    }

    output {
        String md5sum = read_string(stdout())
    }

    runtime {
        memory: 1
    }
}

task mapMd5 {
    Map[String,String] map

    command {
        cat ${write_map(map)} | md5sum - | sed -e 's/  -//'
    }

    output {
        String md5sum = read_string(stdout())
    }

    runtime {
        memory: 1
    }
}

task stringArrayMd5 {
    Array[String] stringArray

    command {
        set -eu -o pipefail
        echo ${sep=',' stringArray} | md5sum - | sed -e 's/  -//'
    }

    output {
        String md5sum = read_string(stdout())
    }

    runtime {
        memory: 1
    }
}

task concatenateTextFiles {
    Array[File] fileList
    String combinedFilePath
    Boolean? unzip=false
    Boolean? zip=false

    command {
        set -e -o pipefail
        ${"mkdir -p $(dirname " + combinedFilePath + ")"}
        ${true='zcat' false= 'cat' unzip} ${sep=' ' fileList} \
        ${true="| gzip -c" false="" zip} > ${combinedFilePath}
    }

    output {
        File combinedFile = combinedFilePath
    }

    runtime {
        memory: 1
    }
}

# inspired by https://gatkforums.broadinstitute.org/wdl/discussion/9616/is-there-a-way-to-flatten-arrays
task flattenStringArray {
    Array[Array[String]] arrayList

    command {
        for line in $(echo ${sep=', ' arrayList}) ; \
        do echo $line | tr -d '"[],' ; done
    }

    output {
        Array[String] flattenedArray = read_lines(stdout())
    }

    runtime {
        memory: 1
    }
}

task appendToStringArray {
    Array[String] array
    String string

    command {
        echo "${sep='\n' array}
        ${string}"
    }

    output {
        Array[String] out_array = read_lines(stdout())
    }

    runtime {
        memory: 1
    }
}

task createLink {
    File inputFile
    String outputPath

    command {
        ln -sf ${inputFile} ${outputPath}
    }

    output {
        File link = outputPath
    }
}