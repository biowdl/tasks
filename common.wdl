task objectMd5 {
    Object the_object
    command {
        cat ${write_object(the_object)} |  md5sum - | sed -e 's/  -//'
    }
    output {
        String md5sum = read_string(stdout())
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
}

task concatenateTextFiles {
    Array[File] fileList
    String combinedFilePath
    Boolean? unzip=false
    command {
        set -e -o pipefail
        ${"mkdir -p $(dirname " + combinedFilePath + ")"}
        ${true='zcat' false= 'cat' unzip} ${sep=' ' fileList} \
        > ${combinedFilePath}
    }
    output {
        File combinedFile = combinedFilePath
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
}