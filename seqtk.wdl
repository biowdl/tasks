version 1.0

task Sample {
    input {
        File sequenceFile
        String outFilePath = "subsampledReads.fq.gz"
        String? preCommand
        Int? seed
        Boolean twoPassMode = false
        Float? fractionOrNumber # when above 1.0 is the number of reads, otherwise it's a fraction
        Boolean zip = true
    }

    command {
        set -e -o pipefail
        mkdir -p $(dirname outFilePath)
        ~{preCommand}
        seqtk sample \
        ~{"-s " + seed} \
        ~{true="-2 " false="" twoPassMode} \
        ~{sequenceFile} \
        ~{fractionOrNumber} \
        ~{true="| gzip" false="" zip} \
        >  ~{outFilePath}
    }

    output {
        File subsampledReads = outFilePath
    }
}