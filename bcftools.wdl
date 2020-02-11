version 1.0

task Bcf2Vcf {
    input {
        File bcf
        String outputPath
        String dockerImage = "quay.io/biocontainers/bcftools:1.9--ha228f0b_3"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        bcftools view ~{bcf} -O v -o ~{outputPath}
    }

    output {
        File outputVcf = outputPath
    }

    runtime {
        docker: dockerImage
    }

    parameter_meta {
        bcf: {description: "The generated BCF from an SV caller", category: "advanced"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}
