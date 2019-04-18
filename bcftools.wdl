version 1.0

task Bcf2Vcf {
    input {
        File bcf
        String outputPath
    }
    
    command <<<
        bcftools view ~{bcf} -O v -o ~{outputPath}
    >>>
    
    output {
        File OutputVcf = "~{outputPath}"
    }
    
    runtime {
        docker: "quay.io/biocontainers/bcftools:1.9--ha228f0b_3"
    }
}
