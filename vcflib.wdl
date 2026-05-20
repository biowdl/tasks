version 1.1

task Vcfuniqalleles{
	input {
		File vcf

		String outputPath = "map.vcf"
		Boolean refConsistent = false

		Int threads = 1
		String memory = "12G"
		Int timeMinutes = 60
		String dockerImage = "quay.io/biocontainers/vcflib:1.0.14--h34261f4_0"
	}

	command {
		set -ex
		mkdir -p "$(dirname ~{outputPath})"

		vcfuniqalleles "~{vcf}" > "~{outputPath}"
	}

	output {
		File out = outputPath
	}

	runtime {
		cpu: threads
		memory: memory
		time_minutes: timeMinutes
		container: dockerImage
	}
}
