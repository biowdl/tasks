version 1.1


task Liftover {
	input {
		File vcf
		File chain
		File referenceFasta

		String outputPath = "map.vcf"

		Int threads = 1
		String memory = "12G"
		Int timeMinutes = 60
		String dockerImage = "quay.io/biocontainers/ucsc-liftover:482--h0b57e2e_0"
	}

	command {
		set -e
		mkdir -p "$(dirname ~{outputPath})"

		CrossMap vcf \
			~{chain} \
			~{vcf} \
			~{referenceFasta} \
			~{outputPath}
	}

	output {
		File out = outputPath
		File unmapped = outputPath + ".unmap"
	}

	runtime {
		cpu: threads
		memory: memory
		time_minutes: timeMinutes
		container: dockerImage
	}
}
