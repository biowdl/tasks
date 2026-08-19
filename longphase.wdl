version 1.0
import "whatshap.wdl" as whatshap
import "samtools.wdl" as samtools

task Modcall {
	input {
		File bam
		File bamIndex
		File referenceFasta
		File referenceFastaFai
		String outputPrefix = "modcall_result"

		Int threads = 8
		String memory = "64GiB"
		Int timeMinutes = 30

		# https://github.com/bioconda/bioconda-recipes/pull/62880
		String dockerImage = "quay.io/biocontainers/longphase:2.0--h13024bc_0"
	}

	command <<<
		set -e
		mkdir -p $(dirname ~{outputPrefix})

		longphase modcall \
			-b ~{bam} \
			--all \
			-r ~{referenceFasta} \
			-t ~{threads} \
			-o ~{outputPrefix}
	>>>

	output {
		File vcf = outputPrefix + ".vcf"
	}

	runtime {
		docker: dockerImage
		cpu: threads
		memory: memory
		time_minutes: timeMinutes
	}
}

task Phase {
	input {
		File bam
		File bamIndex
		File referenceFasta
		File referenceFastaFai

		File snpVcf
		File? svVcf
		File modcallVcf

		String outputPrefix = "phased_prefix"

		String seq = "ont" # Or 'pb'

		Int threads = 8
		String memory = "64GiB"
		Int timeMinutes = 15

		# https://github.com/bioconda/bioconda-recipes/pull/62880
		String dockerImage = "quay.io/biocontainers/longphase:2.0--h13024bc_0"
	}

	command <<<
		set -e
		mkdir -p $(dirname ~{outputPrefix})

		longphase phase \
			-s ~{snpVcf} \
			~{"--sv-file " + svVcf} \
			--mod-file ~{modcallVcf} \
			-b ~{bam} \
			-r ~{referenceFasta} \
			-t ~{threads} \
			-o ~{outputPrefix} \
			~{"--" + seq}
	>>>

	output {
		File vcf = outputPrefix + ".vcf"
	}

	runtime {
		docker: dockerImage
		cpu: threads
		memory: memory
		time_minutes: timeMinutes
	}
}

task Haplotag {
	input {
		File bam
		File bamIndex
		File referenceFasta
		File referenceFastaFai

		File phasedSnpVcf

		String outputPrefix = "phased_prefix"

		String seq = "ont" # Or 'pb'

		Int threads = 8
		String memory = "64GiB"
		Int timeMinutes = 30

		# https://github.com/bioconda/bioconda-recipes/pull/62880
		String dockerImage = "quay.io/biocontainers/longphase:2.0--h13024bc_0"
	}

	command <<<
		set -e
		mkdir -p $(dirname ~{outputPrefix})

		longphase haplotag \
			-r ~{referenceFasta} \
			-s ~{phasedSnpVcf} \
			-b ~{bam} \
			-t ~{threads} \
			-o ~{outputPrefix} \
	>>>

	output {
		File outBam = outputPrefix + ".bam"
	}

	runtime {
		docker: dockerImage
		cpu: threads
		memory: memory
		time_minutes: timeMinutes
	}
}

workflow Longphase {
	input {
		File referenceFasta
		File referenceFastaFai
		File bam
		File bamIndex
		File snpVcf
		File? svVcf

		String sample_id
	}

	call Modcall as modcall { input:
		referenceFasta = referenceFasta,
		referenceFastaFai = referenceFastaFai,
		bam = bam,
		bamIndex = bamIndex,
		outputPrefix = "./~{sample_id}.longphase.modcall",
	}

	call Phase as phase { input:
		referenceFasta = referenceFasta,
		referenceFastaFai = referenceFastaFai,
		bam = bam,
		bamIndex = bamIndex,

		snpVcf = snpVcf,
		svVcf = svVcf,
		modcallVcf = modcall.vcf,
		outputPrefix = "./~{sample_id}.longphase.phase",
	}

	call whatshap.Stats as whatshapStats { input:
		vcf = phase.vcf,
		# Out
		gtf = "~{sample_id}.longphase.stats.gtf",
		tsv = "~{sample_id}.longphase.stats.tsv",
		blockList = "~{sample_id}.longphase.stats.blocklist"
	}

	call Haplotag as haplotag { input:
		referenceFasta = referenceFasta,
		referenceFastaFai = referenceFastaFai,
		bam = bam,
		bamIndex = bamIndex,
		phasedSnpVcf = phase.vcf,
		outputPrefix = "./~{sample_id}.longphase.haplotagged",
	}

	call samtools.Index as index_bam { input:
		bamFile = haplotag.outBam,
	}

	output {
		File mod_5mC_vcf = modcall.vcf
		File phased_vcf = phase.vcf
		File haplotagged_bam = index_bam.indexedBam
		File haplotagged_bam_index = index_bam.index
		File? phasedGTF = whatshapStats.phasedGTF
		File? phasedTSV = whatshapStats.phasedTSV
		File? phasedBlockList = whatshapStats.phasedBlockList
	}
}
