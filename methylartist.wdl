version 1.1

#task Locus {
#	input {
#		Array[File] alignments
#		Array[File] alignmentIndexes
#
#		Region region #= "chr12:49013000-49070000"
#		Region highlight
#
#		String outputPath = "out.png"
#		File annotation_gtf_gz
#		#File annotation_gtf_gz_index
#		File annotation_genome
#		File annotation_genome_index
#
#		Array[String] mods = ["m", "h"]
#
#		Int threads = 2
#		String memory = "8G"
#		Int timeMinutes = 30
#		String dockerImage = "quay.io/biocontainers/methylartist:1.4.0--pyhdfd78af_0"
#	}
#
#	String joined_inputs = sep("\n", alignments)
#
#	String region_fmt = region.chr + ":" + region.start + "-" + region.end
#	String highlight_fmt = highlight.chr + ":" + highlight.start + "-" + highlight.end
#
#	# Use a fasta index file to get the genome sizes. And convert that to the
#	# bedtools specific "genome" format.
#	command {
#		set -e
#		mkdir -p "$(dirname ~{outputPath})"
#
#		echo "~{sep="\n" alignments}" > aln.txt
#
#		methylartist locus \
#			-b aln.txt \
#			-i ~{region_fmt} \
#			-l ~{highlight_fmt} \
#			-g ~{annotation_gtf_gz} \
#			-m ~{sep="," mods} \
#			--ref ~{annotation_genome} \
#			--labelgenes \
#			--motif CG \
#			--outfile "~{outputPath}"
#	}
#		# --skip_align_plot # not available in 1.4.0??
#		# --plot_coverage ~{sep="," alignments} \ # requires samtools
#
#	output {
#		File plot = outputPath
#	}
#
#	runtime {
#		cpu: threads
#		memory: memory
#		time_minutes: timeMinutes
#		docker: dockerImage
#	}
#}
