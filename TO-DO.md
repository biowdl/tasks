#TO DO
This file describes WDL files and tasks within those files which need
more specific attention than just adding outputs to the parameter_meta.

Some tasks have not been updated to match the new SLURM requirements and are
missing a parameter_meta section.

Some tasks are importing other WDL files.

## Out of date with new cluster & parameter_meta:
* bamstats.wdl: `Generate`
* biopet.wdl: `BaseCounter`, `FastqSplitter`, `FastqSync`,
              `ValidateAnnotation`, `ValidateFastq`, `ValidateVcf`, `VcfStats`
* sampleconfig.wdl: `SampleConfig`, `SampleConfigCromwellArrays`, `CaseControl`
* seqstat.wdl: `Generate`
* common.wdl: `AppendToStringArray`, `CheckFileMD5`, `ConcatenateTextFiles`,
              `Copy`, `CreateLink`, `MapMd5`, `StringArrayMd5`
* fastqsplitter.wdl: `Fastqsplitter`
* flash.wdl: `Flash`
* macs2.wdl: `PeakCalling`
* ncbi.wdl: `GenomeDownload`, `DownloadNtFasta`, `DownloadAccessionToTaxId`
* seqtk.wdl: `Sample`
* spades.wdl: `Spades`
* unicycler.wdl: `Unicycler`


## Imports other tasks:
* bamstats.wdl
* biopet.wdl
* sampleconfig.wdl
* seqstat.wdl
* clever.wdl
* strelka.wdl

## Requires input from others:
* somaticseq.wdl
