#TO DO
## Requires parameter_meta:
* biopet.wdl: `ExtractAdaptersFastqc`

## Duplicate tasks:
* 

## Out of date with new cluster & parameter_meta:
* bamstats.wdl: `Generate`
* biopet.wdl: `BaseCounter`, `FastqSplitter`, `FastqSync`,
              `ValidateAnnotation`, `ValidateFastq`, `ValidateVcf`, `VcfStats`
* sampleconfig.wdl: `SampleConfig`, `SampleConfigCromwellArrays`, `CaseControl`
* seqstat.wdl: `Generate`
* common.wdl: `AppendToStringArray`, `CheckFileMD5`, `ConcatenateTextFiles`,
              `Copy`, `CreateLink`, `MapMd5`, `StringArrayMd5`

## Imports other tasks:
* bamstats.wdl
* biopet.wdl
* sampleconfig.wdl
* seqstat.wdl
* clever.wdl
