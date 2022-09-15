# transcript_annotation_2022

## main_script folder

### script_process_RNAseq.sh

Bash script with command lines used for RNAseq data processing: alignment, counting, stringtie, TopHat, etc.

## SCV

### process_SCV.R
### analysis_SCV_GL.R
These two scripts are used to compute transcript annotation using the SCV method and analyse results to get annotations and data.
Outputs can be found in the UTR_annotation directory.

## annotations

Directory containing the transcript annotation file

## TopHat

Directory containing the scripts used to analyse the TopHat output for alternative splicing detection.
The tables containing all results are also included

## IR_results

Directory containing the scripts used to detecte intron retention

## NATs

Directory containing the scripts used to detect antisens transcripts.