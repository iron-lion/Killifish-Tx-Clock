Quality Control checking
=============================================================================================

01_CompareCountMatrices_AcrossBatch.R -  this script checks that the merging of read counts across sequencing lanes worked

02_merge_counts_TPM.R - this script performs read count merging for all batches, TPM merging, and modifies the ExperimentDesign table to censor individuals.

03_checkreads.R - this script inspects read count summaries to remove low count samples

03a_readcountsummary.R - this script does a very simple analysis to determine number of samples with read count > 30M PE reads

04_featureCounts.R - this script inspects mapping summaries to remove low mapped samples

05_checksampleidentity.R - this script inspects sample type and sex, to check identities 

06_outlierdetection.R - performs sample connectivity analysis to detect outliers

07_VariancePartition.R - performs variancePartition analysis

ExtDataFig1b_tissuebreakdown.R - generates a table of sample breakdown to use for filling out ExtDataFig1b
 

Directory Structure:
==========================

------------------------------+
Input   |
------------------------------+
This is where all the data input will be added (like raw counts, TPM, experiment design table, multiqc outputs, etc.)

-----------------+
Output          |
-----------------+
This is where all of the output tables and some plots will be stored.