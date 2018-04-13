If no normal samples are available there may be some samples in the tumour cohort that can be used instead. These would either be cancer samples with structurally and numerically normal genomes, or cancer samples that turn out to be of such extremely low cellularity that they are effectively normal tissues samples.

To identify whether there are suitable samples in the tumour cohort run the script prelimPlots.R for each sample. This takes as input the following arguments (in order): a sample bam file; a headerless bed file containing the first four columns of the input bed file described in the README in the main repository; the sample name.

prelimPlots.R outputs a genome-wide copy number plot for each sample and a text file containing the mean,  standard deviation, maximum and minimum of the autosomal data points. Note that the copy number data will be noisy at this stage as the variability due to different capture efficiencies across the panel has not normalised out.

Select samples with a mean value that is representative of the cohort (in particular avoid samples with a low mean, which would indicate the sample has been under-sequenced) and that have a low standard deviation (which indicates low noise and lower chance that the sample contains copy number alterations). The maximum and minium values can indiciate if there are focal amplifications or deletions in the sample. Review the copy number plots and metrics to select the most suitable samples to make up a pseudo-normal panel. If the bed file contains regions on chrX and chrY then a mixture of male and female samples are needed in the normal panel 

Note: If you select samples with copy number alterations for the pseudo-normal panel there is a risk of low sensitivity and false positive calls.
