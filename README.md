# TFXCan
This page contains code and tutorial for running **TFXCan** models, which use transcription driven gene regulatory scores for weighting genetic variants in a TWAS framework. 

## Requirements
### python(>= 3.6.6)
* numpy(v 1.19.5)
* pandas(v. 1.1.5)
* sci-kit learn(v. 0.22.2)

### R(>= 3.6.0)
* dplyr(0.8.5)
* glmnet(4.0.2)
* reshape2(1.4.4)

### Description
There are two scripts in the **Code** folder. The **get_variant_tfscores.py** generates regulatory TFScores for a set of cis-variants for each gene, while the **run_tfxcan.R** trains TFXcan models using genotype and expression data along with these TFScores. Further details regarding how to run these scripts are provided below:

1) The **get_variant_tfscores.py** generates regulatory tfscores for variants using TFAGNet scores reflecting the influence of variants on TF binding[1] and TF effect estimates representing individual influence of TFs on gene regulation[2]

```
python Code/get_variant_tfscores.py -i $variant_tfbscores -o $output_prefix 

```
Here, **$variant_tfbscores** file contains the influence of genetic variants on TF binding sites, which can be obtained by using our docker container(https://github.com/bushlab-genomics/TFAGNet) for a set of variants. The **$output_prefix** should correspond to the prefix for the output file containing variant tfscores, which will be created in the **$output_prefix_results** folder. 

2) The **run_tfxcan.R** builds and trains the TFXcan models.

```
Rscript Code/run_tfxcan.R $output_prefix $number_of_k-folds $alpha $expression_rds $genotype_file $snp_annotation_rds $gene_annotation_rds  $chromosome_name
```
Here, **$output_prefix** is the same as the one used above for the **get_variant_tfscores.py** script. Thus, the results for the TFXcan models will also be stored in the **$output_prefix_results** folder. The **$number_of_k-folds** and **alpha** represent the number of k-folds used for cross-validation and for calculation of nested cross-validation R2 and p-value, and the alpha parameter used for training the elastic net models respectively. Please refer to the PrediXcan pipeline(https://github.com/hakyimlab/PredictDB-Tutorial) for further details. Additionally, **$expression_rds** should contain the name of the expression file for a dataset saved in the **.RDS** format(you can use **saveRDS** function of R to generate this file). **genotype_file** should contain the genotype information for a set of individuals for a given dataset. **snp_annotation_rds** and **gene_annotation_rds** represent **.RDS** format file names for variant and gene annotations respectively(refer to the PrediXcan pipeline for information on how to generate these files). Lastly, **chromosome_name** should reflect the name of the chromosome. The script iterates over all the genes for the given chromosome within the dataset and produces three files in the **$output_prefix_results** folder: 1) **results file** containing the model accuracy metrics such as R2, p-value, number of snps used in the model etc. for each gene. 2) **covariance file** containing the covariance calculated for each pair of variants for each gene model. 3) **weights file** containing the beta for each variant learned from each gene model. These three files can be used to generate the **.db** files, which can subsequently used to generate S-Predixcan(https://github.com/hakyimlab/MetaXcan). 

### Quick Run
The **Data** folder contains simulated expression, genotype and other files necessary for running TFXcan. You can use the following shell script to quickly build tfxcan models using these files:
```
. run_tfxcan.sh
```



