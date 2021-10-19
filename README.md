# TFXCan
This page contains code and tutorial for running *TFXCan* models, which use transcription driven gene regulatory scores for weighting genetic variants in a TWAS framework. 

Requirements
software
python(>= 3.6.6)
numpy(v 1.19.5)
pandas(v. 1.1.5)
sci-kit learn(v. 0.22.2)

R(>= 3.6.0)
dplyr(0.8.5)
glmnet(4.0.2)
reshape2(1.4.4)

There are two scripts in the *Code* folder. The *get_variant_tfscores.py* generates regulatory TFScores for a set of cis-variants for each gene, while the *run_tfxcan.R* trains TFXcan models using genotype and expression data along with these TFScores. 





