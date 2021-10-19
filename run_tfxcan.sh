#python Code/get_variant_tfscores.py -i Data/Chrom22_tfbs_scores_postprocessed.csv -o Chrom22
Rscript Code/run_tfxcan.R Chrom22 10 0.5 Data/Chrom22_expression.RDS Data/Chrom22_genotypes.txt Data/Chrom22_snp_annotations.RDS Data/Chrom22_gene_annot.RDS 22
