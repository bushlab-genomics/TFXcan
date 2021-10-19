import pandas as pd
import numpy as np                 
from sklearn.preprocessing import MinMaxScaler,StandardScaler
import multiprocessing as mp
import argparse
import matplotlib as plt 
import random   

files_path  = "/mnt/pan/Data14/gene_pred/tfxcan_git/"
weights_df = pd.read_csv(files_path + "Chrom22_results/Chrom22_variant_tfscores.csv")
geno_df = pd.read_csv(files_path + "Data/Chrom22_genotypes.txt", sep = " ")
def get_simulated_values(gene):
    t_gene_df = weights_df.loc[weights_df.hgnc_symbol == gene][["hgnc_symbol","SNP_ID","cube_root_score_scaled"]]
    t_geno_df = geno_df.loc[geno_df.SNP_ID.isin(list(t_gene_df.SNP_ID))]
    t_geno_df = t_geno_df.merge(t_gene_df, on = ["SNP_ID"], how ="inner")
    t_geno_df= t_geno_df.drop("SNP_ID",axis = 1)
    t_geno_df_t = t_geno_df.T
    tw_vector  = list(t_geno_df_t.loc["cube_root_score_scaled"])
    t_geno_df_t = t_geno_df_t.iloc[:2504,:]
    t_sum = t_geno_df_t.sum(axis = 1) + np.random.normal(0,10,len(t_geno_df_t))
    return(pd.DataFrame(t_sum,columns = [gene]))
pool = mp.Pool()
gene_list = sorted(set(weights_df.hgnc_symbol))
gene_expr_df = pd.concat(list(map(get_simulated_values, gene_list)), axis = 1)
#print(gene_df.head(5))
gene_expr_df.insert(0,"IID", list(gene_expr_df.index))
gene_expr_df.to_csv("Data/Chrom22_simulated_expression.csv", index  =False)

