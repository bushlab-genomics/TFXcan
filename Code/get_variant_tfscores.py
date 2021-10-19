import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
import multiprocessing as mp
import scipy.stats
import math
import os
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-i","--input", help="please provide the file name containing variant TFBS TFAGNet scores", action = "store")
parser.add_argument("-o", "--output", help = "please provide the prefix for the output file name",action = "store")

args = parser.parse_args()
out_path = args.output + "_results"

if not os.path.exists(out_path):
	os.mkdir(out_path)

def get_cube_root(x):
    if x>0:
        ans = x**(1./3.)
    else:
        ans = -((-x)**(1./3.))
    return(ans)

res_df = pd.read_csv(args.input)
coef_tf = pd.read_csv("Data/GM12878_tf_coefs.csv")
coef_tf.columns = ["TF_Name", "Enet_Coef"]
ens_hgnc = pd.read_csv("Data/Ensembl_to_hgncsymbol_37.csv")
res_df = res_df.merge(ens_hgnc, on = "Ensembl", how = "inner")
res_df1 = res_df[["TF_Name","hgnc_symbol", "SNP_ID","Log_score"]]
res_df1.columns = ["TF_Name", "hgnc_symbol", "SNP_ID", "Log_score"]

res_df1 = res_df1.drop_duplicates()
o_name = out_path + "/" +  args.output + "_variant_tfscores.csv"
res_coef_df  = res_df1.merge(coef_tf, on = "TF_Name", how = "inner")
res_coef_df["merge_score"] = res_coef_df["Log_score"]*res_coef_df["Enet_Coef"]
res_gn_df = res_coef_df[["hgnc_symbol", "SNP_ID", "merge_score"]]
res_gn_df = res_gn_df.drop_duplicates()
res_panda_gn_df1 = res_gn_df.groupby(["hgnc_symbol", "SNP_ID"]).sum()
ind_list  = list(res_panda_gn_df1.index)
gene_list= [ind_list[i][0] for i in range(len(ind_list))]
snp_id = [ind_list[i][1] for i in range(len(ind_list))]
res_panda_gn_df1.insert(0,"hgnc_symbol", gene_list)
res_panda_gn_df1.insert(1,"SNP_ID", snp_id)
res_panda_gn_df1 = res_panda_gn_df1.reset_index(drop= True)
res_panda_gn_df1["merge_score"] = res_panda_gn_df1["merge_score"]
pool = mp.Pool(20)
res_panda_gn_df1["cube_root_score"] = list(pool.map(get_cube_root,list(np.abs(res_panda_gn_df1["merge_score"]))))
res_panda_gn_df1["cube_root_score_scaled"] =  MinMaxScaler((1,10)).fit_transform(np.array(res_panda_gn_df1["cube_root_score"]).reshape(-1,1)).reshape(-1,)	
res_panda_gn_df1.to_csv(o_name, index = False)




