import pandas as pd
import sqlite3
import argparse
import numpy as np
import os

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--dataset", help = "provide the name of the dataset", action = "store")
parser.add_argument("-m", "--method",  help = "provide the name of the method", action = "store")
args = parser.parse_args()
files_path = "/mnt/pan/Data14/gene_pred/epixcan_results_new/"

ensembl_to_hgnc = pd.read_csv("/mnt/pan/Data14/neural_network/Ensembl_to_hgncsymbol_37.csv")

extra_df = pd.read_csv(files_path + args.dataset + "_" + args.method + "_all_chr_parsed.txt", sep = "\t")
weights_df = pd.read_csv("kunkle_variants_analysis/Kunkle_variants_tfscores.csv")

weights_df = weights_df.dropna()
seed_df = pd.read_csv(files_path + "DGN_seed_all_chr.txt", sep = "\t")

ensembl_to_hgnc.columns = ["gene_name","gene"]
weights_df= weights_df.merge(ensembl_to_hgnc, on = "gene_name", how ="inner")
extra_df.columns = ["gene_name", "n_snps_in_window", "n_snps_in_model", "pred_perf_R2", "pred_perf_pval"]
extra_df =extra_df.merge(ensembl_to_hgnc, on = "gene_name", how = "inner")
extra_df["n_snps_in_model"] = extra_df["n_snps_in_model"].astype("int")
extra_df = extra_df.iloc[:,[5,0,3,2,4]]
weights_df = weights_df[["rsid","gene","weight", "Ref","Alt"]]

def get_column_names_from_db_table(sql_cursor, table_name):
    """
    Scrape the column names from a database table to a list
    :param sql_cursor: sqlite cursor
    :param table_name: table name to get the column names from
    :return: a list with table column names
     """

    table_column_names = 'PRAGMA table_info(' + table_name + ');'
    sql_cursor.execute(table_column_names)
    table_column_names = sql_cursor.fetchall()

    column_names = list()

    for name in table_column_names:
        column_names.append(name[1])

    return column_names

db_file_name = "_".join([args.dataset,args.method,"WB_filtered_tfscores_kunkle.db"])
connection = sqlite3.connect(db_file_name)
cur = connection.cursor()
cur.execute("DROP INDEX IF EXISTS weights_rsid_gene")
cur.execute("DROP INDEX IF EXISTS weights_gene")
cur.execute("DROP INDEX IF EXISTS weights_rsid")
cur.execute('DROP TABLE IF EXISTS weights')
cur.execute('CREATE TABLE weights (rsid TEXT, gene TEXT, weight DOUBLE, ref_allele CHARACTER, eff_allele CHARACTER)')
cur.execute("CREATE INDEX weights_rsid ON weights (rsid)")
cur.execute("CREATE INDEX weights_gene ON weights (gene)")
cur.execute("CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
cur.execute("DROP INDEX IF EXISTS extra_gene")
cur.execute("DROP TABLE IF EXISTS extra")
cur.execute('DROP TABLE IF EXISTS construction')
cur.execute("CREATE TABLE extra (gene TEXT, genename TEXT, " +
            "`pred.perf.R2` DOUBLE, `n.snps.in.model` INTEGER, `pred.perf.pval` DOUBLE)")
cur.execute("CREATE INDEX extra_gene ON extra (gene)")
cur.execute("CREATE TABLE construction (chr INTEGER, " +
            " `cv.seed` INTEGER)")
weights_df.columns = get_column_names_from_db_table(cur, "weights")
extra_df.columns = get_column_names_from_db_table(cur, "extra")
seed_df.columns  =  get_column_names_from_db_table(cur, "construction")
weights_df.to_sql("weights",connection, if_exists='append', index=False)
extra_df.to_sql("extra",connection, if_exists='append', index=False)
seed_df.to_sql("construction",connection, if_exists='append', index=False)
connection.close()


