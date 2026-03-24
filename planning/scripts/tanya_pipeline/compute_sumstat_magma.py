# This script returns mean and spec
# Row is the gene
# Column is the cell type
# Credit: Modified from Rachel's original script single_cell_MAGMA_files_v2.py

import scanpy as sc
import pandas as pd
import argparse
import os
import sys
import anndata
from scipy.sparse import csc_matrix

def main(args):
    """
    :param adata:
    :param ct_colname: this is the column name that define the cell type from the metadata. Typically this could be cell_type, cluster_label, etc...
    :param outdir
    """
    outdir = args.outdir
    adata_original = anndata.read_h5ad(args.h5ad)
    adata_original.var_names_make_unique()
    ct_colname = args.ct_colname

    print("Processing: ", adata_original)

    # check to see if the ct_colname is there
    if not ct_colname in adata_original.obs.columns:
        print("The column name ", ct_colname, " does not exist.")
        exit()

    # because the scRNAseq data contains genes that don't have ensemble, create a filtered adata where genes that were not convertible to ensembl
    # adata_original.var_names_make_unique() #there are no duplicates
    gene_list = adata_original.var[adata_original.var["ensembl"].notnull()].index.values.tolist() #this obtains a list of gene in symbol that has an ensemble conversion
    adata = adata_original[:, gene_list].copy() #subset based on the genes (keep the genes if it was converted to ensemble successfully)

    # now, adata.X is the matrix in numpy format where row is cell ID and column is the gene
    adata_nrow = adata.shape[0]
    adata_ncol = adata.shape[1]
    print(f"adata has {adata_nrow} rows (cells) and {adata_ncol} columns (genes)")

    # first normalize
    # normalize counts per cell (in place)
    sc.pp.normalize_total(adata, target_sum=1e6)

    # define the genes in the dataset
    genes = adata.var["ensembl"].to_list()

    # define the available cell types in the dataset
    cts = adata.obs[ct_colname].dropna().unique()

    # before log-transforming, we determine a filtering on a minimum of 1 transcript per million per celltype
    # start with computing the mean cpM per cell type
    means_cell_counts_pM = pd.DataFrame(data=None, index=cts, columns=genes, dtype=float)  # rows is the cell type; columns is the gene
    for ct in cts:
        means_cell_counts_pM.loc[ct, :] = adata[adata.obs[ct_colname] == ct, :].X.mean(0)

    # the low_filter matrix contains zeros for genes that have < 1 count pM in a certain cell type
    # this matrix will be used to multiply the same data when it is log-transformed
    low_filter = (1 - (means_cell_counts_pM < 1))  # 0 is <1, 1 if >1

    # specificity values
    # set NaNs (i.e. no expression in any of the cells/celltypes) in the mean matrix per cell type to zero
    spec_cell_counts_pM = ((means_cell_counts_pM.mul(low_filter)) / (means_cell_counts_pM.mul(low_filter)).sum(axis=0)).fillna(0)

    # log-transform the original data
    sc.pp.log1p(adata.X, base=2)

    # compute a measure per cell type (mean)
    means_cell_log_counts_pM = pd.DataFrame(data=None, index=cts, columns=genes, dtype=float)

    for ct in cts:
        Y = adata[adata.obs[ct_colname] == ct, :].to_df()
        Y.columns = genes
        means_cell_log_counts_pM.loc[ct, :] = Y.mean(0)

    # filter out the values that have low counts per cell type
    means_cell_log_counts_pM = means_cell_log_counts_pM.mul(low_filter).fillna(0)

    # Format the cell types
    means_cell_log_counts_pM.loc["Average"] = (means_cell_log_counts_pM.mean(axis=0))
    means_cell_log_counts_pM.index = [w.replace(': ', '_') for w in means_cell_log_counts_pM.index.values]
    means_cell_log_counts_pM.index = [w.replace(' ', '_') for w in means_cell_log_counts_pM.index.values]
    means_cell_log_counts_pM.index = [w.replace('/', '_') for w in means_cell_log_counts_pM.index.values]
    means_cell_log_counts_pM.index = [w.replace(':', '_') for w in means_cell_log_counts_pM.index.values]
    means_cell_log_counts_pM.index = [w.replace('-', '_') for w in means_cell_log_counts_pM.index.values]

    spec_cell_counts_pM.index = [w.replace(': ', '_') for w in spec_cell_counts_pM.index.values]
    spec_cell_counts_pM.index = [w.replace(' ', '_') for w in spec_cell_counts_pM.index.values]
    spec_cell_counts_pM.index = [w.replace('/', '_') for w in spec_cell_counts_pM.index.values]
    spec_cell_counts_pM.index = [w.replace(':', '_') for w in spec_cell_counts_pM.index.values]
    spec_cell_counts_pM.index = [w.replace('-', '_') for w in spec_cell_counts_pM.index.values]

    means_cell_log_counts_pM_t = means_cell_log_counts_pM.T
    means_cell_log_counts_pM_t.index.name = "GENE"
    means_cell_log_counts_pM_t.reset_index(inplace=True)
    means_cell_log_counts_pM_out = means_cell_log_counts_pM_t.drop_duplicates(subset=["GENE"])

    spec_cell_counts_pM_t = spec_cell_counts_pM.T
    spec_cell_counts_pM_t.index.name = "GENE"
    spec_cell_counts_pM_t.reset_index(inplace=True)
    spec_cell_counts_pM_out = spec_cell_counts_pM_t.drop_duplicates(subset=["GENE"])

    means_cell_log_counts_pM_out.to_csv(os.path.join(outdir, ct_colname, "means_cell_log_counts_pM.tsv"), sep="\t", index=False) # save to a file)
    spec_cell_counts_pM_out.to_csv(os.path.join(outdir, ct_colname, "spec_cell_log_counts_pM.tsv"), sep="\t", index=False) # save to a file)

def parse_args():
    parser = argparse.ArgumentParser(description='Generate input for magma gene property from scRNAseq data')
    parser.add_argument('--h5ad', required=True, help='Path to the data in h5ad format. This script is to be used with the scRNAseq data after the preprocessing pipeline.')
    parser.add_argument('--ct_colname', required=True, help='Specify the column for the cell type')
    parser.add_argument('--outdir', required=True, help='Path to the output directory')
    return parser.parse_args()


main(parse_args())

