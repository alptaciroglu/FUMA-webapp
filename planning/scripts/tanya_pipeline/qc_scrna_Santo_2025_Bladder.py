# This script processes the scRNAseq data from Santo et al. 2025 (GUDMAP consortium)
# Human urinary bladder — 3 CellxGene h5ad files (urothelial, stromal, immune)
# Author: Alp Taciroglu
# Date: 2026-03-24
# Follows Tanya Phung's preprocessing pipeline (github.com/tanyaphung/scrnaseq_viewer)
#
# Source: CellxGene collection 0e54d4de-44f0-4d50-8649-b5c2bbe8f5d1
# Publication: Santo et al. (2025) iScience, DOI: 10.1016/j.isci.2024.111628
#
# Notes on this dataset:
# - CellxGene h5ad files have Ensembl IDs as var index (no symbol-to-Ensembl conversion needed)
# - adata.X contains normalised/log-transformed values; raw integer counts are in adata.raw.X
# - 3 subsets: Urothelial (14776 cells, 3 types), Stromal (12964 cells, 12 types),
#              Immune (2094 cells, 9 types) — merged into one dataset
# - Single cell type hierarchy level (author_cell_type)
#
# Usage:
#   python qc_scrna_Santo_2025_Bladder.py <base_dir> <id>
#   Example: python qc_scrna_Santo_2025_Bladder.py planning/data/GUDMAP Santo_Human_2025_Bladder

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import os
import sys
from scipy.sparse import csr_matrix

def main():
    # setting up
    base_dir = sys.argv[1]
    id = sys.argv[2]

    # Three h5ad files from CellxGene (Santo et al. 2025, human bladder)
    h5ad_files = [
        os.path.join(base_dir, "e8799823-156b-4d74-b887-2305b86635eb.h5ad"),  # Urothelial
        os.path.join(base_dir, "4feb85d6-7a8e-4d7b-9b3b-5c7aa353a542.h5ad"),  # Stromal
        os.path.join(base_dir, "cdf3c8cc-9eb3-4e5c-9d49-997552ffe55e.h5ad"),  # Immune
    ]

    log = open(os.path.join(base_dir, id, "log.txt"), "w")
    clean_adata_fp = os.path.join(base_dir, id, id + ".h5ad")

    # format table 2
    table2 = open(os.path.join(base_dir, id, "table2.csv"), "w")
    header = ["Description", "Number of cells", "Number of genes"]
    print(",".join(header), file=table2)

    # format table 3
    table3 = open(os.path.join(base_dir, id, "table3.csv"), "w")

    # read and merge h5ad files
    adatas = []
    # Columns to keep in var (shared across all three files)
    keep_var_cols = ["feature_name", "feature_biotype", "feature_type"]

    for fp in h5ad_files:
        print(f"Loading {fp}...", file=log)
        a = anndata.read_h5ad(fp)
        print(f"  {a.n_obs} cells, {a.n_vars} genes", file=log)

        # Use raw counts — adata.X is already normalised in CellxGene files
        print("  Replacing X with raw counts from adata.raw.X", file=log)
        a.X = a.raw.X.copy()

        # Keep only shared var columns so concat works
        a.var = a.var[keep_var_cols]

        adatas.append(a)

    # Save var metadata from first file before concat (concat drops var columns)
    var_metadata = adatas[0].var.copy()

    # Merge (inner join on genes — all three files have the same 35459 genes)
    adata = anndata.concat(adatas, join="inner")
    adata.var_names_make_unique()

    # Restore var metadata
    adata.var = var_metadata.loc[adata.var.index]

    print(f"\nMerged: {adata.n_obs} cells, {adata.n_vars} genes", file=log)

    print("Viewing the adata observations.", file=log)
    print(adata.obs, file=log)

    print("Viewing the adata variables.", file=log)
    print(adata.var, file=log)

    print("Viewing the adata matrix - are these integer counts?", file=log)
    sample = adata.X[:25, :25]
    if hasattr(sample, 'toarray'):
        sample = sample.toarray()
    print(sample, file=log)

    adata_nrow = adata.shape[0]
    adata_ncol = adata.shape[1]
    print(f"adata has {adata_nrow} cells and {adata_ncol} genes.", file=log)
    original = ["Original", str(adata.n_obs), str(adata.n_vars)]
    print(",".join(original), file=table2)

    # mitochondrial genes — CellxGene var has feature_name with gene symbols
    adata.var['mt'] = adata.var["feature_name"].str.startswith("MT-")
    n_mt_genes = np.count_nonzero(adata.var['mt'])
    print(f"adata has {n_mt_genes} mitochondrial genes.", file=log)

    # calculate qc metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # first filtering: accept cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    print(f"After first filtering (see documentation for definition), adata has {adata.n_obs} cells and {adata.n_vars} genes.", file=log)
    first_filter = ["Filter#1", str(adata.n_obs), str(adata.n_vars)]
    print(",".join(first_filter), file=table2)

    # second filtering based on mt percentage
    adata = adata[adata.obs['pct_counts_mt'] < 10, :].copy()
    print(f"After second filtering (see documentation for definition), adata has {adata.n_obs} cells and {adata.n_vars} genes.", file=log)
    second_filter = ["Filter#2", str(adata.n_obs), str(adata.n_vars)]
    print(",".join(second_filter), file=table2)

    # rename cell type columns — this dataset has a single level (author_cell_type)
    adata.obs.rename(columns={"author_cell_type": "cell_type_level_1"}, inplace=True)

    # Gene ID handling:
    # CellxGene var index is already Ensembl IDs (ENSG...), so no symbol-to-Ensembl conversion needed.
    # However, we filter to genes present in gene_names_human.txt to match Tanya's pipeline.
    # In Tanya's standard pipeline, the symbol-to-Ensembl conversion via gene_names_human.txt
    # implicitly filters to known genes. Since CellxGene data already has Ensembl IDs, we
    # replicate this filtering step explicitly to keep the gene set consistent with other
    # FUMA datasets (~20-25k protein-coding + well-annotated genes, rather than ~31k with lncRNAs).
    gene_names_fp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gene_names_human.txt")
    gene_names = pd.read_csv(gene_names_fp, sep="\t")
    known_ensembl = set(gene_names["Ensembl gene ID"].dropna())
    gene_mask = adata.var.index.isin(known_ensembl)
    n_before = adata.n_vars
    adata = adata[:, gene_mask].copy()
    print(f"Filtered to genes in gene_names_human.txt: {n_before} -> {adata.n_vars} genes.", file=log)

    adata.var["ensembl"] = adata.var.index
    ngenes_with_ensembl = adata.var["ensembl"].dropna().shape[0]
    print(f"Ensembl IDs in var index: {ngenes_with_ensembl} genes.", file=log)
    ensembl_convertable = ["Ngenes with ensembl", str(adata.n_obs), str(ngenes_with_ensembl)]
    print(",".join(ensembl_convertable), file=table2)

    # save
    adata_out = adata.copy()
    adata_out.X = csr_matrix(adata_out.X)  # convert to sparse
    adata_out.write_h5ad(filename=clean_adata_fp)

    # print out the saved adata
    print(adata_out, file=log)
    print(adata_out.var, file=log)
    print(adata_out.obs, file=log)

    # tabulate the number of cells per cell type
    level1_cts = adata.obs["cell_type_level_1"].dropna().unique()
    for ct in level1_cts:
        ct_data = adata[adata.obs["cell_type_level_1"] == ct, :].X
        out = ["level_1", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3)

    table2.close()
    table3.close()
    log.close()

    # Also print summary to stdout
    print(f"QC complete: {adata.n_obs} cells, {adata.n_vars} genes, {len(level1_cts)} cell types")
    print(f"Output: {clean_adata_fp}")


main()
