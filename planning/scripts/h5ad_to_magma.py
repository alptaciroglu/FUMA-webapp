#!/usr/bin/env python3
"""
Convert h5ad files to MAGMA cell type expression format for FUMA.

This script produces the same output as Tanya Phung's compute_sumstat_magma.py
(github.com/tanyaphung/scrnaseq_viewer/code/postprocessing/compute_sumstat_magma.py)
but differs in the following ways:

  Differences from Tanya's pipeline:
  - Combined into a single script with two subcommands (qc, magma), whereas Tanya's
    pipeline has separate per-dataset QC scripts + a shared compute_sumstat_magma.py.
  - Optimised for CellxGene h5ad files which already have Ensembl IDs in the var index,
    so gene symbol-to-Ensembl conversion (via gene_names_human.txt) is skipped.
  - For non-CellxGene data or data without Ensembl IDs, use Tanya's helper functions:
    single_cell_helper_functions_v3.py -> add_gene_names_human() / add_gene_names_mouse()
    (requires gene_names_human.txt from genenames.org)
  - Does not handle mouse-to-human homolog conversion (Tanya's pipeline does via MGI IDs).

  What is identical to Tanya's pipeline:
  - QC thresholds: min 200 genes/cell, min 3 cells/gene, MT% < 10%
  - CPM normalisation (target_sum=1e6)
  - Low-expression filter (<1 CPM zeroed out)
  - Log2(x+1) transform
  - Per-cell-type mean computation (pre-log and post-log)
  - Specificity calculation (proportion of total expression per gene)
  - Output format: two TSVs (means + specificity), same column structure

  Scope: standard reference scRNA-seq data only. Case/control data (e.g. disease vs
  healthy donors) requires a different MAGMA analysis approach — see planning docs.

Two-stage process:
  Stage 1 (this script, QC mode): Load h5ad, QC filter, convert gene IDs, save clean h5ad
  Stage 2 (this script, MAGMA mode): CPM normalise, log2 transform, compute per-cell-type
          means, apply low-expression filter, output MAGMA-format TSV + specificity TSV

Output format: tab-delimited .tsv with columns:
  GENE | CellType_1 | CellType_2 | ... | Average

Usage:
  # Stage 1: QC and gene conversion (for CellxGene h5ad with Ensembl IDs already present)
  python3 h5ad_to_magma.py qc \
    --input file1.h5ad file2.h5ad ... \
    --output clean_merged.h5ad \
    --cell_type_col author_cell_type

  # Stage 2: Generate MAGMA-format output from clean h5ad
  python3 h5ad_to_magma.py magma \
    --input clean_merged.h5ad \
    --outdir output_directory/ \
    --cell_type_col author_cell_type \
    --dataset_name Santo_Human_2025_Bladder_level2
"""

import argparse
import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp


# =============================================================================
# Stage 1: QC and preprocessing
# =============================================================================

def load_and_merge(input_files):
    """Load h5ad files and concatenate into a single AnnData object."""
    adatas = []
    for f in input_files:
        print(f"Loading {f}...")
        adata = ad.read_h5ad(f)
        adatas.append(adata)
        print(f"  {adata.n_obs} cells, {adata.n_vars} genes")

    if len(adatas) == 1:
        merged = adatas[0]
    else:
        merged = ad.concat(adatas, join="inner")

    merged.var_names_make_unique()
    print(f"\nMerged: {merged.n_obs} cells, {merged.n_vars} genes")
    return merged


def qc_filter(adata):
    """Apply standard QC filters following Tanya's preprocessing pipeline.

    - Filter cells with < 200 detected genes
    - Filter genes expressed in < 3 cells
    - Remove cells with > 10% mitochondrial counts
    """
    print(f"\nBefore QC: {adata.n_obs} cells, {adata.n_vars} genes")

    # Mitochondrial genes
    if "feature_name" in adata.var.columns:
        adata.var["mt"] = adata.var["feature_name"].str.startswith("MT-")
    else:
        # CellxGene format may use var index as gene symbols or Ensembl IDs
        # Try to detect MT genes from gene names
        adata.var["mt"] = adata.var.index.str.startswith("MT-")

    n_mt = adata.var["mt"].sum()
    print(f"  Mitochondrial genes: {n_mt}")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Filter 1: min genes per cell, min cells per gene
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"After filter 1 (min 200 genes/cell, min 3 cells/gene): {adata.n_obs} cells, {adata.n_vars} genes")

    # Filter 2: mitochondrial percentage
    adata = adata[adata.obs["pct_counts_mt"] < 10, :].copy()
    print(f"After filter 2 (MT < 10%): {adata.n_obs} cells, {adata.n_vars} genes")

    return adata


def resolve_ensembl_ids(adata):
    """Ensure adata has an 'ensembl' column in var for downstream processing.

    CellxGene h5ad files typically have Ensembl IDs as the var index.
    Other sources may need conversion via gene_names_human.txt (see Tanya's helper functions).
    """
    # Check if var index is already Ensembl IDs
    sample_ids = adata.var.index[:10].tolist()
    ensembl_in_index = all(str(g).startswith("ENSG") for g in sample_ids if pd.notna(g))

    if ensembl_in_index:
        adata.var["ensembl"] = adata.var.index
        n_ensembl = adata.var["ensembl"].notna().sum()
        print(f"\nGene IDs: Ensembl IDs found in var index ({n_ensembl} genes)")
    elif "ensembl" in adata.var.columns:
        n_ensembl = adata.var["ensembl"].notna().sum()
        print(f"\nGene IDs: 'ensembl' column already present ({n_ensembl} genes)")
    else:
        print("\nWARNING: No Ensembl IDs found in var index or 'ensembl' column.")
        print("For non-CellxGene data, use Tanya's helper functions for gene conversion:")
        print("  github.com/tanyaphung/scrnaseq_viewer/code/preprocessing/single_cell_helper_functions_v3.py")
        print("  Requires: gene_names_human.txt from genenames.org")
        sys.exit(1)

    # Filter to genes with valid Ensembl IDs
    before = adata.n_vars
    gene_mask = adata.var["ensembl"].notna() & adata.var["ensembl"].str.startswith("ENSG")
    adata = adata[:, gene_mask].copy()
    print(f"Filtered to valid Ensembl IDs: {before} -> {adata.n_vars} genes")

    return adata


def run_qc(args):
    """Stage 1: QC, filter, gene conversion, save clean h5ad."""
    adata = load_and_merge(args.input)

    # Validate cell type column exists
    if args.cell_type_col not in adata.obs.columns:
        available = [c for c in adata.obs.columns if "cell" in c.lower() or "type" in c.lower()]
        print(f"ERROR: '{args.cell_type_col}' not found. Candidates: {available}")
        sys.exit(1)

    # Report cell types
    cts = sorted(adata.obs[args.cell_type_col].dropna().unique())
    print(f"\nCell types in '{args.cell_type_col}' ({len(cts)}):")
    for ct in cts:
        n = (adata.obs[args.cell_type_col] == ct).sum()
        print(f"  {ct}: {n} cells")

    # QC filtering
    adata = qc_filter(adata)

    # Resolve Ensembl gene IDs
    adata = resolve_ensembl_ids(adata)

    # Save
    adata.write_h5ad(args.output)
    print(f"\nClean h5ad saved to: {args.output}")
    print(f"Final: {adata.n_obs} cells, {adata.n_vars} genes, {len(cts)} cell types")


# =============================================================================
# Stage 2: MAGMA format generation
# =============================================================================

def sanitise_names(names):
    """Replace characters that MAGMA can't handle in variable names.
    Follows Tanya's convention: replace ': ', ' ', '/', ':', '-' with '_'.
    """
    sanitised = []
    for w in names:
        w = w.replace(": ", "_")
        w = w.replace(" ", "_")
        w = w.replace("/", "_")
        w = w.replace(":", "_")
        w = w.replace("-", "_")
        sanitised.append(w)
    return sanitised


def run_magma(args):
    """Stage 2: Generate MAGMA-format output from preprocessed h5ad.

    Follows compute_sumstat_magma.py from Tanya Phung's pipeline:
    1. Filter to genes with Ensembl IDs
    2. CPM normalise (target_sum=1e6)
    3. Compute mean CPM per cell type (pre-log)
    4. Apply low-expression filter (<1 CPM zeroed out)
    5. Log2(x+1) transform
    6. Compute mean log2-CPM per cell type
    7. Apply low-expression filter to log means
    8. Compute specificity values
    9. Add Average row, sanitise names, output TSV files
    """
    print(f"Loading {args.input}...")
    adata_original = ad.read_h5ad(args.input)
    adata_original.var_names_make_unique()
    ct_colname = args.cell_type_col

    if ct_colname not in adata_original.obs.columns:
        available = [c for c in adata_original.obs.columns if "cell" in c.lower() or "type" in c.lower()]
        print(f"ERROR: '{ct_colname}' not found. Candidates: {available}")
        sys.exit(1)

    # Filter to genes with Ensembl IDs
    if "ensembl" not in adata_original.var.columns:
        # CellxGene: Ensembl IDs are the var index
        if adata_original.var.index[0].startswith("ENSG"):
            adata_original.var["ensembl"] = adata_original.var.index
        else:
            print("ERROR: No 'ensembl' column found. Run 'qc' stage first.")
            sys.exit(1)

    gene_list = adata_original.var[adata_original.var["ensembl"].notna()].index.tolist()
    adata = adata_original[:, gene_list].copy()

    print(f"adata: {adata.n_obs} cells, {adata.n_vars} genes")

    # CPM normalise
    sc.pp.normalize_total(adata, target_sum=1e6)

    # Gene IDs for output
    genes = adata.var["ensembl"].tolist()

    # Cell types
    cts = adata.obs[ct_colname].dropna().unique()
    print(f"\nCell types ({len(cts)}):")
    for ct in sorted(cts):
        print(f"  {ct}: {(adata.obs[ct_colname] == ct).sum()} cells")

    # Compute mean CPM per cell type (before log transform)
    means_cell_counts_pM = pd.DataFrame(data=None, index=cts, columns=genes, dtype=float)
    for ct in cts:
        subset = adata[adata.obs[ct_colname] == ct, :]
        if sp.issparse(subset.X):
            means_cell_counts_pM.loc[ct, :] = np.asarray(subset.X.mean(axis=0)).flatten()
        else:
            means_cell_counts_pM.loc[ct, :] = subset.X.mean(axis=0).flatten()

    # Low-expression filter: zero out genes with < 1 CPM per cell type
    low_filter = (1 - (means_cell_counts_pM < 1))  # 0 if <1 CPM, 1 if >=1 CPM

    # Specificity values (proportion of total expression per gene)
    filtered_means = means_cell_counts_pM.mul(low_filter)
    spec_cell_counts_pM = (filtered_means / filtered_means.sum(axis=0)).fillna(0)

    # Log2(x+1) transform the normalised data
    sc.pp.log1p(adata.X, base=2)

    # Compute mean log2-CPM per cell type
    means_cell_log_counts_pM = pd.DataFrame(data=None, index=cts, columns=genes, dtype=float)
    for ct in cts:
        Y = adata[adata.obs[ct_colname] == ct, :].to_df()
        Y.columns = genes
        means_cell_log_counts_pM.loc[ct, :] = Y.mean(0)

    # Apply low-expression filter to log means
    means_cell_log_counts_pM = means_cell_log_counts_pM.mul(low_filter).fillna(0)

    # Add Average row
    means_cell_log_counts_pM.loc["Average"] = means_cell_log_counts_pM.mean(axis=0)

    # Sanitise cell type names
    means_cell_log_counts_pM.index = sanitise_names(means_cell_log_counts_pM.index.values)
    spec_cell_counts_pM.index = sanitise_names(spec_cell_counts_pM.index.values)

    # Transpose: genes as rows, cell types as columns
    means_t = means_cell_log_counts_pM.T
    means_t.index.name = "GENE"
    means_t.reset_index(inplace=True)
    means_out = means_t.drop_duplicates(subset=["GENE"])

    spec_t = spec_cell_counts_pM.T
    spec_t.index.name = "GENE"
    spec_t.reset_index(inplace=True)
    spec_out = spec_t.drop_duplicates(subset=["GENE"])

    # Create output directory
    outdir = os.path.join(args.outdir, ct_colname)
    os.makedirs(outdir, exist_ok=True)

    # Write output files
    means_fp = os.path.join(outdir, "means_cell_log_counts_pM.tsv")
    spec_fp = os.path.join(outdir, "spec_cell_log_counts_pM.tsv")
    means_out.to_csv(means_fp, sep="\t", index=False)
    spec_out.to_csv(spec_fp, sep="\t", index=False)

    # QC report
    ct_cols = [c for c in means_out.columns if c != "GENE" and c != "Average"]
    print("\n" + "=" * 60)
    print(f"QC REPORT: {args.dataset_name}")
    print("=" * 60)
    print(f"Total genes: {len(means_out)}")
    print(f"Cell types: {len(ct_cols)}")
    print(f"Cell type names: {ct_cols}")

    # Expression stats (exclude GENE and Average columns)
    expr_df = means_out[ct_cols].astype(float)
    print(f"\nExpression statistics (log2-CPM, low-expression filtered):")
    print(f"  Overall min: {expr_df.min().min():.4f}")
    print(f"  Overall max: {expr_df.max().max():.4f}")
    print(f"  Overall median: {expr_df.median().median():.4f}")
    print(f"\nPer cell type:")
    for ct in ct_cols:
        vals = expr_df[ct]
        nonzero = (vals > 0).sum()
        print(f"  {ct}: median={vals.median():.3f}, max={vals.max():.3f}, "
              f"nonzero={nonzero}/{len(vals)} ({100*nonzero/len(vals):.1f}%)")

    avg_col = means_out["Average"].astype(float)
    print(f"\nAverage column: median={avg_col.median():.3f}, max={avg_col.max():.3f}")
    print(f"\nNA values: {means_out[ct_cols].isna().sum().sum()}")
    print(f"Negative values: {(expr_df < 0).sum().sum()}")

    n_low_filtered = int((low_filter == 0).sum().sum())
    total_entries = low_filter.shape[0] * low_filter.shape[1]
    print(f"Low-expression entries zeroed (<1 CPM): {n_low_filtered}/{total_entries} "
          f"({100*n_low_filtered/total_entries:.1f}%)")

    print("=" * 60)
    print(f"\nOutput files:")
    print(f"  Means: {means_fp}")
    print(f"  Specificity: {spec_fp}")
    print(f"  Dimensions: {means_out.shape[0]} genes x {len(ct_cols)} cell types + Average")


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Convert h5ad to MAGMA cell type format for FUMA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Stage 1: QC and merge CellxGene h5ad files
  python3 h5ad_to_magma.py qc \\
    --input file1.h5ad file2.h5ad \\
    --output clean.h5ad \\
    --cell_type_col author_cell_type

  # Stage 2: Generate MAGMA-format output
  python3 h5ad_to_magma.py magma \\
    --input clean.h5ad \\
    --outdir output/ \\
    --cell_type_col author_cell_type \\
    --dataset_name Santo_Human_2025_Bladder_level2
        """)

    subparsers = parser.add_subparsers(dest="stage", help="Pipeline stage")
    subparsers.required = True

    # QC subcommand
    qc_parser = subparsers.add_parser("qc", help="Stage 1: QC filter, gene conversion, save clean h5ad")
    qc_parser.add_argument("--input", nargs="+", required=True, help="Input h5ad file(s)")
    qc_parser.add_argument("--output", required=True, help="Output clean h5ad file path")
    qc_parser.add_argument("--cell_type_col", default="author_cell_type",
                           help="obs column with cell type annotations (default: author_cell_type)")

    # MAGMA subcommand
    magma_parser = subparsers.add_parser("magma", help="Stage 2: Generate MAGMA-format TSV from clean h5ad")
    magma_parser.add_argument("--input", required=True, help="Input clean h5ad file (from QC stage)")
    magma_parser.add_argument("--outdir", required=True, help="Output directory")
    magma_parser.add_argument("--cell_type_col", default="author_cell_type",
                              help="obs column with cell type annotations (default: author_cell_type)")
    magma_parser.add_argument("--dataset_name", default="dataset",
                              help="Dataset name for QC report header")

    args = parser.parse_args()

    if args.stage == "qc":
        run_qc(args)
    elif args.stage == "magma":
        run_magma(args)


if __name__ == "__main__":
    main()
