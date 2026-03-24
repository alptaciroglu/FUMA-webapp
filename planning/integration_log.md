# FUMA Update Log

Running log of all updates to FUMA, maintained by Alp and Yigit.
Covers three work areas: single-cell data additions, bug fixes / sanity checks, and the GRCh38 assembly update.
Detailed technical notes live here (reviewed by Tanya via PR); high-level status and to-dos are on Notion.

## Work areas

### 1. Addition of single-cell data
Two types of data to add:

**Reference scRNA-seq (standard):** Expression from healthy donors. Produces average
expression per cell type. FUMA uses this to test "which cell types are enriched for my
GWAS genes?" Use our existing script (`h5ad_to_magma.py`) which follows Tanya's pipeline.

**Case/control scRNA-seq:** Expression from disease patients AND healthy controls.
Enables testing which cell types show *disease-specific* expression changes for GWAS genes.
Doug and Emil's AD paper found results differ from reference-only data. Requires a different
MAGMA command (`--model` with different covariate structure) — a separate script is appropriate
here (per Tanya's advice). To be scoped later; follow up with Doug on methodology.

### 2. Bug fixes / sanity checks
Yigit noticed FUMA sometimes annotates exonic variants as intronic and vice versa.
FUMA also doesn't handle extremely low p-values below ~1e-300 well.
We will document and fix such issues as we encounter them, in coordination with Tanya.

### 3. GRCh38 update
Lower priority for now — defer until we have more experience with the codebase.
Detailed notes on what needs updating (eQTL, chromatin, MAGMA, etc.) are on Notion (GRCh38 Update page).

---

## Single-cell data integration

### Summary

| Date | Dataset | Tissue | Source | Status |
|---|---|---|---|---|
| 2026-03-09 | Santo et al. 2025 (GUDMAP) | Human Bladder | CellxGene | Processed, output verified (2026-03-24) |

### 2026-03-09: Human Bladder (GUDMAP/Santo et al. 2025)

### Source
- CellxGene collection: 0e54d4de-44f0-4d50-8649-b5c2bbe8f5d1
- Publication: Santo et al. (2025) iScience, DOI: 10.1016/j.isci.2024.111628
- Raw data: GSE267964

### Datasets downloaded (h5ad from CellxGene)
- e8799823-156b-4d74-b887-2305b86635eb.h5ad — Urothelial subset (14,776 cells)
- 4feb85d6-7a8e-4d7b-9b3b-5c7aa353a542.h5ad — Stromal subset (12,964 cells)
- cdf3c8cc-9eb3-4e5c-9d49-997552ffe55e.h5ad — Immune subset (2,094 cells)

### Plan
- Merge three subsets into one unified bladder dataset
- Process following Tanya's pipeline: QC -> CPM normalise -> log2 -> pseudo-bulk means
- Output MAGMA-format TSV + specificity TSV
- Test with FUMA Celltype command line tool

### Status
- [x] Identified dataset on CellxGene
- [x] Downloaded h5ad files to planning/data/GUDMAP/
- [x] Created README under planning/data/GUDMAP/
- [x] Confirmed no existing preprocessing script in FUMA repo
- [x] Reviewed Tanya's preprocessing pipeline (scrnaseq_viewer repo)
- [x] Cloned Tanya's pipeline to planning/scripts/tanya_pipeline/ (compute_sumstat_magma.py, helper functions, gene_names_human.txt)
- [x] Updated workflow docs with Tanya's instructions and priorities
- [x] Downloaded 876 example celltype files from fuma.ctglab.nl/downloadPage to downloads/celltype/
- [x] Verified example file format: GENE (Ensembl) | cell_type_1 | ... | Average, ~24750 genes, tab-delimited
- [x] Wrote per-dataset QC script following Tanya's pattern (qc_scrna_Santo_2025_Bladder.py)
- [x] Ran QC: 29,834 cells merged -> 22,284 after filtering (MT<10%), 24 cell types, 22,839 genes
- [x] Ran compute_sumstat_magma.py -> means + specificity TSVs
- [x] Compared output against FUMA examples: format, gene count (~22k), value ranges all match
- [x] Created placeholder FUMA-named file: `XXX_Santo_Human_2025_Bladder_Bladder_level1.txt` (ID TBD from Tanya)
- [x] Copied to shared output folder: `planning/data/fuma_celltype_output/`
- [x] Set up Snakemake workflow (Snakefile + config.json) for reproducible Stage 2 processing
- [ ] Test with FUMA Celltype command line (github.com/tanyaphung/FUMA_Celltype_cmd)
- [ ] Register dataset in blade file (bladder_options.blade.php)
- [ ] Create PR to vufuma/FUMA-webapp
- [ ] Verify on local FUMA server (setup target: before end of May)
- [ ] Send processed data to Tanya (only after local verification)

### Notes
- HuBMAP bladder data was tried first but is behind controlled access (dbGaP request pending via supervisor)
- CellxGene data used as demo to validate the integration workflow
- FUMA currently has zero human bladder datasets (only mouse: MouseCellAtlas, TabulaMuris)
- Tanya's preprint added 31 datasets; check supplementary tables for overlap before adding more
- Spatial transcriptomics (brain) and case/control data are next priorities after bladder demo
- For case/control data, check Doug/Emil AD paper for the separate MAGMA command used
- Follow up with Doug on the scatter plot methodology (gene selection decisions)
- 876 example celltype files downloaded (matches all blade-registered datasets). Format confirmed:
  GENE | cell_type columns | Average. Values: log2(CPM+1), tab-delimited, ~24750 genes per file.
- Tanya sharing preprocessing metadata docs (week of 2026-03-23) — metadata needed for publishing
- 2TB hard drive ordered ~2026-03-23 for Tanya's reference file transfer
- Switched from custom h5ad_to_magma.py to Tanya's actual pipeline code (consistency > convenience)
- CellxGene gotcha: adata.X is pre-normalised; must use adata.raw.X for raw counts
- CellxGene gotcha: var index has Ensembl IDs already, so symbol conversion is skipped.
  But must still filter genes against gene_names_human.txt to match the implicit filtering
  in Tanya's conversion step (otherwise ~31k genes including lncRNAs vs expected ~20-25k).
- anndata.concat drops var columns when they differ across files; save var metadata before merge
- anndata.read() deprecated in newer versions; use anndata.read_h5ad()

---

## Bug fixes / sanity checks

No entries yet. Issues will be logged here as they are encountered and investigated.

---

## GRCh38 update

Deferred. See Notion (Yigit's GRCh38 Update page) for detailed component-by-component notes.
