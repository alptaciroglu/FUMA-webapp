# FUMA Dataset Integration Workflow

Generated from sessions: 2026-03-06, updated 2026-03-09

## Key references

- **Tanya's preprocessing code**: github.com/tanyaphung/scrnaseq_viewer/tree/main/code/preprocessing
- **Tanya's MAGMA output code**: github.com/tanyaphung/scrnaseq_viewer/code/postprocessing/compute_sumstat_magma.py
- **FUMA Celltype command line**: github.com/tanyaphung/FUMA_Celltype_cmd (for testing locally)
- **Tanya's biorxiv preprint** (31 datasets added): biorxiv.org/content/10.64898/2025.12.05.692533v1.abstract
- **FUMA webapp PR workflow**: github.com/vufuma/FUMA-webapp (Tanya reviews and deploys)

## Step-by-step integration plan

For each new dataset, follow these steps in order. Do not skip any.

| Step | Who | What | Details |
|---|---|---|---|
| 1. Identify dataset | Alp | Find dataset + accession | Check CellxGene, UCSC Cell Browser, DISCO v1.0, HuBMAP. Cross-ref existing 876 FUMA datasets + Tanya's 31 new datasets to avoid duplicates. |
| 2. Download data | Alp | Get h5ad / count matrices | See Data Access section below. |
| 3. Write QC script | Claude | Per-dataset QC script | Follow Tanya's pattern (see `planning/scripts/tanya_pipeline/`). One script per dataset. |
| 4. Run QC (Stage 1) | Alp/Claude | QC filter, gene conversion | Output: clean h5ad with `ensembl` column + log.txt + table2.csv + table3.csv. |
| 5. Run MAGMA (Stage 2) | Alp/Claude | `compute_sumstat_magma.py` | Output: `means_cell_log_counts_pM.tsv` + `spec_cell_log_counts_pM.tsv`. |
| 6. Verify output | Alp/Claude | Compare against FUMA examples | Check: gene count (~20-25k), value range ([0, ~16]), GENE col, Average col, no NAs/negatives. |
| **7. Rename output** | **Alp/Claude** | **Rename means TSV to FUMA convention** | **`{ID}_{Author}_{Species}_{Year}_{Tissue}_{SubRegion}_level{N}.txt` — see naming section below. This is easy to forget!** |
| 8. Test locally | Alp | Test with FUMA Celltype CLI | Use github.com/tanyaphung/FUMA_Celltype_cmd (after local FUMA setup). |
| 9. Register in blade | Claude | Add `<option>` to blade file | `value` must exactly match the `.txt` filename from step 7. |
| 10. Create PR | Alp/Claude | PR to vufuma/FUMA-webapp | Code only (blade file). No data files in git. |
| 11. Send data to Tanya | Alp | Slack or email | Only after verifying on local FUMA. Include: renamed .txt file + metadata (await Tanya's docs). |
| 12. Tanya deploys | Tanya | Places data on server | File goes to `/data/MAGMA/celltype/{dataset_name}.txt`. |

### Gotchas learned during bladder dataset (2026-03-24)

These are easy to miss. Check each one when processing a new dataset:

- **CellxGene raw counts**: `adata.X` is pre-normalised. Use `adata.raw.X.copy()` for raw integer counts.
- **Gene filtering for CellxGene data**: Ensembl IDs are already in the var index, so symbol conversion
  is skipped. But you MUST still filter genes against `gene_names_human.txt` to match the implicit
  filtering in Tanya's conversion step. Without this you get ~31k genes (including lncRNAs) instead of ~20-25k.
- **anndata.concat drops var columns**: Save var metadata from the first file before merging.
- **anndata.read() is deprecated**: Use `anndata.read_h5ad()`.
- **compute_sumstat_magma.py output dir**: Script writes to `{outdir}/{ct_colname}/` — create this directory first.
- **Renaming**: `means_cell_log_counts_pM.tsv` is NOT the final filename. It must be renamed
  to the FUMA convention (step 7) before sending to Tanya or registering in the blade file.

## Data priorities (from Tanya, 2026-03-09)

1. **Spatial transcriptomics data** — priority on brain tissues. Find publicly available datasets.
2. **Additional scRNA-seq data** — many tissues still missing. Sources to check:
   - CellxGene: cellxgene.cziscience.com/collections
   - UCSC Cell Browser: cells.ucsc.edu
   - DISCO v1.0: disco.cs.ucdavis.edu
3. **Case/control RNA-seq data** — scRNA-seq from both disease patients and healthy controls,
   not just healthy reference data. Allows testing which cell types show *disease-specific*
   expression changes for GWAS genes (vs. reference data which only tests general enrichment).
   Doug and Emil's AD paper found results differ when using case/control vs. reference-only data.
   Requires a different MAGMA command (`--model` with different covariate structure), so a
   separate script is appropriate here (per Tanya). Follow up with Doug on scatter plot
   methodology and gene selection decisions.

## Required output file format (MAGMA --gene-covar input)

Tab-delimited .tsv file placed at: `{ref_data_path}/celltype/{dataset_name}.txt`

| GENE | CellType_1 | CellType_2 | ... | Average |
|---|---|---|---|---|
| ENSG00000001 | 0.34 | 1.21 | ... | 0.87 |

- GENE column: Ensembl gene IDs, GRCh37/hg19, Ensembl v92
- One column per cell type (these become the VARIABLE names in MAGMA output)
- Average column: mean expression across all cell types (used as conditioning covariate)
- Values: mean log2(CPM+1) per gene per cell type, with low-expression filter (<1 CPM zeroed)
- Mouse datasets require gene conversion to human homologs first (mm2hs)
- Two output files: means TSV (for MAGMA) + specificity TSV

## Dataset naming convention

Pattern: `{ID}_{Author}_{Species}_{Year}_{Tissue}_{SubRegion}_level{N}`

- **ID**: Numeric identifier. Existing FUMA datasets go up to 595. Some older datasets have no numeric
  prefix (e.g. `MouseCellAtlas_Bladder`, `TabulaMuris_FACS_Bladder`). ID assignment likely handled by
  Tanya — confirm with her metadata documentation (expected week of 2026-03-23). Use `XXX` as placeholder
  until ID is assigned.
- **Author**: First author surname (e.g. `Santo`, `Xu`, `Siletti`)
- **Species**: `Human` or `Mouse`
- **Year**: Publication year
- **Tissue / SubRegion**: Anatomical location (e.g. `Bladder_Bladder`, `Kidney_CortexOfKidney`)
- **level{N}**: Cell type hierarchy depth
  - Level1 = broad cell types (e.g. epithelial, immune, stromal)
  - Level2 = fine-grained subtypes (e.g. proximal tubule S1, S2, S3)
  - Level3 = cluster-level (if available)

Examples from existing FUMA datasets:
- `563_Xu_Human_2023_Kidney_CortexOfKidney_level1`
- `540_Travaglini_Human_2020_Lung_level1`
- `443_Braun2023_Human_2023_FirstTrimester_Brain_CarnegieStage18_level1`
- `3_Gabitto_MTG_Human_2023_level1` (no SubRegion, MTG is the tissue)

**How to rename**: Copy `means_cell_log_counts_pM.tsv` to `{name}.txt`. The `.txt` extension is
required (not `.tsv`). The specificity file (`spec_cell_log_counts_pM.tsv`) is kept as-is for reference
but is not deployed to the FUMA server (only the means file is used by MAGMA).

## UI registration (blade file)

Each dataset becomes one `<option>` entry in `resources/views/celltype/celltype_options/{tissue}_options.blade.php`:

```html
<option value="Santo_Human_2025_Bladder_Bladder_level2"
data-section="Bladder/Human" data-key="0">
Santo_Human_2025_Bladder_Bladder_level2</option>
```

The `data-section` attribute defines the tree path shown to users in the UI.
New tissues also need a new `@include` line added to `resources/views/pages/celltype.blade.php`.

## Verified: what needs to change when adding a dataset

**Code change (goes in PR):** Only the blade file. The rest of the system is dynamic:

- `CellController.php` — receives dataset names from the form, passes them through. No hardcoded dataset list.
- `magma_celltype.R` — reads dataset names from job parameters (`params.config`), not from any static config.
- `celltype.js` — UI logic only, no dataset-specific references.
- Database — the old `celltype` table was dropped (migration `2023_06_15_133548`). No dataset metadata stored in DB.

**Data file (sent to Tanya):** The `.txt` file at `/data/MAGMA/celltype/{dataset_name}.txt` on the server (mounted into Docker as `/data`).

**Critical link:** The `value` attribute in the blade `<option>` must exactly match the `.txt` filename on the server. This is the only connection between the code and the data.

*Verified 2026-03-10 by reading CellController.php, magma_celltype.R, celltype.js, and database migrations.*

## PR and deployment workflow

**What goes in the PR (code only — no data files):**
- Blade file changes: `resources/views/celltype/celltype_options/{tissue}_options.blade.php`
  (add `<option>` entries for new datasets)
- If adding a new tissue: add `@include` line in `resources/views/pages/celltype.blade.php`
- Any controller or script changes (e.g. `app/Jobs/CelltypeProcess.php`, `scripts/magma_celltype/`)

**What does NOT go in the PR (data files — sent separately to Tanya):**
- The processed TSV files (means + specificity) produced by our preprocessing script
- These live on the server at `{ref_data_path}/celltype/{dataset_name}.txt`
  (mounted into Docker containers as `/data/celltype/`)
- Data files are not tracked in the git repo — the repo only has code
- Send data files to Tanya directly; she places them on the server

**Steps:**
1. Process data locally, test with FUMA Celltype command line tool
2. Create PR to github.com/vufuma/FUMA-webapp with code changes only
3. Send processed TSV files to Tanya separately (email or shared drive)
4. Tanya reviews PR, merges, and places data files on the server
5. For local development: set up a local FUMA server (meeting with Tanya scheduled 20 May)

---

## Data Access Strategy

### Databases to check (in order of preference)

**1. CellxGene (try first, no access needed)**
- URL: https://cellxgene.cziscience.com/
- Fully public, downloadable h5ad with cell type annotations
- Programmatic access via `cellxgene-census` Python package

**2. UCSC Cell Browser**
- URL: https://cells.ucsc.edu/
- Public scRNA-seq datasets with cell type annotations

**3. DISCO v1.0**
- URL: https://disco.cs.ucdavis.edu/
- Integrated single-cell omics database

**4. HuBMAP public processed data**
- portal.hubmapconsortium.org
- Some processed h5ad files may be publicly downloadable
- Note: As of 2026-03-09, bladder Salmon-processed files are still behind controlled access

**5. HuBMAP dbGaP controlled access (fallback)**
- Apply via https://dbgap.ncbi.nlm.nih.gov/
- Requires: institutional sign-off, PI approval, research use statement
- Timeline: typically 2-8 weeks for approval

### Legal note
FUMA only needs pseudo-bulk aggregated expression per cell type — not individual-level data.
The NIH GDS policy permits publishing derived/aggregated results that cannot re-identify individuals.

---

## Preprocessing pipeline

### Reference implementation

Tanya Phung's established pipeline (github.com/tanyaphung/scrnaseq_viewer):

**Stage 1 — QC & gene conversion** (per-dataset scripts in `code/preprocessing/`):
1. Load h5ad (raw integer counts)
2. QC filter: min 200 genes/cell, min 3 cells/gene, MT% < 10%
3. Convert gene symbols to Ensembl IDs via `gene_names_human.txt` (genenames.org)
4. For mouse data: convert via MGI IDs to human homologs
5. Save clean h5ad + log file + cell count tables

**Stage 2 — MAGMA format generation** (`code/postprocessing/compute_sumstat_magma.py`):
1. Filter to genes with valid Ensembl IDs
2. CPM normalise (`sc.pp.normalize_total(target_sum=1e6)`)
3. Compute mean CPM per cell type (pre-log)
4. Low-expression filter: zero out entries with < 1 CPM per cell type
5. Log2(x+1) transform (`sc.pp.log1p(base=2)`)
6. Compute mean log2-CPM per cell type (post-log)
7. Apply low-expression filter to log means
8. Compute specificity values (proportion of total expression per gene)
9. Add Average row, sanitise cell type names, transpose, output two TSV files:
   - `means_cell_log_counts_pM.tsv` — expression levels (this is the MAGMA input)
   - `spec_cell_log_counts_pM.tsv` — specificity values

### Our local copy

All of Tanya's pipeline code is cloned to `planning/scripts/tanya_pipeline/`:

| File | Purpose |
|---|---|
| `compute_sumstat_magma.py` | Stage 2: CPM normalise -> log2 -> means per cell type -> TSV (Tanya's code, only change: `anndata.read` -> `anndata.read_h5ad`) |
| `compute_sumstat_magma.smk` | Snakemake workflow for batch processing (reference only) |
| `single_cell_helper_functions_v3.py` | Gene symbol -> Ensembl conversion + mouse -> human homolog mapping |
| `gene_names_human.txt` | HGNC gene name reference (from Tanya's repo, originally genenames.org) |
| `fumacelltype_datasets_master.csv` | Master list of all 388 FUMA celltype datasets with IDs, levels, cell counts |
| `scrna_fumacelltype_master.csv` | Detailed metadata per dataset (tissue, author, year, cell count, PMID) |
| `scRNAseq_categorized_by_brain_regions.csv` | Brain region categorisation for brain datasets |
| `scrna_env.yml` | Conda environment spec (reference for dependency versions) |
| `Snakefile` | **Our Snakemake workflow**: runs Stage 2 for all levels, renames output to FUMA convention (adapted from `compute_sumstat_magma.smk`) |
| `config.json` | Dataset definitions for Snakefile (id, h5ad path, fuma_name, levels) |
| `qc_scrna_*.py` | Per-dataset QC scripts (Stage 1) — one per dataset, following Tanya's pattern |

### Usage

```bash
# Stage 1: Run per-dataset QC script (writes clean h5ad + log + tables)
python3 planning/scripts/tanya_pipeline/qc_scrna_Santo_2025_Bladder.py \
  planning/data/GUDMAP Santo_Human_2025_Bladder

# Stage 2 + rename: Run Snakemake (handles all datasets x levels from config.json)
cd planning/scripts/tanya_pipeline
snakemake --configfile config.json -c1
```

The Snakefile reads `config.json` which lists every dataset with its h5ad path, FUMA name,
and cell type levels. It runs `compute_sumstat_magma.py` for each dataset x level, then
copies the output to `planning/data/fuma_celltype_output/` with FUMA-convention filenames.

For a dataset with 2 levels, set `"levels": [1, 2]` in `config.json`. The QC script must
create both `cell_type_level_1` and `cell_type_level_2` columns in obs.

### Shared output folder

All finalised FUMA celltype files are collected in `planning/data/fuma_celltype_output/`.
This is the folder to send to Tanya — one folder with all processed datasets.
The Snakefile writes directly to this folder (configured via `outdir` in `config.json`).

### Cell type level conventions (from fumacelltype_datasets_master.csv)

- 164 datasets have 1 level, 185 have 2, 39 have 3
- Level count depends on what the original study authors annotated
- Datasets with 3 levels only publish levels 1+2 to FUMA (never level 3)
- Do not construct levels that the original authors didn't define

### QC of output
- Compare gene count against expected ~20-25k (existing FUMA datasets range from ~15k to ~25k)
- Value range should be [0, ~16] (log2 CPM+1 scale)
- Verify: GENE column (Ensembl IDs), cell type columns, Average column, tab-delimited
- No NAs, no negative values, zeros expected for low-expression filtered entries
- Compare against example files in `downloads/celltype/`
- Test with FUMA Celltype command line tool (github.com/tanyaphung/FUMA_Celltype_cmd)

---

## First target dataset: Human Bladder (GUDMAP/Santo et al. 2025)

Chosen because:
- FUMA currently has ZERO human bladder scRNA-seq (only mouse: MouseCellAtlas, TabulaMuris)
- RNA-only (no ATAC complexity) = simplest pipeline
- Publicly available on CellxGene

Source: CellxGene collection 0e54d4de-44f0-4d50-8649-b5c2bbe8f5d1
Publication: Santo et al. (2025) iScience, DOI: 10.1016/j.isci.2024.111628

FUMA UI location:
- Blade file: `resources/views/celltype/celltype_options/bladder_options.blade.php`
- data-section: `Bladder/Human` (new sub-section, currently file only has Mouse entries)

---

## Timeline and coordination

- Tanya is working on FUMA v2.0, unavailable until ~end of May 2026
- FUMA v2.0 goes live ~end of May; Tanya available to help with local install after that
- **Meeting scheduled: 20 May 2026, 10am** — Tanya walks through local FUMA installation
- Until then: process data, download example files to verify format, document everything
- Tanya will share preprocessing metadata documentation (week of 2026-03-23)
- **2TB external hard drive ordered** (~March 23, 1-2 week delivery) — Tanya will copy reference files to it
- **Local FUMA server setup target: before end of May** — Tanya shared Windows install guide (`How_to_install_FUMA_locally.pdf`), needs macOS adaptation. Document the macOS process (Tanya's request).
- **Only send processed data to Tanya after verifying on local FUMA server**
- When ready to share: message Tanya on Slack or email

## Resources from Tanya (updated 2026-03-23)

- Preprocessing code: github.com/tanyaphung/scrnaseq_viewer/tree/main/code/preprocessing
- FUMA contributing docs: fuma-docs.readthedocs.io/en/latest/contributing.html
- Example single-cell data files: fuma.ctglab.nl/downloadPage (scroll to bottom)
- Local FUMA install guide (Windows): `How_to_install_FUMA_locally.pdf` (shared on Slack)
