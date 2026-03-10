# FUMA Dataset Integration Workflow

Generated from sessions: 2026-03-06, updated 2026-03-09

## Key references

- **Tanya's preprocessing code**: github.com/tanyaphung/scrnaseq_viewer/tree/main/code/preprocessing
- **Tanya's MAGMA output code**: github.com/tanyaphung/scrnaseq_viewer/code/postprocessing/compute_sumstat_magma.py
- **FUMA Celltype command line**: github.com/tanyaphung/FUMA_Celltype_cmd (for testing locally)
- **Tanya's biorxiv preprint** (31 datasets added): biorxiv.org/content/10.64898/2025.12.05.692533v1.abstract
- **FUMA webapp PR workflow**: github.com/vufuma/FUMA-webapp (Tanya reviews and deploys)

## Step-by-step integration plan

| Step | Who does it | Notes |
|---|---|---|
| 1. Identify exact datasets/accessions | Alp | Check databases: CellxGene, UCSC Cell Browser, DISCO v1.0, HuBMAP |
| 2. Download raw or processed count matrices | Alp | See Data Access section below for strategy |
| 3. Preprocess: QC filter, gene conversion | Alp/Claude | Stage 1 of h5ad_to_magma.py or Tanya's per-dataset QC scripts |
| 4. Generate MAGMA-format output | Alp/Claude | Stage 2 of h5ad_to_magma.py (CPM normalise -> log2 -> pseudo-bulk) |
| 5. Test with FUMA Celltype on command line | Alp | Use github.com/tanyaphung/FUMA_Celltype_cmd |
| 6. Register dataset in FUMA (blade file entry) | Claude | Pure code change in celltype_options blade files |
| 7. Create PR to vufuma/FUMA-webapp | Alp/Claude | Tanya reviews and makes live |
| 8. Place processed data on server reference path | Tanya | Tanya adds to reference database (await her detailed instructions) |

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

Pattern: `{NextID}_{Author}_{Species}_{Year}_{Tissue}_{SubRegion}_level{N}`

Examples from existing FUMA datasets:
- `563_Xu_Human_2023_Kidney_CortexOfKidney_level1`
- `540_Travaglini_2020_Lung_level1`
- `443_Braun2023_Human_2023_FirstTrimester_Brain_CarnegieStage18_level1`

Level1 = broad cell types (e.g. epithelial, immune, stromal)
Level2 = fine-grained subtypes (e.g. proximal tubule S1, S2, S3)

## UI registration (blade file)

Each dataset becomes one `<option>` entry in `resources/views/celltype/celltype_options/{tissue}_options.blade.php`:

```html
<option value="Santo_Human_2025_Bladder_Bladder_level2"
data-section="Bladder/Human" data-key="0">
Santo_Human_2025_Bladder_Bladder_level2</option>
```

The `data-section` attribute defines the tree path shown to users in the UI.
New tissues also need a new `@include` line added to `resources/views/pages/celltype.blade.php`.

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

### Our script

`planning/scripts/h5ad_to_magma.py` — follows Tanya's pipeline, adapted for CellxGene h5ad
files which already have Ensembl IDs in the var index (skipping symbol-to-Ensembl conversion).

For non-CellxGene data or data without Ensembl IDs, use Tanya's helper functions:
- `single_cell_helper_functions_v3.py` → `add_gene_names_human()` / `add_gene_names_mouse()`
- Requires `gene_names_human.txt` from genenames.org

### Usage

```bash
# Stage 1: QC and merge CellxGene h5ad files
python3 planning/scripts/h5ad_to_magma.py qc \
  --input planning/data/GUDMAP/*.h5ad \
  --output planning/data/GUDMAP/clean_merged.h5ad \
  --cell_type_col author_cell_type

# Stage 2: Generate MAGMA-format output
python3 planning/scripts/h5ad_to_magma.py magma \
  --input planning/data/GUDMAP/clean_merged.h5ad \
  --outdir planning/data/GUDMAP/output/ \
  --cell_type_col author_cell_type \
  --dataset_name Santo_Human_2025_Bladder_level2
```

### QC of output
- Compare gene count against expected ~20k protein-coding genes
- Check expression distribution per cell type (min/median/max reported by script)
- Verify no NA or negative values (reported by script)
- Check low-expression filter rate (reported by script)
- Compare structure against existing FUMA celltype files on server
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

- Tanya is working on a big FUMA update, unavailable until early May 2026
- FUMA update goes live early May; Tanya available for meetings after that
- **Meeting scheduled: 20 May 2026, 10am** — Tanya walks through local FUMA installation
- Until then: process data, test with FUMA Celltype command line, document everything
- Tanya will provide detailed instructions on adding processed data to the server
