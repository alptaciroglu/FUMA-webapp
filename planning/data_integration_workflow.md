# FUMA Dataset Integration Workflow
Generated from session: 2026-03-06

## Step-by-step integration plan

| Step | Who does it | Notes |
|---|---|---|
| 1. Identify exact datasets/accessions | You | Check paper supplementary tables or browse HuBMAP/CellxGene portal directly |
| 2. Download raw or processed count matrices | You | See Data Access section below for strategy |
| 3. Write preprocessing R/Python script (normalize -> log2 -> compute average -> MAGMA format) | Claude | Written based on existing FUMA magma_celltype.R patterns |
| 4. Register dataset in FUMA (blade file entry + dataset name) | Claude | Pure code change in celltype_options blade files |
| 5. Place processed .txt file on the server reference data path | You | Requires server access to FUMA production server |

## Required output file format (MAGMA --gene-covar input)

Tab-delimited .txt file placed at: `{ref_data_path}/celltype/{dataset_name}.txt`

| GENE | CellType_1 | CellType_2 | ... | Average |
|---|---|---|---|---|
| ENSG00000001 | 0.34 | 1.21 | ... | 0.87 |

- GENE column: Ensembl gene IDs, GRCh37/hg19, Ensembl v92
- One column per cell type (these become the VARIABLE names in MAGMA output)
- Average column: mean expression across all cell types (used as conditioning covariate)
- Values: average log2-normalised expression per gene per cell type, non-negative, shifted by minimum
- Mouse datasets require gene conversion to human homologs first (mm2hs)

## Dataset naming convention

Pattern: `{NextID}_{Author}_{Species}_{Year}_{Tissue}_{SubRegion}_level{N}`

Examples from existing FUMA datasets:
- `563_Xu_Human_2023_Kidney_CortexOfKidney_level1`
- `540_Travaglini_2020_Lung_level1`
- `443_Braun2023_Human_2023_FirstTrimester_Brain_CarnegieStage18_level1`

Level1 = broad cell types (e.g. epithelial, immune, stromal)
Level2 = fine-grained subtypes (e.g. proximal tubule S1, S2, S3)

Proposed HuBMAP pattern: `HuBMAP_Human_2026_{Tissue}_{SubRegion}_level{N}`

## UI registration (blade file)

Each dataset becomes one `<option>` entry in `resources/views/celltype/celltype_options/{tissue}_options.blade.php`:

```html
<option value="HuBMAP_Human_2026_Bladder_Bladder_level1"
data-section="Bladder/Human" data-key="0">
HuBMAP_Human_2026_Bladder_Bladder_level1</option>
```

The `data-section` attribute defines the tree path shown to users in the UI.
New tissues (Knee, Fallopian Tube) also need a new `@include` line added to `resources/views/pages/celltype.blade.php`.

---

## Data Access Strategy

### The core legal distinction
HuBMAP controlled access applies to raw genomic FASTQ/BAM files (individual-level, HIPAA risk).
FUMA only needs pseudo-bulk aggregated expression per cell type - not individual-level data.
The NIH GDS policy (which governs dbGaP/HuBMAP access) permits publishing derived/aggregated results
that cannot re-identify individuals. The aggregated profiles used by FUMA meet this standard.

### Three paths (in order of preference)

**Path 1 - CellxGene Census (try first, no access request needed)**
- URL: https://cellxgene.cziscience.com/
- Chan Zuckerberg Initiative portal with millions of publicly available scRNA-seq cells
- Has human bladder, ovary, uterus, and other tissues already processed with cell type annotations
- No data access request required
- Can be queried and downloaded programmatically via the `cellxgene-census` Python package
- If CellxGene has the tissue, this is the fastest path to full integration

**Path 2 - HuBMAP public processed data**
- HuBMAP releases processed h5ad/AnnData files separately from protected raw FASTQs
- The 401 "public" datasets in the metadata (assay_type = N/A) may be processed derivatives
- Check portal.hubmapconsortium.org for each tissue to see if h5ad files are publicly downloadable
- No access request needed if the processed file is already public

**Path 3 - HuBMAP dbGaP controlled access (fallback)**
- Apply via https://dbgap.ncbi.nlm.nih.gov/
- Requires: institutional sign-off, PI approval, research use statement
- FUMA as a VU Amsterdam academic tool qualifies for research use
- Timeline: typically 2-8 weeks for approval
- Publishing aggregated profiles in FUMA is permitted under NIH GDS policy

### Current decision
Check CellxGene first for Bladder data. If found, proceed with full integration immediately.

---

## First target dataset: Human Bladder (snRNA-seq)

Chosen because:
- Priority 1 in HuBMAP integration plan
- FUMA currently has ZERO human bladder scRNA-seq (only mouse: MouseCellAtlas, TabulaMuris)
- RNA-only (no ATAC complexity) = simplest pipeline
- High Q30 in HuBMAP metadata (94.6%)

FUMA UI location:
- Blade file: `resources/views/celltype/celltype_options/bladder_options.blade.php`
- data-section: `Bladder/Human` (new sub-section, currently file only has Mouse entries)

ATAC data (SNP2GENE): N/A for bladder (assay was snRNAseq-10xGenomics-v3, no paired ATAC)
