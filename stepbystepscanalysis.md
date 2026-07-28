what is the vars and obs in this: n_obs × n_vars = 3153 × 37764
That line is AnnData telling you the shape of your data matrix.

n_obs = 3153 → observations = cells (barcodes). Rows of the matrix.
n_vars = 37764 → variables = genes (features). Columns of the matrix.
So adata.X is a 3153 × 37764 matrix where entry [i, j] is the expression count of gene j in cell i.

The naming is generic on purpose — AnnData is a general-purpose data container, so it says "observations" and "variables" rather than "cells" and "genes". In single-cell work it's always cells × genes.

Two companion tables hang off those axes:

adata.obs — a pandas DataFrame with 3153 rows, one per cell: per-cell metadata like n_genes_by_counts, total_counts, pct_counts_mt, sample/condition labels.
adata.var — a DataFrame with 37764 rows, one per gene: gene metadata like gene_ids, mt (is it mitochondrial), n_cells_by_counts, and later highly_variable.
Quick way to see it yourself:


adata.shape          # (3153, 37764)
adata.obs.head()     # cell metadata
adata.var.head()     # gene metadata
adata.obs_names[:5]  # cell barcodes
adata.var_names[:5]  # gene names
One thing worth noticing: 37764 genes is the full reference annotation, so most of those columns are genes detected in zero or almost zero cells. That's exactly why the QC step you're on filters genes (e.g. sc.pp.filter_genes(adata, min_counts=1)) — n_vars will drop substantially after it, while n_obs drops when you filter low-quality cells.

why have they done gene annotation after this?
Bash List demonstration notebooks
IN
ls "e:/RNA-seq-courses/scRNA-seq_python_course_2024/Demonstrations/"

OUT
00_iris_setup_demo.ipynb
01_jiawei_cellranger_count.sh
01_jiawei_cellranger_demo.md
01_jiawei_cellranger_demo.pdf
01_jiawei_cellranger_mkref.sh
01_jiawei_cellranger_setup.md
01_jiawei_cellranger_setup.pdf
02_iris_qc_normalisation_demo.ipynb
02_iris_qc_normalisation_yuyao_mod_Feb2025_demo.ipynb
03_yuyao_feature_selection_dim_reduction_demo.ipynb
04_anna_clustering_cell_annotation_demo.ipynb
05_yuyao_data_integration_demo.ipynb
06_andrian_differential_expression_demo.ipynb
07_andrian_differential_abundance_demo.ipynb
08_nadav_01_b2c_mouse_brain_demo.ipynb
08_nadav_02_course_mouse_brain_demo.ipynb
D:RNA-seq-coursesscRNA-seq_python_course_2024Data
esults
D:RNA-seq-coursesscRNA-seq_python_course_2024Dataresults02
Data
Grep "n_obs|annotation|biomart|pyensembl|gene_ids|var\[|mygene" (in e:/RNA-seq-courses/scRNA-seq_python_course_2024/Demonstrations/02_iris_qc_normalisation_yuyao_mod_Feb2025_demo.ipynb)
36 lines of output
Read 02_iris_qc_normalisation_yuyao_mod_Feb2025_demo.ipynb
Read 02_iris_qc_normalisation_yuyao_mod_Feb2025_demo.ipynb (lines 281-310)
Because everything downstream depends on it — specifically, the mitochondrial QC metric can't be computed until the genes are annotated.

What CellRanger gives you

When sc.read_10x_mtx() reads features.tsv.gz, .var only has two columns:


                        gene_ids    feature_types
MIR1302-2HG      ENSG00000243485  Gene Expression
ENSG00000238009  ENSG00000238009  Gene Expression
No gene symbols for many genes (they fall back to the Ensembl ID), and crucially no chromosome. So at that point the object has no idea which of the 37764 genes are mitochondrial.

Why annotation has to come before QC

The very next thing the notebook does is:


sc.pp.calculate_qc_metrics(ETV6_RUNX1_1, qc_vars=["mt"], ...)
That qc_vars=["mt"] argument means "give me pct_counts_mt per cell", and it works by looking up a boolean column called mt in .var. That column only exists because of the annotation step. High mitochondrial fraction is one of the main signals of a dying/damaged cell, so it drives the barcode filtering in section 4.1 — no annotation, no pct_counts_mt, no mito-based filtering.

The intended chain is:

Read → you now know which gene IDs are present (you can't merge annotation before this — the merge is onto adata.var).
Annotate → query Biomart for gene_name + chromosome, merge onto .var.
Use the annotation → flag mt, drop genes on unassembled scaffolds.
QC metrics → needs mt.
Filter cells/genes → needs the QC metrics.
The scaffold removal also belongs before step 4, since dropping those genes changes each cell's total_counts, and the QC thresholds are computed from those totals.

One thing to notice in this Feb 2025 version

The chrom column never actually arrived — there's a note about a pybiomart bug, and the Biomart query only asks for ensembl_gene_id and external_gene_name, using chromosome as a filter rather than retrieving it as a column. Consequences:

The scaffold-filtering cell is a raw cell, not a code cell, so it doesn't run. That's why n_vars stays at 37764 the whole way through.
The mt column is instead built the long way round: a second Biomart query for MT genes only, then ETV6_RUNX1_1.var["mt"] = False followed by setting True where gene_ids is in that MT list.
The reusable sc_preprocess() function at the end of the notebook still uses adata.var["chrom"], so it assumes the annotation DataFrame you pass in does carry a chrom column — it's written for the un-bugged version.
Also worth flagging: sc.pp.filter_genes(ETV6_RUNX1_1, min_counts=0) in section 4.2 removes nothing (every gene has ≥ 0 counts), which is why n_vars prints 37764 straight after the text saying it's useful to remove the 11296 undetected genes. It needs min_counts=1.

e#1-1
this code gave feature types: ETV6_RUNX1_1.var why is it important
feature_types is the third column of CellRanger's features.tsv.gz. It tells you what kind of measurement each row of .var is — not every row in a 10x matrix is necessarily a gene.

In your data every value is Gene Expression, so right now it carries no information:


ETV6_RUNX1_1.var["feature_types"].value_counts()
# Gene Expression    37764
That's expected — this is a plain 3' single-cell gene expression run, one modality only.

When it matters

CellRanger writes multimodal data into a single matrix, stacking the different assays as extra rows in .var. The possible values include:

feature_types	What it is
Gene Expression	mRNA counts — the usual scRNA-seq data
Antibody Capture	CITE-seq surface protein tags, or hashtag oligos (HTOs)
CRISPR Guide Capture	Perturb-seq guide RNA assignments
Multiplexing Capture	10x CellPlex sample barcodes
Peaks	ATAC peaks in Multiome data
Custom	user-defined features
So if you ran CITE-seq with 30 antibodies, read_10x_mtx() would hand you an object with ~37800 vars where 30 of them are protein counts sitting in the same matrix as the genes. Everything in the QC and normalisation workflow would then be wrong if you ignored it:

total_counts per cell would mix mRNA counts and antibody counts, so your MAD-based filtering would be filtering on a meaningless quantity.
pct_counts_mt would be a fraction of that contaminated total.
sc.pp.normalize_total() assumes one modality with a shared size factor — antibody counts are orders of magnitude higher per feature and need their own normalisation (CLR, typically).
Highly variable gene selection and PCA in the next notebook would be dominated by the handful of protein features.
The habitual check

This is why it's worth glancing at the column right after reading data, even when you expect one modality:


adata.var["feature_types"].value_counts()
If you ever see more than one row in that output, split before doing anything else:


prot = adata[:, adata.var["feature_types"] == "Antibody Capture"].copy()
adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
and process each separately (or use MuData from the muon package, which is built for holding multiple modalities together properly).

There's also a convenience function that does the split for you when reading:


adata = sc.read_10x_mtx(path, gex_only=True)   # this is the default
gex_only=True is why many people never see the problem — it silently keeps only the Gene Expression rows. The notebook doesn't pass the argument, so it's relying on that default. Since your data is GEX-only anyway, it makes no difference here.

So: for this dataset, feature_types is a constant column you can ignore. Its importance is as a guard rail — it's the field that tells you whether the object you just loaded is what you assume it is.
