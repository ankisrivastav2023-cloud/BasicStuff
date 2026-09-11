The big picture first
Before any code: what is the object you're manipulating?

A droplet scRNA-seq experiment (10x Genomics) works like this. You put a cell suspension into a microfluidic chip. Each cell gets captured in an oil droplet along with a gel bead. That bead carries millions of copies of a DNA barcode — and every bead has a different barcode. Inside the droplet, the cell lyses, its mRNA binds to the bead's barcodes, and reverse transcription tags every mRNA molecule with two things:

a cell barcode (which droplet it came from)
a UMI — Unique Molecular Identifier (which individual mRNA molecule it was)
Then everything is pooled and sequenced together. CellRanger (the step before this notebook) demultiplexes the reads and counts, for each cell barcode and each gene, how many distinct UMIs were seen. That gives you the counts matrix: rows = barcodes (~ cells), columns = genes, values = number of mRNA molecules of that gene detected in that cell.

That matrix is what this notebook receives. Everything in the notebook is about answering one question: which of these barcodes are actually real, single, healthy cells, and how do I make their counts comparable to each other?

Two things you should be able to say out loud in an interview:

The matrix is ~90–95% zeros. This is "sparsity" — partly biological (the gene really isn't expressed) and partly technical dropout (only ~10–30% of a cell's mRNA is actually captured, so lowly-expressed transcripts are often missed entirely). You cannot tell these two apart from the data alone.
A "cell barcode" is not a cell. It's a droplet. It might contain a real cell, no cell (just floating RNA), a dying cell, or two cells. QC is the process of figuring out which is which.
Cell 7 — the imports

import os
from glob import glob
import numpy as np
from scipy.stats import median_abs_deviation
from scipy.sparse import csr_matrix
import pandas as pd
import scanpy as sc
import anndata as ad
import pybiomart as bm
import seaborn as sns
sns.set_theme()
Since you're new to Python, plainly:

Import	What it is	Why it's here
os	talks to the file system	making output folders
numpy (np)	fast numerical arrays	the maths on the count matrix
median_abs_deviation	one statistic from scipy	the outlier rule in §4.1
csr_matrix	"compressed sparse row" format	stores a 95%-zero matrix by only recording the non-zeros — huge memory saving
pandas (pd)	spreadsheet-like tables (DataFrame)	the cell and gene metadata tables
scanpy (sc)	the Python single-cell toolkit	does almost every analysis step
anndata (ad)	the data container scanpy uses	holds matrix + metadata together
pybiomart (bm)	queries Ensembl's gene database	converting gene IDs to gene names
seaborn (sns)	plotting	the QC histograms
sc.pp. = preprocessing, sc.pl. = plotting, sc.tl. = tools (analysis). Learning that prefix convention makes scanpy much less intimidating.

The dataset (cell 5)
Caron et al. 2020 — childhood acute lymphoblastic leukemia (ALL), a cancer of immature lymphocytes. Four groups:

ETV6-RUNX1 (4 replicates) — B-cell precursor ALL with a t(12;21) translocation fusing the ETV6 and RUNX1 genes
HHD (2) — "high hyperdiploid" B-ALL, extra whole chromosomes
PRE-T (2) — T-cell precursor ALL
PBMMC (3) — healthy paediatric bone marrow mononuclear cells, the controls
The notebook develops the pipeline on one sample (SRR9264343, an ETV6-RUNX1 replicate) and then applies it to all of them. That's a good working habit and worth saying in an interview: develop interactively on one sample, then wrap it in a function and loop.

§1 (cells 8–16) — Ambient RNA correction
The problem. When you prepare a cell suspension, some cells break. Their mRNA spills into the buffer. Every droplet you then make is filled with a bit of that soup. So every cell barcode picks up counts for genes it never expressed — contamination that looks like real expression.

This matters enormously in this dataset: if the leukemic blasts are fragile and lyse, haemoglobin or blast-specific transcripts smear across all cells, and you'd conclude cell types express genes they don't.

The barcode-rank plot (the image in cell 9). Sort all barcodes by total UMI count, highest first, and plot on log-log axes. You see:

a plateau on the left — real cells, high counts
a sharp knee/cliff
a long tail on the right — droplets with only a few hundred counts
That tail is empty droplets. But they're not empty of RNA — they contain pure ambient soup. That's the key insight: empty droplets give you a free profile of the contamination, because whatever's in them is by definition not from a cell.

CellBender learns that ambient profile from the empty droplets and subtracts the expected contamination from every real cell, using a deep generative model. Crucially it needs the raw matrix (all ~700,000 barcodes including empties), not the filtered one — the empties are the training signal.


ETV6_RUNX1_1 = sc.read_10x_mtx(f"{prefix_inputs}/SRR9264343/outs/raw_feature_bc_matrix/")
ETV6_RUNX1_1.write(f"{prefix_outputs}/cellranger_h5ad/SRR9264343.raw.h5")
Then in bash:


singularity exec cellbender.sif cellbender remove-background --cuda \
    --input ".../SRR9264343.raw.h5ad" --output ".../SRR9264343.denoised.h5"
singularity is a container system (like Docker but usable on shared HPC clusters where you don't have root). --cuda uses a GPU. Without a GPU this takes ~5 hours per sample, which is exactly why the notebook skips this step and the cells are marked raw rather than code — raw cells don't execute, so you can't run them by accident.

Alternatives to know: SoupX (R, simpler, uses marker genes), DecontX, CellBender (best-performing but slowest).

§2 (cells 17–21) — Reading the data

ETV6_RUNX1_1 = sc.read_10x_mtx(f"{prefix_inputs}/SRR9264343/outs/filtered_feature_bc_matrix/")
filtered vs raw: CellRanger writes both. raw_feature_bc_matrix = every one of the ~737,000 possible barcodes. filtered_feature_bc_matrix = only barcodes CellRanger's EmptyDrops-style algorithm called as real cells (~thousands). Since CellBender was skipped, the notebook starts from CellRanger's filtered output.

The folder contains three files, which read_10x_mtx stitches together:

matrix.mtx.gz — the counts in Matrix Market sparse format: a list of (row, column, value) triplets, only for non-zeros
barcodes.tsv.gz — the row names (cell barcodes, e.g. AAACCTGAGACAATAC-1)
features.tsv.gz — the column names (gene IDs)
The AnnData object
This is the thing to understand. AnnData = "Annotated Data". One Python object holding:


        genes (var) →
      ┌──────────────────┐
cells │                  │   .X      the count matrix
(obs) │       .X         │   .obs    a DataFrame, one row per CELL
  ↓   │                  │   .var    a DataFrame, one row per GENE
      └──────────────────┘   .layers alternative versions of .X (same shape)
                             .obsm   per-cell matrices (PCA, UMAP coords)
.X — the matrix. Rows = cells, columns = genes. (Note: this is the transpose of the convention in R/Seurat and bulk RNA-seq, where genes are rows. A very common confusion — worth flagging if asked.)
.obs — observations = cells. Every QC metric you compute becomes a new column here.
.var — variables = genes. Gene names, chromosome, mitochondrial flag.
.layers — dictionary of same-shaped matrices; used later to hold raw counts and two normalisations side by side.
The genius of it: when you subset (adata[keep_cells, :]), .X, .obs and .var all stay in sync automatically. You can never accidentally misalign your metadata with your matrix.

Printing the object (cell 20) shows n_obs × n_vars — cells × genes.

§3 (cells 22–47) — Gene annotation
ETV6_RUNX1_1.var shows the genes are labelled with Ensembl IDs like ENSG00000139618. These are stable, unambiguous machine identifiers. But you need two extra things:

Gene symbols (BRCA2) so results are human-readable
Which chromosome each gene is on — specifically to find the mitochondrial (MT) genes, which are the workhorse QC metric

h38_genes = pd.read_csv(".../hg38_genes_dataframe.csv", index_col=0)
h38_genes = h38_genes.rename(columns={"Gene stable ID": "gene_ids", "Gene name": "gene_name"})
Why the rename matters: to merge two tables, the shared column must have the same name in both. adata.var calls it gene_ids; Ensembl calls it Gene stable ID. Renaming makes the join possible.

Note that cells 25 and 41 are commented out — the live pybiomart queries to Ensembl. Cell 40 explains why: a bug in the pybiomart–Ensembl connection meant it had to be run twice and couldn't return chromosomes cleanly. So someone ran it once, saved to CSV, and the notebook now reads the CSV. This is good practice regardless — pinning annotation to a saved file makes the analysis reproducible instead of depending on whatever Ensembl serves today. If an interviewer asks about reproducibility, this is a clean example.

The merge (cells 33–36)

gene_annot = (
  ETV6_RUNX1_1.var
  .merge(h38_genes, how="left", on="gene_ids")
  .set_axis(ETV6_RUNX1_1.var.index)
)
ETV6_RUNX1_1.var = gene_annot.loc[ETV6_RUNX1_1.var_names]
Step by step:

.merge(..., how="left", on="gene_ids") — a left join: keep every row of adata.var, and pull in matching annotation where it exists. Genes with no match get NaN rather than being dropped. Critical: an inner join would silently delete genes and break alignment with .X.
.set_axis(...) — pandas merge throws away the row index and replaces it with 0,1,2,3… This line puts the original gene names back as the index.
.loc[ETV6_RUNX1_1.var_names] — reorders to exactly match the column order of .X.
This is the single most dangerous point in the notebook, and a fantastic interview answer. If a merge duplicates or reorders rows, your gene metadata no longer lines up with your count matrix — every downstream result is silently wrong, with no error message. Hence three separate defensive moves: left join, restore index, explicit reorder. (Notice the workflow function in cell 103 adds a fourth: gene_annot.drop_duplicates(subset="gene_ids"), guaranteeing the merge can't duplicate rows.)

Chromosome filtering (cell 38)

vars_to_keep = ETV6_RUNX1_1.var["chrom"].isin([str(i) for i in range(1, 23)] + ["X", "Y", "MT"])
ETV6_RUNX1_1 = ETV6_RUNX1_1[:, vars_to_keep].copy()
[str(i) for i in range(1,23)] is a list comprehension — builds ["1","2",...,"22"]. Plus X, Y, MT.
.isin(...) gives True/False per gene.
adata[:, mask] — all rows, selected columns. In [rows, columns] notation, : means "everything".
This drops genes on unplaced scaffolds — bits of sequence in the genome build that couldn't be assigned to a chromosome. They're poorly annotated and mostly noise.

.copy() — cell 39 flags this explicitly and it's a classic gotcha. Subsetting an AnnData returns a View — a lightweight window onto the original object, not new data. Views are memory-efficient but you get warnings (or wrong behaviour) when you try to modify them. .copy() forces a real, independent object. Rule of thumb: whenever you subset and intend to keep working with the result, add .copy().

Flagging mitochondrial genes (cells 45–47)

ETV6_RUNX1_1.var["mt"] = False
ETV6_RUNX1_1.var.loc[ETV6_RUNX1_1.var.gene_ids.isin(h38_MT['Gene stable ID']), "mt"] = True
Create a column that's all False, then set it True for the 13 protein-coding mitochondrial genes. (Cell 45, the clean one-liner var["mt"] = var["chrom"] == "MT", is a raw cell — it's the version that works when the chromosome column is available.)

Why mitochondria? This is the most important biological idea in the QC section. Mitochondrial transcripts are inside the mitochondrion, not the cytoplasm. When a cell is dying or its membrane is damaged, cytoplasmic mRNA leaks out of the cell but mitochondrial mRNA is retained inside the organelles. So the proportion of counts that are mitochondrial rises sharply. High %MT ⇒ stressed, dying, or broken cell.

Caveats you should mention: the threshold is tissue-dependent. Cardiomyocytes, hepatocytes and muscle are genuinely mitochondria-rich (20–30% can be normal). And in single-nucleus RNA-seq you expect near-zero %MT, so the metric flips meaning entirely.

§4 (cells 48–54) — Computing QC metrics

sc.pp.calculate_qc_metrics(ETV6_RUNX1_1, qc_vars=["mt"], inplace=True,
                           percent_top=[20], log1p=True)
qc_vars=["mt"] — "also compute stats for the gene set flagged by the mt column"
inplace=True — modify the object directly rather than returning a copy
percent_top=[20] — compute what fraction of a cell's counts come from its top 20 genes
log1p=True — also add log-transformed versions of the metrics
This adds columns to .obs (per cell):

Column	Meaning	What extreme values suggest
total_counts	total UMIs in the cell — its "library size"	Low: empty droplet / dying cell / poor capture. High: doublet, or a genuinely large/active cell
n_genes_by_counts	how many distinct genes were detected (count > 0)	Low: empty or degraded. High: doublet
pct_counts_mt	% of counts from mitochondrial genes	High: dying / membrane-compromised
pct_counts_in_top_20_genes	% of the library taken by the 20 most-expressed genes	High: low-complexity library — either a dying cell whose transcriptome collapsed, or a specialist cell (e.g. a red cell full of haemoglobin, a plasma cell full of immunoglobulin)
log1p_*	log(x + 1) versions of the above	used for the statistics below
And to .var (per gene): n_cells_by_counts, total_counts, mean_counts, pct_dropout_by_counts.

Why log1p and not plain log? Because log(0) is undefined (negative infinity). log1p(x) = log(x+1) maps 0 → 0 cleanly. This tiny trick appears everywhere in single-cell analysis.

§4.1 (cells 55–76) — Filtering barcodes
The plots (cells 57–61)

sns.displot(ETV6_RUNX1_1.obs, x="total_counts", bins=100)
sns.displot(ETV6_RUNX1_1.obs, x="pct_counts_mt", bins=100)
sns.scatterplot(ETV6_RUNX1_1.obs, x="total_counts", y="n_genes_by_counts", hue="pct_counts_mt")
Because .obs is just a pandas DataFrame, ordinary plotting libraries work on it directly. That's a real advantage of the AnnData design.

How to read that scatterplot — this is the single most informative QC plot and you should be able to interpret it cold:

x = total counts, y = genes detected, colour = %MT
Healthy cells form a tight, slightly saturating curve (more counts → more genes, with diminishing returns as you saturate the transcriptome)
Points bottom-left = low counts and few genes = empty droplets / debris
Points below the curve (many counts but few genes) = low-complexity libraries; usually red/high %MT = dying cells whose diverse mRNA has leaked away leaving only mitochondrial and a few abundant transcripts
Points far top-right = suspiciously high on both = candidate doublets
sc.pl.violin (cell 59, 61) is scanpy's own version. A violin plot is a histogram mirrored and rotated — the width at any height is the density of cells at that value. multi_panel=True puts three metrics side by side.

The outlier function (cell 63) — the statistical core

def is_outlier(adata, metric: str, nmads: int):
  M = adata.obs[metric]
  outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
      np.median(M) + nmads * median_abs_deviation(M) < M
  )
  return outlier
Reading the Python:

def name(args): defines a function; metric: str is a type hint — documentation for humans, not enforced
M = adata.obs[metric] grabs that column
| is element-wise OR on a whole column at once (vectorised — no loop). In pandas you must use | and &, never the words or/and
The condition is: below median − n·MAD OR above median + n·MAD
What is MAD, and why not standard deviation?

MAD = Median Absolute Deviation = median(|xᵢ − median(x)|). It's a measure of spread, like standard deviation — but built from medians.

The point is robustness. Standard deviation is computed from the mean, and both are dragged around by the very outliers you're trying to detect. A handful of enormous doublets inflates the SD, which widens your threshold, which lets the doublets through. Median and MAD barely move no matter how extreme the outliers are — MAD has a 50% breakdown point, meaning up to half your data could be garbage and the estimate still holds. You're using the good cells to define what "normal" is, instead of letting the bad ones define it.

This is a strong interview answer: "data-driven, robust thresholds instead of hard-coded ones, so the filter adapts per sample rather than assuming every sample has the same quality."

Applying it (cells 65–72)

ETV6_RUNX1_1.obs["counts_outlier"]   = is_outlier(ETV6_RUNX1_1, "log1p_total_counts", 5)
ETV6_RUNX1_1.obs["genes_outlier"]    = is_outlier(ETV6_RUNX1_1, "log1p_n_genes_by_counts", 5)
ETV6_RUNX1_1.obs["topgenes_outlier"] = is_outlier(ETV6_RUNX1_1, "pct_counts_in_top_20_genes", 5)
ETV6_RUNX1_1.obs["mito_outlier"]     = is_outlier(ETV6_RUNX1_1, "pct_counts_mt", 3) | (ETV6_RUNX1_1.obs["pct_counts_mt"] > 8)
Three design decisions worth understanding:

1. Why log-transformed counts (cell 66)? Library sizes are extremely right-skewed — most cells have a few thousand UMIs, a few have tens of thousands. A symmetric ± threshold on a skewed distribution is wrong: it would flag far too many cells on the high side and almost none on the low side. Logging makes the distribution roughly symmetric, so a symmetric rule becomes appropriate.

2. Why 5 MADs for most, but 3 for mitochondria? 5 MADs is deliberately permissive — you only remove genuinely extreme barcodes. The philosophy is that over-filtering is worse than under-filtering: aggressive QC can delete an entire rare cell type (small cells with genuinely low RNA content, like resting lymphocytes or platelets) and you'd never know it was missing. High %MT, though, is a much more reliable indicator of a truly dead cell, so it gets the stricter 3 MADs.

3. The extra | (pct_counts_mt > 8) — a hybrid rule: adaptive or an absolute biological ceiling. Reason: if a sample is uniformly bad, the median %MT itself is high, and the purely relative MAD rule would happily keep dying cells because they look "normal for this sample". The hard 8% cap is a floor of biological sanity that the relative rule can't argue away. Combining adaptive and absolute thresholds like this is a genuinely sophisticated touch.


ETV6_RUNX1_1.obs["outlier"] = (
    ETV6_RUNX1_1.obs["counts_outlier"] | ETV6_RUNX1_1.obs["genes_outlier"]
    | ETV6_RUNX1_1.obs["topgenes_outlier"] | ETV6_RUNX1_1.obs["mito_outlier"]
)
Union, not intersection — fail any test and you're out. Note the most recent commit on this repo (a4bc3e4b Correct omission of counts_outlier to overall outliers) fixed exactly this line: counts_outlier had been left out of the OR, so those cells were flagged but never actually removed. A nice illustration of how quietly these bugs hide.


ETV6_RUNX1_1 = ETV6_RUNX1_1[~ETV6_RUNX1_1.obs["outlier"], :].copy()
~ is element-wise NOT — "keep the rows where outlier is False". [rows, :] = selected rows, all columns. .copy() again.

Cell 76 makes the conceptual promotion explicit: after this filtering, we stop calling them barcodes and start calling them cells.

§4.2 (cells 77–81) — Filtering genes

ETV6_RUNX1_1.var["total_counts"].eq(0).value_counts()   # how many genes are never seen
sc.pp.filter_genes(ETV6_RUNX1_1, min_counts=0)
Cell 78 explains the asymmetry: gene filtering is far less critical than cell filtering, because downstream steps (highly-variable-gene selection, PCA) will effectively ignore uninformative genes anyway. But dropping never-detected genes saves memory and avoids divide-by-zero problems later.

The reference genome has ~36,000 genes; in any one sample maybe 15,000–20,000 are detected at all. The rest are genes for cell types not present in your tissue.

⚠️ Something to flag if you present this notebook. Cell 81 uses min_counts=0, which means "keep genes with ≥ 0 counts" — i.e. it filters nothing. Compare cell 103, where the workflow function uses min_cells=1 (keep genes seen in at least one cell), which is the intended behaviour. The commit 8f0184a3 Change gene filtering minimum count to 1 shows this was being corrected, and cell 81 appears to have been left behind. Spotting an inconsistency like this in your own pipeline is exactly the kind of thing that lands well in an interview.

§4.3 (cells 82–94) — Doublet removal
The problem. Two cells occasionally end up in the same droplet, get the same barcode, and are reported as one "cell" with a merged transcriptome. Rate scales with loading concentration — roughly 0.8% per 1,000 cells recovered, so at 10,000 cells you expect ~8% doublets. That is a lot.

Why this is genuinely dangerous, not just noisy: a T-cell/B-cell doublet expresses both CD3 and CD19. It doesn't sit inside either cluster — it forms its own little cluster in between, and you might publish it as a novel intermediate or transitional cell type. Doublets are a well-known source of spurious "novel cell states" in the literature. That's the answer to give if asked why it matters.

Also note simple filtering cannot catch them all: a doublet of two small cells has a perfectly ordinary total count. You need a dedicated algorithm.


sc.pp.scrublet(ETV6_RUNX1_1)
How Scrublet works (Wolock et al. 2019) — the logic is elegant and easy to explain:

Take your real cells and randomly add pairs of them together to create thousands of simulated doublets. You know these are doublets because you made them.
Merge simulated doublets into the real dataset, do PCA and build a k-nearest-neighbour graph.
For each real cell, ask: what fraction of my nearest neighbours are simulated doublets? That's the doublet_score.
A real cell sitting in a neighbourhood dense with simulated doublets is probably itself a doublet. Threshold the bimodal score distribution → predicted_doublet.
Adds two columns to .obs: doublet_score (continuous) and predicted_doublet (boolean).


ETV6_RUNX1_1.obs["predicted_doublet"].value_counts()
ETV6_RUNX1_1 = ETV6_RUNX1_1[~ETV6_RUNX1_1.obs["predicted_doublet"], :].copy()
Critical methodological point (cell 92): doublet detection must be run per sample, before merging samples. A "doublet" can only form within a single droplet emulsion — pairing cells from two different 10x runs is physically meaningless and would produce nonsense simulations. Notice that in the workflow function (cell 103), sc.pp.scrublet sits inside the per-sample function, before sc.concat. That ordering is deliberate.

Cell 92 also raises the alternative: instead of deleting doublets now, cluster first and then drop whole clusters with high mean doublet score. That's often safer — a cluster-level signal is more robust than a per-cell call. And you can run several detectors and take the consensus.

Cell 94 lists alternatives: DoubletFinder (R/Seurat — best performer in the Xi & Li 2021 benchmark), Solo (deep learning, scvi-tools), DoubletDetection.

§5 (cells 95–100) — Normalisation
Why you must normalise
Cell A has 20,000 total UMIs; cell B has 5,000. Gene X shows 100 counts in A and 50 in B. Is X higher in A?

No — it's 0.5% of A's library versus 1.0% of B's. In relative terms it's twice as high in B. The raw difference in library size is technical: differences in cell size, lysis efficiency, capture efficiency, sequencing depth. If you don't correct for it, your first principal component will just be "library size" and your clusters will separate by sequencing depth rather than by cell type.

Layers (cell 96–97)

ETV6_RUNX1_1.layers["counts"] = ETV6_RUNX1_1.X.copy()
.layers is a dictionary of matrices, all the same shape as .X. This preserves the raw integer counts under the name "counts" before anything overwrites .X.

Always keep the raw counts. Differential expression tools that model counts properly (DESeq2, edgeR, and scanpy's own count-based methods) require raw integers — feeding them log-normalised values is statistically invalid, because the whole point of those models is the mean–variance relationship of counts. If asked "why keep raw counts after normalising?", that's the answer.

Method 1 — Pearson residuals

ETV6_RUNX1_1.layers["pearson"] = ETV6_RUNX1_1.X.copy()
sc.experimental.pp.normalize_pearson_residuals(ETV6_RUNX1_1, layer="pearson")
ETV6_RUNX1_1.layers["pearson"] = csr_matrix(ETV6_RUNX1_1.layers["pearson"])
Model each count with a negative binomial (the standard distribution for over-dispersed count data — the variance exceeds the mean, unlike Poisson) whose expected value depends on the cell's total counts and the gene's overall abundance. Then replace each count with its residual:

$$z = \frac{\text{observed} - \text{expected}}{\sqrt{\text{variance}}}$$

So the value is no longer "how much expression" but "how many standard deviations above or below what I'd expect for this gene in a cell of this depth". Values can be negative. Lause et al. 2021 argue this captures biological signal better, especially for lowly-expressed genes that log-normalisation tends to over-shrink.

csr_matrix(...) converts back to sparse: the residual computation produces a dense matrix (every entry filled, including all the zeros that now have non-zero residuals), which is a memory disaster at this scale.

Method 2 — Shifted logarithm (the one actually used)

ETV6_RUNX1_1.layers["logcounts"] = ETV6_RUNX1_1.X.copy()
sc.pp.normalize_total(ETV6_RUNX1_1, layer="logcounts", target_sum=10000)
sc.pp.log1p(ETV6_RUNX1_1, layer="logcounts")
Two steps:

1. normalize_total(target_sum=10000) — divide every count in a cell by that cell's total, then multiply by 10,000. Now every cell has exactly 10,000 total counts, and values are "counts per 10 thousand" (CP10K). Library-size differences are gone. (The commit 5de999e3 Add a target sum when normalising added this explicitly — without target_sum, scanpy defaults to the median library size of the dataset, which makes values dataset-dependent and harder to compare across analyses. Being explicit is better.)

2. log1p — take log(x + 1). Three reasons:

Variance stabilisation. In count data, highly-expressed genes have much larger variance than lowly-expressed ones. Downstream methods (PCA, Euclidean distances, clustering) assume roughly comparable variance across features; without logging, a handful of very high genes dominate everything.
Compresses the dynamic range. Counts span 0 to tens of thousands; logging pulls extreme values in so single outlier cells don't drive the whole analysis.
Differences become ratios. A difference of 1 on the log scale is a constant fold-change regardless of baseline — which is what "fold change" means biologically. This is also why later log-fold-changes are interpretable.
And the +1 handles the zeros, as before.

Comparing them (cell 99)

sns.pairplot(
  pd.DataFrame({"Raw counts": np.nansum(ETV6_RUNX1_1.layers["counts"].toarray(), 1),
                "Log-normalised counts": np.nansum(ETV6_RUNX1_1.layers["logcounts"].toarray(), 1),
                "Pearson residuals": np.nansum(ETV6_RUNX1_1.layers["pearson"].toarray(), 1)})
)
.toarray() converts sparse → dense (needed because np.nansum doesn't handle sparse). np.nansum(M, 1) sums along axis 1 = across genes = per-cell totals. pairplot shows every variable against every other, with distributions on the diagonal.

The result (cell 100) and why it makes sense: log-normalised totals correlate well with raw totals; Pearson residuals do not. That's not a bug — it's the definition. Residuals measure deviation from expectation given depth, so they've deliberately removed the depth information that raw totals consist of.

The notebook then chooses log-normalisation because it's much cheaper computationally and performs well in practice. Cell 100's framing is the honest one and worth repeating in an interview: there is no single correct normalisation; try more than one and check whether your conclusions survive.

§6 (cells 101–118) — Wrapping it into a workflow

def sc_preprocess(path: str, gene_annot: pd.DataFrame):
    print("1. Reading data matrix")
    adata = sc.read_10x_mtx(path)
    ...
    return adata
Everything above, in order, in one function. The print() statements are progress logging — useful because each sample takes minutes and you want to see where you are.

Note the improvements over the interactive version:

gene_annot.drop_duplicates(subset="gene_ids") — the fourth merge safeguard
sc.pp.filter_genes(adata, min_cells=1) — the corrected gene filter
Scrublet runs inside the per-sample function (per the point above)
Normalisation writes to .X directly here (with raw preserved in layers["counts"])
Why functionise? Guarantees every sample gets identical treatment. If you copy-paste per sample, you will eventually change a threshold in one and not the others, and the resulting batch artefact will be indistinguishable from biology. This is a real reproducibility argument, not just tidiness.

Sample sheet and loop (cells 105–107)

sample_info = pd.read_table(".../sample_sheet.tsv")
sample_info = sample_info.rename(columns={"Sample": "sample_id", "SampleName": "sample_name",
                                          "SampleGroup": "sample_group"})
sample_info = sample_info[sample_info["sample_id"].isin(["SRR9264343", "SRR9264347"])]
Reads the metadata table (which sample is which ALL subtype), standardises column names, and — for tractability on a course VM — restricts to two samples.


samples = {}
for s in sample_info["sample_id"]:
    print(f"\nReading {s}")
    samples[s] = sc_preprocess(f"{prefix_inputs}/{s}/outs/filtered_feature_bc_matrix/", h38_genes)
samples = {} is an empty dictionary — a lookup table of key → value. Here: sample ID → its processed AnnData. f"..." is an f-string: anything in {} gets substituted with the variable's value.

Rescuing gene metadata (cells 108–109)

all_var = [x.var for x in samples.values()]
all_var = pd.concat(all_var, join="outer")
all_var = all_var[["gene_ids", "gene_name", "chrom"]]
all_var = all_var[~all_var.duplicated()]
Cell 108 explains the reason: sc.concat() discards .var. It's a known scanpy behaviour (the linked scverse forum thread) — because different samples may have different gene sets after filtering, and there's no unambiguous way to merge conflicting per-gene metadata. So you save it first and put it back afterwards.

The four lines: collect every sample's .var; stack them; keep only the three columns that are sample-independent (a gene's name and chromosome don't vary by sample — unlike total_counts, which does); drop duplicate rows so each gene appears once.

Concatenation (cells 111–115)

adata = sc.concat(samples, join="outer", label="sample_id", index_unique="-")
join="outer" — keep the union of genes across samples, filling missing with zeros. "inner" would keep only genes present in all, losing information.
label="sample_id" — creates an .obs column recording which sample each cell came from. Essential — without it you cannot assess or correct batch effects later.
index_unique="-" — appends the sample name to each barcode (AAACCTGAGACAATAC-1-SRR9264343). Necessary because the same barcode sequence occurs in every 10x run; without this you'd get duplicate cell names.

adata.var = all_var.loc[adata.var_names]       # restore gene metadata, correctly ordered
adata.obs = (adata.obs
             .merge(sample_info, how="left", on="sample_id")
             .set_axis(adata.obs.index))       # attach sample-level metadata to every cell
The second block joins the sample sheet onto the cells, so each cell now knows its subtype (ETV6-RUNX1, PBMMC, …). Cell 116 spells out the .set_axis() requirement again — merge destroys the index, and an AnnData whose .obs index doesn't match its barcodes is invalid and will error or corrupt downstream.


adata.write(f"{prefix_outputs}/caron_filtered_full.h5ad")
.h5ad is HDF5-based — the standard on-disk format for AnnData. Saves the matrix, all layers, .obs, .var, everything, in one file. Next notebook loads it with sc.read_h5ad() and starts from there.

Things likely to come up in an interview
"Why filter cells before normalising?"
Normalisation factors are estimated from the data. Dying cells and doublets have distorted library sizes, so including them biases the size factors for everyone. Clean first, then normalise.

"Why is this order — ambient → QC → doublets → normalise?"
Each step assumes the previous one is done. Doublet detection uses PCA/kNN, which is distorted by debris, so debris goes first. Normalisation assumes the cells are real cells. Ambient correction needs the raw matrix including empties, so it must come before anything filters them away.

"What's the risk of over-filtering?"
Losing real biology. Small quiescent cells legitimately have low RNA content; a strict total_counts floor deletes them, and you'll conclude that cell type isn't present. Undetected, irreversible, and invisible in your final figures. This is why 5 MADs rather than something aggressive.

"How would you check your QC was reasonable?"
Plot the metrics before and after. Compare how many cells each sample lost — if one sample loses 60% and the others 10%, that sample has a problem and may need excluding entirely. And check that cells removed by QC don't form a coherent cluster with sensible marker genes, which would mean you just deleted a real cell type.

"What comes next?"
Feature selection (highly variable genes) → scaling → PCA → neighbourhood graph → UMAP → clustering → annotation → batch integration → DE. That's notebooks 03 onwards in this repo.

"What are the limitations of this pipeline as written?"
Ambient RNA correction was skipped for compute reasons — real contamination remains. Doublets are removed by a single method rather than a consensus. Only two samples are processed here. And the hard 8% MT ceiling is a human-blood-appropriate choice that would need revisiting in another tissue.

Two concrete things you might want to fix before showing this notebook: min_counts=0 in cell 81 (filters nothing — should be min_cells=1 to match the workflow function), and the hard-coded /mnt/e/... WSL paths in cells 19, 26, 42 and 105, which won't run on anyone else's machine.

Want me to fix those, or turn this walkthrough into a formatted study page you can pull up on your phone before the interview?

good give me more questions that they can ask in my interview especiall;y considering they have two bioinformaatians
Let me check what the rest of the course covers, so the questions match what you'd be expected to know.

Bash Extract headers and scanpy calls from downstream notebooks
IN
cd "e:/RNA-seq-courses/scRNA-seq_python_course_2024" && for f in Demonstrations/03_*.ipynb Demonstrations/04_*.ipynb Demonstrations/05_*.ipynb Demonstrations/06_*.ipynb Demonstrations/07_*.ipynb; do echo "##### $f"; python -c "
import json,sys,re
nb=json.load(open(sys.argv[1],encoding='utf-8'))
for c in nb['cells']:
    s=''.join(c['source'])
    if c['cell_type']=='markdown':
        for line in s.split('\n'):
            if line.strip().startswith('#'): print(line.strip())
    else:
        for m in re.findall(r'(sc\.[a-z]+\.[a-z_0-9]+|sce?\.[A-Za-z_]+)\(', s): print('   >',m)
" "$f" 2>/dev/null | awk '!seen[$0]++'; done

OUT
##### Demonstrations/03_yuyao_feature_selection_dim_reduction_demo.ipynb
## Feature Selection and Dimensionality Reduction Demonstrations
### scRNA-seq analysis in Python course 2024, EBI
#### Preparation
#### Gene variance & HVGs
   > sc.pp.normalize_total
   > sc.pp.log1p
   > sc.pp.highly_variable_genes
   > sc.pl.highly_variable_genes
   > sc.pl.violin
#### Dimensionality reduction by PCA
   > sc.tl.pca
   > sc.pl.pca_variance_ratio
   > sc.pl.pca
#### Visualisation by UMAP
   > sc.pp.neighbors
   > sc.tl.umap
   > sc.pl.umap
#### We can see that the different samples have quite a large batch effect
#### So lets proceed to batch correction and data integration
### Exercises
##### Demonstrations/04_anna_clustering_cell_annotation_demo.ipynb
# Clustering tutorial
   > sc.pl.embedding
   > sc.pp.neighbors
   > sc.tl.leiden
   > sc.tl.umap
   > sc.pl.umap
   > sc.pp.normalize_total
   > sc.pp.log1p
   > sc.pl.violin
# Compute marker genes and annotate clusters
   > sc.tl.rank_genes_groups
   > sc.pl.rank_genes_groups_dotplot
   > sc.get.rank_genes_groups_df
##### Demonstrations/05_yuyao_data_integration_demo.ipynb
## Batch correction and data integration Demonstrations
### scRNA-seq analysis in Python course 2024, EBI
#### Preparation
#### Load files and create concatenated object
### Data integration
#### First check unintegrated data
   > sc.pp.normalize_total
   > sc.pp.log1p
   > sc.pp.highly_variable_genes
   > sc.tl.pca
   > sc.pl.pca_variance_ratio
   > sc.pp.neighbors
   > sc.tl.umap
   > sc.pl.umap
#### Batch-aware feature selection
#### Example: integration with harmony
### Exercises
##### Demonstrations/06_andrian_differential_expression_demo.ipynb
# Differential expression analysis with edgeR
## Setup
   > sc.pl.embedding
## Create pseudobulk
## Differential expression analysis
   > sc.pp.normalize_total
   > sc.pp.log1p
   > sc.pp.scale
   > sc.tl.pca
   > sc.pl.pca
### Optional: Gene set analysis
## Putting it all together
### Optional exercises
##### Demonstrations/07_andrian_differential_abundance_demo.ipynb
# Differential abundance analysis with Milo
## Setup
   > sc.pl.embedding
## Differential abundance analysis
   > sc.pp.neighbors
#### Setting up contrast for milo
## Optional exercises
Good — the course runs the full pipeline: HVGs → PCA → UMAP → Leiden clustering → marker genes → Harmony integration → pseudobulk + edgeR DE → Milo differential abundance. Two bioinformaticians will almost certainly push past QC into that territory, and especially into why rather than which function. Here's a bank organised by what they're actually testing.

How bioinformaticians interview differently
Non-specialists ask "what did you do?". Bioinformaticians ask three things, usually disguised:

"Do you understand the assumption, or did you just run the tutorial?" — signalled by why, what happens if, what does that assume.
"Would you notice if it went wrong?" — signalled by how would you check, how do you know that's real.
"Can you actually write code, or only run notebooks?" — signalled by anything about scale, memory, reproducibility, or version control.
The single most valuable habit: when unsure, say what you'd check rather than guessing an answer. "I don't know offhand, but I'd test it by…" scores far better with bioinformaticians than a confident wrong answer. They're evaluating a future colleague, not a quiz contestant.

1. Statistics and methods (their favourite territory)
Q: Your QC uses median ± 5×MAD. Why 5? Justify it.
It's a convention, not a derivation — and say so. For a normal distribution, MAD × 1.4826 ≈ SD, so 5 MADs ≈ 3.4 SD ≈ the outer ~0.07%. The choice is deliberately permissive because the asymmetry of costs is asymmetric: removing a real rare cell type is unrecoverable and invisible, while keeping a few bad cells adds noise that downstream steps tolerate. I'd sanity-check by plotting the retained/removed cells and confirming the removed ones don't form a coherent cluster with sensible markers.

Q: MAD is robust — but what's its weakness?
Two. First, it assumes a unimodal distribution. If your sample genuinely contains two populations with very different RNA content (say lymphocytes and large blasts), the median sits between them and MAD is inflated, so the filter becomes meaningless. Second, if more than half the barcodes are bad — a badly failed run — the median itself describes garbage, and the filter happily certifies garbage as normal. That's precisely why the pipeline adds the absolute pct_counts_mt > 8 backstop.

Q: Why is scRNA-seq count data negative binomial rather than Poisson?
Poisson has variance = mean, which is what pure sampling noise would give. Real data is over-dispersed: variance > mean, because on top of sampling there's genuine biological variation in each cell's expression level (transcriptional bursting, cell cycle, cell size). The NB adds a dispersion parameter to absorb that extra variance. If you model over-dispersed data as Poisson, you underestimate variance and every DE test becomes wildly anti-conservative.

Q: What is a Pearson residual actually telling you?
Observed minus expected, divided by the expected standard deviation — where "expected" comes from an NB model conditioned on the cell's depth and the gene's overall abundance. So the value answers "is this gene higher or lower than a cell of this sequencing depth should show?" rather than "how much of it is there?". This is why residuals don't correlate with raw totals — depth information has been deliberately removed, not lost.

Q: log1p adds a pseudocount of 1. Is that arbitrary?
Yes, and it's a known distortion. A pseudocount shrinks fold-changes for lowly-expressed genes more than for highly-expressed ones, so log-normalisation systematically under-represents low-expression biology. The size of the pseudocount interacts with the target_sum you chose too. This is exactly the critique in Lause et al. 2021 motivating Pearson residuals. In practice log1p works well enough for clustering and visualisation, which is why the field still uses it — but I wouldn't rely on it for quantitative claims about lowly-expressed genes.

Q: You set target_sum=10000. What if you set it to 1e6, or left the default?
target_sum is just the scale, so after logging it's a constant shift — it barely affects clustering. But it interacts with the pseudocount: a larger target sum makes +1 relatively smaller, reducing shrinkage. The default (median library size) makes values dataset-dependent, so the same cell normalised in two different analyses gets two different values. Setting it explicitly makes results comparable and reproducible — which is why the repo's commit history shows it being added deliberately.

Q: Multiple testing — you test 20,000 genes. What do you do?
Benjamini-Hochberg FDR, not Bonferroni; Bonferroni controls family-wise error rate and is far too conservative when you expect many true positives. But the deeper point: in single-cell DE on individual cells, p-values are inflated regardless of correction, because thousands of cells from one patient are pseudoreplicates, not independent samples. Correcting the p-values doesn't fix a wrong null.

2. Pseudobulk and DE (notebook 06 — expect real scrutiny here)
Q: Why does your course use pseudobulk + edgeR rather than a per-cell test like Wilcoxon?
This is the flagship question in modern scRNA-seq and worth memorising. Cells within one patient are not independent replicates — they share that patient's genotype, treatment history, batch, and dissociation. A per-cell test treats 5,000 cells as n=5,000, so the null hypothesis is wrong and the false positive rate is enormous (Squair et al. 2021, Nature Communications, showed most single-cell DE methods have FDRs far above nominal). Pseudobulk sums counts per cell type per sample, giving you n = number of patients, which is the real unit of replication. Then edgeR/DESeq2 apply well-validated bulk RNA-seq NB models. You lose per-cell resolution and gain valid statistics.

Q: So per-cell tests are always wrong?
No — they're fine for marker gene identification within a dataset, where you're ranking genes descriptively to annotate clusters, not making a population-level inference. They're wrong for comparing conditions or patients. That distinction is the answer they want.

Q: You have 4 ETV6-RUNX1 and 3 PBMMC samples. Is that enough for DE?
Barely, and I'd say so explicitly. n=3–4 per group is at the low end; edgeR can fit it, but power for modest effects is poor and the dispersion estimate is unstable (which is why edgeR shares information across genes via empirical Bayes moderation). I'd report effect sizes, not just p-values, be cautious about negative results, and treat findings as hypothesis-generating pending validation.

Q: Why sum counts for pseudobulk rather than averaging normalised values?
Because edgeR/DESeq2 model raw counts and use the library size to weight confidence — a pseudobulk from 2,000 cells should carry more weight than one from 50. Summing preserves that; averaging normalised values throws it away and breaks the NB model's variance assumptions. Same reason you keep layers["counts"] in notebook 02.

Q: A gene comes out top-significant. How do you check it's not an artefact?
Several things: is it mitochondrial or ribosomal (residual QC signal)? Is it driven by one sample — plot per-sample pseudobulk values, not just the p-value. Is it a known ambient contaminant (haemoglobin, in blood data — and remember CellBender was skipped here, so ambient RNA is uncorrected). Is the cell-type composition confounded — if one group has more of a cell type, "DE" within a cluster can reflect a shifted mixture. And does it make biological sense given the ETV6-RUNX1 fusion?

3. Clustering, UMAP and annotation (notebooks 03–04)
Q: Leiden has a resolution parameter. How do you choose it?
There's no ground truth, so I'd be honest: I run a range, and check stability. Practical criteria — do clusters have distinct, interpretable marker genes; do they split/merge sensibly across resolutions (a clustering tree / clustree-style view); does a cluster split produce two genuinely different marker sets or just a gradient? For annotation I'd rather over-cluster and merge based on markers than under-cluster and hide a rare population.

Q: Why Leiden and not Louvain?
Louvain can produce badly connected or even internally disconnected communities — a known failure mode. Leiden (Traag et al. 2019) adds a refinement step that guarantees connected communities and converges to better partitions. Scanpy has deprecated Louvain in favour of Leiden for this reason.

Q: Can you interpret distances in a UMAP?
Mostly no, and this is a favourite trap. UMAP preserves local neighbourhood structure; global inter-cluster distances and the arrangement of clusters on the page are largely arbitrary, and change with n_neighbors, min_dist, and the random seed. So "cluster A is closer to B than to C" is not a supportable claim from a UMAP. It's a visualisation, not a result — I'd verify relationships with PCA distances, correlation of pseudobulk profiles, or a proper trajectory method.

Q: Then why not just cluster on the UMAP coordinates?
Because you'd be clustering on a lossy 2D distortion. Clustering is done on the kNN graph built in PCA space (typically 30–50 PCs), which retains far more structure. UMAP is computed from the same graph purely for display. Note the ordering in the notebooks: sc.pp.neighbors → then both sc.tl.leiden and sc.tl.umap. They're siblings, not sequential.

Q: How many PCs, and why?
Elbow plot (sc.pl.pca_variance_ratio) as a starting point, typically 30–50. But the elbow is subjective, and the honest answer is that clustering results are fairly robust to this within a sensible range — I'd check that my conclusions don't change between, say, 30 and 50. Too few PCs loses rare cell types; too many adds noise.

Q: Why select highly variable genes at all — why not use all 20,000?
Most genes are either not expressed or vary only by technical noise; including them dilutes biological signal and inflates the distance computation with noise. HVG selection (~2,000–3,000) improves signal-to-noise and speeds everything up. The risk: HVG selection is done on the whole dataset, so genes marking a very rare cell type may not make the cut — a real cause of missed rare populations.

Q: Notebook 05 mentions "batch-aware feature selection". Why does that matter?
If you pick HVGs across the pooled dataset, genes that vary because of batch look highly variable and get selected — so you actively select the features that drive your batch effect, then integrate to remove it. Selecting HVGs per batch and taking genes variable in many batches finds genes that vary biologically in a consistent way, which is a much better feature set to integrate on.

Q: How would you annotate a cluster you don't recognise?
Marker genes via rank_genes_groups, cross-referenced against literature and reference atlases (Human Cell Atlas, PanglaoDB, CellTypist for automated label transfer). Then sanity checks: is it a doublet cluster (high doublet score, co-expression of two exclusive lineage markers)? Is it a low-quality cluster (high %MT, low genes — QC leakage)? Is it cell-cycle driven? Only after excluding those would I call it a novel state — and in leukemia I'd also consider it may be a patient-specific malignant clone rather than a cell type, which is a genuinely different kind of entity.

4. Batch effects and integration (notebook 05)
Q: Notebook 03 says "the different samples have quite a large batch effect". In this dataset, how do you know it's batch and not biology?
This is the hardest and most important question in the whole pipeline, and the honest answer is that in this design you partly can't. Each patient is a separate 10x run, so "patient" and "batch" are completely confounded. What you can reason about: technical batch effects tend to affect all cell types similarly and show up in QC-correlated ways (depth, %MT, dissociation-stress genes), whereas biology should be cell-type specific. The strongest evidence is the control cells: PBMMC samples from different individuals should contain the same normal cell types, so residual separation among those is likely technical. Malignant blasts separating by patient is expected biology — each patient's leukemia is a distinct clone.

Q: So what's the risk of integrating this dataset with Harmony?
Over-correction. Integration methods are told "remove differences between these groups", and if the groups differ biologically, the method will obediently delete real signal. Here you could erase the very inter-patient tumour heterogeneity that's the point of the study. The mitigation is to integrate for shared normal cell types and be careful about interpreting malignant populations, and to always compare integrated vs unintegrated embeddings.

Q: How do you evaluate whether integration worked?
Not by eye. Quantitative metrics: batch mixing (kBET, iLISI, graph connectivity) and biological conservation (cLISI, ARI/NMI against pre-integration labels, cell-type silhouette). These trade off against each other — perfect mixing means you deleted biology. scIB (Luecken et al. 2022, Nature Methods) is the standard benchmarking framework, and its main finding is worth citing: no method wins everywhere, and the right choice depends on how strong your batch effect is relative to your biological signal.

Q: Harmony vs scVI vs BBKNN vs Seurat CCA — pick one and defend it.
Harmony is fast, operates on the PCA embedding, and scales well — a good default. But it returns a corrected embedding, not a corrected expression matrix, so you can't do DE on its output. scVI is a deep generative model, handles counts natively, and can output corrected expression, but needs more data and tuning. BBKNN is very fast and modifies only the graph. The key point: integration output is for clustering and visualisation; DE goes back to raw counts with batch as a covariate, or to pseudobulk. Never run DE on integrated values.

Q: What would you have changed about the experimental design?
Multiplexing — cell hashing or genetic demultiplexing (souporcell/vireo) to run multiple patients in one 10x lane. That breaks the patient/batch confound, cuts cost, and lets you distinguish technical from biological variation. Also a shared reference sample across all runs as a technical anchor. Saying this shows you think about design, not just analysis — bioinformaticians rate that highly.

5. Differential abundance / Milo (notebook 07)
Q: What does Milo do that a cluster-based composition test doesn't?
Cluster-based DA tests counts per cluster per sample — which forces you to commit to a clustering, and can only detect changes that align with cluster boundaries. Milo tests overlapping neighbourhoods on the kNN graph instead, so it can detect a shift in a sub-region of a cluster or along a continuum, without discretising first. It uses a NB GLM (edgeR machinery again) on neighbourhood counts, with a spatial FDR correction because neighbourhoods overlap and aren't independent tests.

Q: Why is differential abundance statistically awkward in general?
Because compositions are constrained to sum to 1. If one population genuinely expands, every other population's proportion falls even though nothing happened to it. So you cannot interpret proportion changes independently — this is the compositional data problem. It's also why scRNA-seq can't tell you about absolute cell numbers at all: you sequenced a fixed number of cells, not a tissue.

6. Python and software engineering (the "can you actually code" probe)
Since you're new to Python, don't oversell — but these are all learnable answers, and admitting a level while showing you understand the concepts works well.

Q: What's a sparse matrix and why does it matter here?
~90–95% of a scRNA-seq matrix is zeros. Dense storage of 10,000 cells × 20,000 genes at 4 bytes is ~800 MB for almost nothing. CSR stores only non-zero values plus index arrays — often 10–20× smaller. It also explains a real failure mode: operations that densify (like .toarray(), or Pearson residuals) can blow up memory, which is exactly why cell 97 wraps the result back in csr_matrix().

Q: What's the difference between a view and a copy in AnnData, and why do you care?
Subsetting returns a view — a window onto the original, no data duplicated. Modifying a view either warns or silently fails to propagate, and it keeps the (large) parent object alive in memory. .copy() materialises an independent object. It's the most common source of confusing scanpy bugs.

Q: Your dataset is 500,000 cells and won't fit in RAM. What do you do?
Options ladder: back the AnnData with HDF5/Zarr and use backed mode or AnnData's on-disk access; use Dask-backed arrays; subsample for exploratory work and run the final pipeline on the full set on an HPC node; or use out-of-core-friendly tooling (scVI trains in minibatches; scanpy has incremental PCA). Also: do the per-sample preprocessing in parallel as separate jobs, since it's embarrassingly parallel — exactly the structure of sc_preprocess().

Q: How would you make this analysis reproducible?
Version control (git — this repo already has it), a pinned environment (conda environment.yml or the Singularity container that's sitting in this repo as scrnaseq2024.sif), set random seeds — Leiden, UMAP and Scrublet are all stochastic — pin the reference annotation to a saved file rather than a live Ensembl query (as this notebook does), and avoid hard-coded absolute paths. I'd note the notebook currently has /mnt/e/... paths that would break on any other machine, and I'd parameterise those.

Q: Notebooks or scripts?
Notebooks for exploration and communication; scripts/workflows (Snakemake, Nextflow) for anything you'll run more than once or on more than a few samples. The sc_preprocess() function in cell 103 is the natural boundary — that's the bit that should become a script. Notebooks are bad at version control (JSON diffs) and let you run cells out of order, which silently breaks reproducibility.

Q: R or Python for single cell?
Both, and the honest answer is you go where the method is. Python/scanpy for scalability, scVI, and the scverse ecosystem; R/Bioconductor for edgeR/DESeq2, Milo, and much of the statistical tooling — note this course itself calls edgeR and Milo from a Python workflow. Interoperability via anndata2ri, zellkonverter, or just writing .h5ad/.rds between steps.

7. Biology of this dataset (don't get caught out)
They chose a leukemia dataset. Expect at least one question showing whether you engaged with the biology or just the pipeline.

Q: What is ETV6-RUNX1 and why does it matter?
A t(12;21) translocation fusing ETV6 (a transcriptional repressor) to RUNX1 (a master haematopoietic transcription factor). The fusion acts as a dominant repressor of normal RUNX1 targets, blocking B-cell differentiation. It's the most common recurrent translocation in childhood B-ALL (~25%), often arises in utero as a pre-leukemic clone, requires secondary hits to become overt leukemia, and carries a relatively good prognosis.

Q: What would you expect to see when comparing ALL samples to PBMMC controls?
Controls should show a normal spread of bone marrow mononuclear cells — T cells, B cells, NK, monocytes, some progenitors. The ALL samples should be dominated by a large, relatively homogeneous population of blasts arrested at a specific differentiation stage, plus residual normal cells. Practically: the malignant cells will separate strongly by patient, because each clone is genetically distinct, while the normal cells should mix across patients. That per-patient separation is the signature that tells you which population is malignant.

Q: HHD is "high hyperdiploid". Could you detect that from scRNA-seq?
Yes, indirectly — via inferred copy number from expression, using tools like inferCNV or CopyKAT, which average expression across genomic windows and look for chromosome-scale shifts. Extra chromosome copies produce elevated expression across that whole chromosome. This is also a standard way to distinguish malignant from normal cells when you don't have a clean marker. Notice this requires the chromosome annotation the notebook attached in §3 — a nice link back to the QC step.

Q: Why is %MT thresholding at 8% a reasonable choice here specifically?
Blood and bone marrow cells are relatively mitochondria-poor and come as a native suspension needing no harsh enzymatic dissociation, so baseline %MT is low and a tight threshold is safe. In solid tissue requiring collagenase digestion, or in mitochondria-rich tissue like heart, kidney or muscle, 8% would delete healthy cells. The threshold is a property of the tissue, not a universal constant.

Q: Any concern about the blasts specifically in QC?
Yes — a good point to volunteer. Leukemic blasts are large, highly proliferative cells with high RNA content, so they sit at the high end of total_counts and n_genes. An aggressive upper threshold intended to catch doublets could preferentially delete the malignant population you're studying. Similarly, real red cells / erythroid precursors have genuinely low-complexity, haemoglobin-dominated libraries and get caught by pct_counts_in_top_20_genes. Both are arguments for permissive filtering plus inspecting what you removed.

8. Curveballs and judgement questions
Q: Your collaborator says "I ran the pipeline and got 47 clusters — that's 47 cell types!" What do you say?
That cluster count is a function of resolution, not biology, and clusters aren't cell types until they're annotated and validated. I'd ask what resolution was used, check whether clusters have distinct markers or just gradients, check for QC-driven and doublet clusters, and reframe: the question isn't how many clusters, it's which populations are reproducible and interpretable.

Q: You get a beautiful result. How do you try to break it?
Re-run with a different seed, a different resolution, a different number of PCs, and with/without integration. Hold out a sample and see if it replicates. Check the result isn't driven by one patient. Check it against an independent dataset if one exists. If the finding survives all that, it's worth believing.

Q: What's the biggest limitation of scRNA-seq as a technology?
Pick one and go deep rather than listing. Good options: (a) you measure a dissociated snapshot — you lose spatial context and the dissociation itself induces stress-response genes, so some of what you measure is an artefact of the protocol; (b) dropout means absence of a transcript is not evidence of absence of expression; (c) mRNA is a poor proxy for protein, with correlation often ~0.4.

Q: A reviewer asks why you didn't correct for ambient RNA. Answer.
Straight and unspun: it was skipped for compute reasons (CellBender needs hours per sample and ideally a GPU), which is a real limitation given this is blood data where lysed cells contaminate every droplet. I'd note which conclusions are most vulnerable — genes from abundant fragile cell types appearing in other populations — and say I'd rerun with CellBender or SoupX before publication.

Q: What would you do differently if you designed this study?
Multiplex patients per lane to break the patient/batch confound; include matched remission or germline samples per patient; add CITE-seq surface markers for confident immune cell annotation; consider single-nucleus for frozen samples. And more replicates per subtype — 2 replicates for HHD and PRE-T is too few for confident group-level inference.

9. Questions to ask them
Bioinformaticians notice when you ask technical questions back. These signal you're evaluating the environment like a colleague:

What compute do PhD students have access to — local cluster, cloud, GPU?
Does the group work with wet-lab data generated in-house, or public/collaborator data? (This determines whether you'll have any say in experimental design.)
How are analyses shared and reviewed within the group — is there a code review culture, shared pipelines, a group GitHub?
Are projects mostly method development or method application? Where would this project sit?
Is the group moving towards spatial or multi-omic data?
What does the balance look like between wet and dry work, and would I be expected to generate my own data?
The three sentences to have ready
If you take nothing else in: bioinformaticians repeatedly probe the same three ideas, so have a crisp version of each.

Independence — "cells are not replicates; patients are." (Pseudobulk, DE, why n matters.)
Confounding — "in this design, patient and batch are the same variable, so I have to be careful about what integration removes."
Robustness — "I don't trust a result until I've changed the seed, the resolution, and the parameters and it's still there."
Almost every hard question above is one of those three wearing a different costume.
