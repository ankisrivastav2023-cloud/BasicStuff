Good instinct — this is the notebook that matters most. Feature selection and dimensionality reduction are where scRNA-seq analysis is actually decided; QC is hygiene, clustering is downstream of the choices made here. Let me go cell by cell.

What this notebook is for, in one sentence
You start with ~10,000 cells × ~20,000 genes and you cannot do anything with that. This notebook reduces it to ~2,000 informative genes → 50 meaningful dimensions → a neighbour graph → a 2D picture, and each of those reductions is a decision you have to be able to defend.

Setup (cells 2–5)

from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning, NumbaWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
...
Numba is a just-in-time compiler that converts Python into fast machine code. UMAP and Leiden depend on it heavily — without it they'd be unusably slow. It's noisy with deprecation warnings, so this silences them.

Purely cosmetic, but worth a caveat if asked: blanket warning suppression is a habit that eventually hides something real. Suppressing specific named warning classes, as here, is the defensible version.


prefix_output = r"E:\RNA-seq-courses\scRNA-seq_python_course_2024\Data\results\03"
os.makedirs(prefix_output, exist_ok=True)
The r"..." is a raw string — backslashes are treated literally rather than as escape characters. Without it, \R and \n inside a Windows path would be interpreted as escapes (\n = newline) and the path would silently break. Any Windows Python user needs this reflex.

exist_ok=True means "don't error if the folder is already there" — makes the notebook safe to re-run.

⚠️ Inconsistency to be aware of: this cell uses a Windows path (E:\...) while cell 6 uses a WSL path (/mnt/e/...). Those are two different ways of addressing the same drive from two different environments. It works for you because you've run it, but it means the notebook only runs in one specific setup. In an interview, mentioning that you'd parameterise paths rather than hard-code them reads as real experience.


pd.set_option('display.max_columns', 40)
Pandas truncates wide tables by default. This shows up to 40 columns — necessary because .obs and .var here have many QC columns you actually want to see.


print("cwd:", os.getcwd())
print(glob.glob("../Data/**/Caron_normalized.500.h5ad", recursive=True))
Debugging: print the working directory, then search for the file. ** with recursive=True searches all subdirectories. This exists because relative paths in notebooks depend on where Jupyter was launched from — a constant source of confusion.

Loading and inspecting the data (cells 6–10)

adata = sc.read_h5ad(".../Data_from_rds/Caron_normalized.500.h5ad")
print(adata)
Data_from_rds tells you the provenance: this was converted from an R SingleCellExperiment object (.rds is R's save format). The course runs in both R and Python, and this is the R object exported for Python use. .500 means 500 cells per sample, subsampled so it runs fast.

You can confirm the R origin from the column names later: SampleGroup, subsets_Mito_percent are scater/scran (Bioconductor) naming conventions. Scanpy would have called them sample_group and pct_counts_mt. Small detail, but noticing it is exactly the kind of thing that marks someone who reads data rather than assuming.

Cells 8–10: the most underrated cells in the notebook

adata.X.sum(axis=1)     # cell 8
adata.X.max()           # cell 9
Learn this pattern — it's a genuinely valuable skill and a great interview answer.

The problem: .X is just a matrix of numbers. Nothing in the file tells you reliably whether it's raw counts, normalised, log-transformed, or scaled. The filename says "normalized", but filenames lie — and here it did. So you interrogate the matrix:

Check	What it tells you	Result here
.X.sum(axis=1) — row sums = total per cell	If size-factor normalised, every cell sums to the same value. Different sums = not normalised.	Different → not normalised
.X.max() — largest value	log1p data maxes out around 5–10. Hundreds or thousands = raw counts.	Large → not log-transformed
Two more you should know: are the values integers? Then they're raw counts. Are there negative values? Then they've been scaled (z-scored) or are model residuals.

The comment in cell 8 diagnoses it correctly: "perhaps the conversion method we used did not pick up the normalized data slot and just kept the raw data." An R SingleCellExperiment holds multiple assays — typically counts and logcounts. The R→Python converter grabbed counts and left logcounts behind. Hence cell 12 re-does the normalisation.

Interview question this prepares you for: "You're handed someone else's .h5ad. How do you know what state the data is in?" Answer with the checks above. It's practical, it shows you don't trust metadata blindly, and it's a real thing that goes wrong.

The consequence of not checking is severe and silent: run log1p on already-logged data and you get log-of-log values. Nothing errors. Your clustering just quietly becomes wrong.

EXPLAIN ED IN SIMPLER TERMS: 

Start with what adata.X is
It's just a big table of numbers. Rows = cells, columns = genes.


              CD3E   MS4A1   ACTB
Cell 1         100      0    2000
Cell 2          40      5     800
Cell 3           0     90    1500
That's it. A grid of numbers.

The problem
Those numbers might be raw (straight from the machine), or they might have been cleaned up already by someone else.

Looking at the table, you can't tell which. And your next step depends on knowing.

It's like being handed a pot of soup and asked to add salt. Did someone already salt it? If you can't tell, you might double-salt it and ruin it.

Two ways to tell
Check 1: add up each row

adata.X.sum(axis=1)

Cell 1:  100 + 0 + 2000 = 2100
Cell 2:   40 + 5 +  800 =  845
Cell 3:    0 + 90 + 1500 = 1590
Different totals → raw. Nobody's touched it.

If someone had cleaned it up, every row would add to the same number (usually 10,000), because that's exactly what the cleanup does — it makes all cells equal.

Check 2: find the biggest number

adata.X.max()
Here it's 2000 — big.

There's a second cleanup step called "log" that shrinks big numbers down. After it, nothing is bigger than about 10.

Big number → the log step hasn't happened.

What this told them
Both checks said: raw, untouched.

Even though the file was named Caron_normalized.h5ad — the word "normalized" right there in the name.

The name was just wrong. Somebody named the file, and the actual contents didn't match.

So they fixed it

sc.pp.normalize_total(adata, target_sum=1e4)   # make every row add to 10,000
sc.pp.log1p(adata)                             # shrink the big numbers
Two lines. Now the data really is what the filename claimed.

Why bother checking
Because doing it twice breaks everything — and Python won't warn you.

Back to the soup: salt it twice and it's inedible, but nothing stops you. No error message. You only find out at the end, when your results are wrong and you don't know why.

What to remember
Before you start work on any dataset someone gives you, run these two lines:


adata.X.sum(axis=1)[:5]    # all the same number? → already normalised
adata.X.max()              # smaller than 10?     → already logged
Ten seconds. Then you know what you're holding.

Don't trust the filename. Look at the numbers.

adata.var_names_make_unique()
The genes are indexed by gene symbol, and symbols are not unique — several Ensembl IDs can map to the same symbol (paralogs, patches, readthrough transcripts). AnnData needs a unique index or lookups like adata[:, "CD79A"] become ambiguous. This appends -1, -2 to duplicates.

This is precisely the trade-off from notebook 02: Ensembl IDs are unique but unreadable; symbols are readable but not unique. Real pipelines keep both — IDs as the index, symbols as a column.

Feature selection: highly variable genes (cells 12–16)

sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
sc.pp.log1p(adata)
Same two steps as notebook 02 (CP10K then log), redone because the conversion dropped them. Refresher on why: normalize_total removes library-size differences so cells are comparable; log1p stabilises variance so a few high-expression genes don't dominate.


sc.pp.highly_variable_genes(adata, n_top_genes=int(0.1 * adata.n_vars), flavor="seurat")
Why select genes at all
Of ~20,000 genes, only a minority distinguish cell types:

Many aren't expressed in this tissue at all
Housekeeping genes (ACTB, GAPDH, ribosomal proteins) are similar everywhere
What's left — lineage transcription factors, surface markers, effector proteins — is where cell identity lives
Include everything and the informative genes are diluted by thousands of noisy ones. Distance between cells becomes dominated by noise. Selecting HVGs raises signal-to-noise, cuts compute, and improves every downstream step.

int(0.1 * adata.n_vars) = top 10% of genes. (n_vars = number of genes; int() because it must be a whole number.) Typical practice is 2,000–3,000.

The mean–variance problem — the key concept here
You cannot simply rank genes by variance, because in count data, variance scales with the mean. A gene averaging 100 counts will have far larger absolute variance than one averaging 2, purely as a property of counting statistics — nothing biological.

Rank by raw variance and you just recover "the most highly expressed genes", which are largely housekeeping genes. Exactly the wrong answer.

So you need genes that are variable relative to what their expression level predicts. Cell 13's comment describes how flavor="seurat" does it:

Compute dispersion (variance / mean) for every gene
Bin genes by mean expression
Within each bin, standardise dispersion (subtract the bin's mean, divide by its SD)
Take the top genes by that standardised score
The binning is the trick: a gene only competes against others at a similar expression level. So a moderately expressed transcription factor that's ON in T cells and OFF elsewhere can outrank a highly expressed housekeeping gene.

Flavours — know the difference, it's a common gotcha:

Flavour	Expects	Method
seurat (default)	log-normalised	binned dispersion, as above
seurat_v3	raw counts	variance-stabilising transform, picks top N by standardised variance
cell_ranger	log-normalised	similar binning, different implementation (this is what nb 05 uses with batch_key)
Feed seurat_v3 log data (or seurat raw counts) and it won't error — it'll just return a bad gene list. Silent failure again.


plt.figure(figsize=(10, 6))
sc.pl.highly_variable_genes(adata)
The diagnostic plot. x-axis = mean expression, y-axis = dispersion, with selected genes highlighted. Cell 14's comment tells you what you're checking: "not all HVGs are highly expressed but they are roughly drawn from different expression levels, while having a high dispersion."

Good: selected genes spread across the whole x-axis — you've found variable genes at every expression level
Bad: selected genes clumped at high expression — the mean-correction failed, and you've just selected abundant genes
That's the visual proof the normalisation and the binning did their job.


hvgs = adata.var.index[adata.var.highly_variable]
print(f"Number of HVGs: {len(hvgs)}")
highly_variable_genes doesn't delete anything — it adds a boolean column highly_variable to .var. This line uses that column to pull out the names. Non-destructive by design: all genes stay in the object, and PCA simply chooses to use only the flagged ones.


sc.pl.violin(adata, hvgs[:5], jitter=0.4, alpha=0.05, multi_panel=True)
Sanity check on the top 5. Cell 16's comment states the expected shape: "much higher expression in some of the cells, but not in the majority."

That bimodal pattern — a fat blob near zero plus a separate population high up — is the signature of a genuine cell-type marker: off in most cells, strongly on in one type. A gene that's uniformly moderate in everything is variable for the wrong reasons.

jitter=0.4 scatters points horizontally so they don't overlap; alpha=0.05 makes them 95% transparent so density is visible with thousands of points.

Have ready: "What's the risk of HVG selection?" You can miss rare cell types — a marker for 20 cells out of 10,000 may not register as globally variable. Mitigations: select more genes, sub-cluster and re-select HVGs within a population of interest, or add known markers back manually. Also note HVGs are dataset-dependent: the same cell in a different experiment gets a different feature set, which is one reason results don't always replicate.

PCA (cells 18–21)

sc.tl.pca(adata, use_highly_variable=True, n_comps=100)
What PCA actually does
You have 2,000 genes = 2,000 dimensions. But genes don't vary independently — all the B-cell genes rise and fall together, all the myeloid genes together. The data really lives on a much lower-dimensional surface inside that 2,000-dimensional space.

PCA finds new axes that are weighted combinations of genes, chosen so that:

PC1 captures as much variance as possible
PC2 captures the most of what's left, and is uncorrelated (orthogonal) to PC1
and so on
Each PC has loadings — a weight per gene. In practice PC1 might load heavily on B-cell genes, so a cell's PC1 score is effectively its "B-cell-ness". The PCs are interpretable composite signals, not arbitrary.

PCA is linear — every PC is a straight-line combination of genes. That's its limitation (it can't unfold curved structure, which is why UMAP follows) and its strength (it's deterministic, interpretable, and distances in PCA space are meaningful).

Why it's mandatory, not optional
The curse of dimensionality. In very high dimensions, the distances between all pairs of points converge towards being equally large. "Nearest neighbour" stops meaning anything — and every step after this is built on nearest neighbours. PCA restores meaningful distance. This is the deepest answer to "why reduce dimensions?" and the one that impresses.
Denoising. Real biological signal is correlated across many genes and lands in early PCs. Random technical noise is uncorrelated and spreads thinly across late PCs. Keeping the first 50 keeps signal and discards noise.
Speed. 50 dimensions instead of 2,000.
use_highly_variable=True restricts to the HVGs — the two steps are designed to work together.

n_comps=100 deliberately over-computes (default is 50), as the comment says, so you can see where the variance plateaus rather than assume it.


print(adata.obsm["X_pca"][:10, :5])
.obsm = "observation matrices" — per-cell matrices with more than one column. X_pca is cells × 100. The slice prints the first 10 cells × first 5 PCs, just to show it's an ordinary numeric matrix.

Where things live in AnnData — worth keeping straight:

Slot	Holds	Example
.X	the main matrix	expression
.obs	per-cell columns	SampleGroup, leiden
.var	per-gene columns	highly_variable
.obsm	per-cell matrices	X_pca, X_umap
.varm	per-gene matrices	PCs (the loadings)
.uns	unstructured extras	variance_ratio, plot colours
.layers	alternative .X	counts, logcounts

sc.pl.pca_variance_ratio(adata, n_pcs=70, log=True)
The scree plot — variance explained per PC, in descending order. You look for the elbow: the point where the curve flattens, after which each extra PC adds only noise. The comment reads it as plateauing around 50–60.

log=True puts the y-axis on a log scale, which makes the elbow much easier to see — on a linear scale the first few PCs dwarf everything and the tail looks flat immediately.

Be honest about this if asked: the elbow is subjective, and there's no principled cutoff. What matters is that results are robust across a sensible range — if your clusters change completely between 30 and 50 PCs, that's the finding, not the parameter.


sc.pl.pca(adata,
    color=["SampleGroup", "SampleGroup", "subsets_Mito_percent", "subsets_Mito_percent"],
    dimensions=[(0,1), (2,3), (0,1), (2,3)], ncols=2, size=4)
This is the best cell in the notebook to talk about in an interview.

It plots four panels: SampleGroup on PC1–2 and PC3–4, then %MT on the same. The question being asked is: what is actually driving the largest variation in my data?

The comment records the answer: "different SampleGroup showed significant variation on the PC space, but pct_mito didn't."

Read that properly — it's two separate findings:

%MT does not structure the PCs. This is a QC pass. If it had, it would mean dying cells were driving your principal axes of variation, your notebook-02 filtering was insufficient, and any clusters you found downstream might just be "healthy vs dying". Always check this. Also worth checking: total counts, number of genes detected, and cell-cycle phase, all of which can dominate PC1 if uncontrolled.

SampleGroup does structure the PCs. This is ambiguous, and saying so is the sophisticated answer. It could be real biology (ETV6-RUNX1 blasts genuinely differ from healthy PBMMCs) or technical batch (each sample is its own 10x run). In this dataset those are completely confounded — patient, condition and batch are the same variable — so PCA alone cannot separate them. That ambiguity is exactly what drives you into notebook 05.

dimensions=[(0,1),(2,3)] matters too: always look past PC1–2. Structure that doesn't appear in the first two PCs is common — a rare cell type may only separate on PC7.

The neighbour graph (cell 24)

sc.pp.neighbors(adata, n_pcs=50)
Do not gloss over this line — it's the hinge of the entire pipeline.

For each cell, find its k nearest cells in 50-dimensional PCA space. Build a graph: cells are nodes, "is a near neighbour of" is an edge.

Everything downstream is computed from this one object:

UMAP lays this graph out in 2D
Leiden (notebook 04) finds communities in this graph
Milo (notebook 07) tests neighbourhoods of this graph
So UMAP and clustering are siblings computed from a shared parent, not sequential steps. That's why clusters usually look coherent on a UMAP — same underlying graph — and it's also why you must never cluster on UMAP coordinates. Being able to state this relationship cleanly is a strong signal you understand the pipeline rather than the button order.

n_pcs=50 comes straight from the scree plot — the comment makes the reasoning explicit, including the trade-off: more PCs capture more structure in complex tissue, but late PCs are mostly noise.

Scanpy's default n_neighbors is 15. It controls local-versus-global sensitivity: small values (5) give fragmented, very local structure; large values (50+) smooth everything and can absorb rare populations.

The graph is also weighted and pruned — scanpy uses UMAP's fuzzy simplicial set construction, so edges carry a similarity weight rather than being simply on/off. That's why the graph is a reasonable object for community detection.

UMAP (cells 25–27)

sc.tl.umap(adata)
What UMAP does, conceptually: treat the kNN graph as a set of springs — connected cells attract, unconnected cells repel — and let it settle in 2D. More formally it builds a fuzzy topological representation of the high-dimensional data and finds a low-dimensional layout with a similar structure, optimised by cross-entropy.

It is non-linear, unlike PCA — it can unfold curved manifolds that PCA would flatten. That's why it produces interpretable pictures where PCA plots look like a blob.


sc.pl.umap(adata, color="SampleName", size=3)   # cell 26
sc.pl.umap(adata, color="SampleGroup", size=3)  # cell 27
Colour by individual sample, then by group. Cell 28's conclusion — "the different samples have quite a large batch effect" — comes from seeing samples occupy separate territory rather than intermingling within shared clusters.

The visual rule: if a cell type is present in all samples, you expect one cluster containing all colours mixed. Getting one cluster per sample instead means something is separating them — batch, or biology, and in this design you can't tell which from the picture alone. That's why rigorous work uses quantitative metrics (kBET, LISI) rather than eyeballing.

The thing you must be able to say about UMAP
Cell 23 states it explicitly: "The results of a UMAP projection should be used for visualisation only and not for downstream analysis (such as cell clustering)." This is the single most reliable trap question in single-cell interviews, and the notebook's Exercise 2 is exactly it:

If the Euclidean distance between cluster A's centroid and B's is smaller than to C's on the UMAP, are A and B transcriptomically more similar?

No. UMAP optimises local neighbourhood preservation. Between-cluster distances, and the arrangement of clusters on the page, are largely arbitrary — they change with n_neighbors, min_dist, and the random seed. To make a similarity claim you'd use distances in PCA space, correlation between cluster pseudobulk profiles, or a dedicated trajectory method.

Cell 23 also flags the honest caveat about the marketing: UMAP claims better global structure preservation than t-SNE, "but there has been debates about whether this property is as good as it seems." That's correct — Chari & Pachter's "The Specious Art of Single-Cell Genomics" argues 2D embeddings distort so severely that the distortion often exceeds the signal. Knowing that debate exists is a genuine differentiator.

Exercise 1 (n_neighbors = 5, 50, 500) is the practical demonstration: at 5 the UMAP shatters into fragments, at 500 it becomes smooth blobs. The picture changes dramatically and the underlying data never did. Run it — it's the fastest way to internalise why you shouldn't over-read a UMAP.

Worth knowing the alternatives: t-SNE (older, purely local, slower, no global claims), PCA (linear, honest distances, poor visual separation), and the emerging view that you should verify any UMAP-based claim in the space you actually computed.

Saving (cell 29)

adata.write_h5ad(f'{prefix_output}/03_dim_reduced_Caron_normalized.500.h5ad')
Writes the whole object — .X, .obs, .var, .obsm (with X_pca and X_umap), .uns — to HDF5. The expensive computations are now cached, so the next notebook doesn't recompute them.

(As established earlier: the course's next notebook actually loads its own provided checkpoint rather than this file. Your output is for your own continuity.)

The five things to have ready for the interview
"Why select highly variable genes?" — Most genes don't distinguish cell types; including them buries signal in noise. And you can't rank by raw variance because variance scales with the mean, so you rank by dispersion within expression bins.

"Why PCA before clustering?" — The curse of dimensionality: in 2,000 dimensions all distances converge, so nearest-neighbour becomes meaningless — and everything downstream is nearest-neighbour based. PCA also denoises, since signal is correlated and lands in early PCs while noise is uncorrelated and spreads across late ones.

"How many PCs?" — Scree plot elbow, ~50 here. But it's subjective, so what matters is that conclusions are robust across a reasonable range.

"Can I interpret this UMAP?" — Local structure yes, between-cluster distances no. Changes with parameters and seed. It's a visualisation, not a result.

"What did the PCA diagnostic tell you?" — %MT didn't drive the PCs, so QC was adequate. SampleGroup did — but patient, condition and batch are confounded in this design, so PCA can't tell you which, and that's what motivates integration.

A good way to signal depth without overclaiming: when you describe this notebook, frame it as three deliberate reductions, each with a cost — genes (lose rare markers), dimensions (lose late-PC structure), 2D (lose almost everything quantitative). Someone who names the cost of each step is obviously thinking, not reciting.

Want me to fix that cell 27 syntax error so you can run the n_neighbors exercise cleanly? It's genuinely the fastest way to make the UMAP point stick.
