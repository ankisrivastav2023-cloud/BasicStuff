e:\RNA-seq-courses\scRNA-seq_python_course_2024\Demonstrations\02_iris_qc_normalisation_yuyao_mod_Feb2025_demo.ipynb
what is hue?
hue is a seaborn parameter that means: colour the plot by this column.

The word "hue" just means colour. You hand it a column name, and seaborn gives each distinct value in that column its own colour, plus a legend.

From your own notebook

sns.displot(ETV6_RUNX1_1.obs, x="pct_counts_mt", hue="mito_outlier", multiple="stack")
                              └── what to plot ──┘  └─ what to colour by ─┘
mito_outlier holds only True or False, so you get two colours:


      │██                          ██ = False (kept)
count │██                          ▓▓ = True  (flagged)
      │██▓
      │██▓▓▓▓▓▓▓▓▓
      └──────────────── pct_counts_mt
        low        high
Now you can see where the cutoff fell — the colour change marks the boundary. Without hue you'd get one solid grey histogram and no idea which cells your threshold caught.

The other example you had

sns.scatterplot(ETV6_RUNX1_1.obs, x="total_counts", y="n_genes_by_counts",
                hue="pct_counts_mt")
Here pct_counts_mt is a continuous number, not True/False. Seaborn notices and switches from distinct colours to a colour gradient — light for low mito %, dark for high — with a colour bar instead of a legend.

Column type	What you get
True/False, or categories	separate colours + legend
Continuous numbers	gradient + colour bar
That plot is doing real work: it puts three QC metrics on one picture. Position tells you counts and genes, colour tells you mitochondrial content. The bad cells cluster in one corner — low counts, low genes, dark colour.

The general idea
hue adds a third variable to a 2D plot by using colour as the extra dimension. It works across most seaborn functions — scatterplot, displot, violinplot, boxplot.

Scanpy's equivalent is spelled color=:


sc.pl.scatter(ETV6_RUNX1_1, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
sc.pl.umap(adata, color="leiden")
Same concept, different keyword — seaborn says hue, scanpy says color.

e:\RNA-seq-courses\scRNA-seq_python_course_2024\Demonstrations\02_iris_qc_normalisation_yuyao_mod_Feb2025_demo.ipynb
explaain grap[hs in step 25 and what are they comntributing
Bash List notebook cells with index and execution count
IN
cd "e:/RNA-seq-courses/scRNA-seq_python_course_2024" && python -c "
import json
nb=json.load(open('Demonstrations/02_iris_qc_normalisation_yuyao_mod_Feb2025_demo.ipynb',encoding='utf8'))
for i,c in enumerate(nb['cells']):
    src=''.join(c['source']).strip().replace(chr(10),' | ')[:110]
    ec=c.get('execution_count')
    print(f'{i:3d} [{str(ec):>4}] {c[\"cell_type\"][:4]}: {src}')
" 2>&1 | head -70

OUT
  0 [None] mark: # Counts Matrix QC and Normalisation
  1 [None] mark: This notebook starts from a high-throughput single-cell RNA-seq data counts matrix, the product of a previous 
  2 [None] mark: ## Table of contents   | 0. [Background](#background) | 1. [Ambient RNA correction](#rna_correction) | 2. [Rea
  3 [None] mark: ## 0. Background
  4 [None] mark: ### Data
  5 [None] mark: We are working with data from [Caron et al. 2020](https://doi.org/10.1038/s41598-020-64929-x), investigating t
  6 [None] mark: ### Setting up
  7 [   1] code: import os | from glob import glob | import numpy as np | from scipy.stats import median_abs_deviation | from s
  8 [None] mark: <a id="rna_correction"></a> | ## 1. Ambient RNA correction
  9 [None] mark: We remove ambient RNA from the dataset, typically using the raw dataset with empty droplets.  The empty drople
 10 [None] mark: However, because of VM and time limitations, we will not run this step.  To give you an idea of runtime, the P
 11 [None] mark: First, we create an anndata h5ad file from the raw CellRagner output.  In this example, let us use one of the 
 12 [None] raw: # These paths are specific to this demo only, where the root path is the root of the course materials repo. | 
 13 [None] mark: We then go to a bash terminal to run CellBender.  CellBender usage is described [here](https://cellbender.read
 14 [None] raw: #!/usr/bin/env bash |  | PREFIX_INPUTS="../Data/CellRanger_Outputs" | PREFIX_OUTPUTS="../Data/results/02" |  |
 15 [None] mark: Several outputs will be written, but below are main outputs. | 1. `$PREFIX_OUTPUTS/cellbender/SRR9264343.denoi
 16 [None] raw: ETV6_RUNX1_1 = sc.read_10x_h5("$PREFIX_OUTPUTS/cellbender/SRR9264343.denoised_filtered.h5")
 17 [None] mark: <a id="reading_data"></a> | ## 2. Reading data
 18 [None] mark: As we do not run CellBender in this demo, we use the filtered outputs of CellRanger, instead.  We start by rea
 19 [   2] code: # These paths are specific to this demo only, where root is the root of the course materials repo. | # You may
 20 [   3] code: ETV6_RUNX1_1
 21 [None] mark: This object stores several pieces of information, including:  |  | - A matrix of raw counts, with samples as r
 22 [None] mark: <a id="gene_annotation"></a> | ## 3. Gene Annotation
 23 [   4] code: ETV6_RUNX1_1.var
 24 [None] mark: We can see that our gene ids use ENSEMBL identifiers.  | For each of these identifiers, we would like to know 
 25 [   5] code: # connect to the Human genes database (GRCh38.p14) | h38_mart = bm.Dataset(name="hsapiens_gene_ensembl", |    
 26 [   6] code: # retrieve gene information | h38_genes = h38_mart.query(attributes=["ensembl_gene_id", "external_gene_name"],
 27 [   7] code: h38_genes
 28 [   8] code: # rename the columns | h38_genes = h38_genes.rename(columns={"Gene stable ID": "gene_ids", "Gene name": "gene_
 29 [   9] code: h38_genes.to_csv("hg38_genes_dataframe.csv")
 30 [None] mark: This returns a Pandas DataFrame object:
 31 [  10] code: h38_genes
 32 [None] mark: We can now merge this DataFrame with the DataFrame from our AnnData metadata:
 33 [  11] code: gene_annot = ( |   ETV6_RUNX1_1.var |   .merge(h38_genes, how="left", on="gene_ids") |   .set_axis(ETV6_RUNX1_
 34 [  12] code: gene_annot.head()
 35 [None] mark: Finally, we re-annotate our AnnData, being careful to ensure the order of our genes is the same as in the orig
 36 [  13] code: ETV6_RUNX1_1.var = gene_annot.loc[ETV6_RUNX1_1.var_names]
 37 [None] mark: We will keep only genes in the autosomes, X, Y and MT chromosome (i.e. remove genes in unassembled scaffolds):
 38 [None] raw: vars_to_keep = ( |   ETV6_RUNX1_1.var["chrom"] |   .isin([str(i) for i in range(1, 23)] + ["X", "Y", "MT"]) | 
 39 [None] mark: Note the use of the `.copy()` method. This ensures that we make a new copy of the object (which we replace bac
 40 [None] mark: Since we couldn't directly get a dataframe with chromosoms, we need to use this command to obtain MT gene ids 
 41 [  14] code: # retrieve MT genes | h38_MT = h38_mart.query(attributes=["ensembl_gene_id", "external_gene_name"], |         
 42 [  15] code: h38_MT.head()
 43 [  16] code: h38_MT.to_csv("hg38_MT_genes_dataframe.csv")
 44 [None] raw: # code for if we do obtain the chromosomes in the dataframe  | ETV6_RUNX1_1.var["mt"] = ETV6_RUNX1_1.var["chro
 45 [  17] code: ETV6_RUNX1_1.var["mt"] = False
 46 [  18] code: ETV6_RUNX1_1.var.loc[ETV6_RUNX1_1.var.gene_ids.isin(h38_MT['Gene stable ID']), "mt"] = True
 47 [None] mark: <a id="filtering"></a> | ## 4. Filtering
 48 [None] mark: We start by doing some exploratory analysis of our raw count data, namely in terms of:   |  | - number of tota
 49 [  19] code: sc.pp.calculate_qc_metrics( |     ETV6_RUNX1_1,  |     qc_vars=["mt"],  |     inplace=True,  |     percent_top
 50 [None] mark: This function added several metrics for each barcode, i.e. our observations:
 51 [  20] code: ETV6_RUNX1_1.obs
 52 [None] mark: And also to our genes, i.e. variables:
 53 [  21] code: ETV6_RUNX1_1.var
 54 [None] mark: <a id="filtering_barcodes"></a> | ### 4.1 Filtering barcodes
 55 [None] mark: Since `ETV6_RUNX1_1.obs` is a regular DataFrame, we can use standard plotting libraries to visualise these sta
 56 [  22] code: sns.displot(ETV6_RUNX1_1.obs, x="total_counts", bins=100) | sns.displot(ETV6_RUNX1_1.obs, x="pct_counts_mt", b
 57 [None] mark: Alternatively, we can use `scanpy`'s own plotting functions (histogram is not available, but we can do violin 
 58 [  23] code: sc.pl.violin(ETV6_RUNX1_1, "total_counts") | sc.pl.violin(ETV6_RUNX1_1, "pct_counts_mt") | sc.pl.scatter(ETV6_
 59 [None] mark: We can even do several violin plots at once:
 60 [  24] code: sc.pl.violin( |   ETV6_RUNX1_1, |   ["n_genes_by_counts", "total_counts", "pct_counts_mt"], |   multi_panel=Tr
 61 [None] mark: We can filter our object based on hard thresholds set manually.  | Alternatively, we can define a function tha
 62 [  25] code: def is_outlier(adata, metric: str, nmads: int): |   M = adata.obs[metric] |    |   outlier = (M < np.median(M)
 63 [None] mark: The function returns `True` or `False` depending on whether the value exceeds the specified value of median ab
 64 [  26] code: # create barcode metadata column indicating counts outliers | ETV6_RUNX1_1.obs["counts_outlier"] = is_outlier(
 65 [None] mark: Note that we used the log-transformed counts, as its distribution is less skewed and therefore more suitable f
 66 [  27] code: ETV6_RUNX1_1.obs["genes_outlier"] = is_outlier(ETV6_RUNX1_1, "log1p_n_genes_by_counts", 5) | sns.displot(ETV6_
 67 [  28] code: ETV6_RUNX1_1.obs["topgenes_outlier"] = is_outlier(ETV6_RUNX1_1, "pct_counts_in_top_20_genes", 5) | sns.displot
 68 [None] mark: We also check for outliers with regards to percentage of mitochondrial counts, where we use more strict filter
 69 [  29] code: ETV6_RUNX1_1.obs["mito_outlier"] = is_outlier(ETV6_RUNX1_1, "pct_counts_mt", 3) | (ETV6_RUNX1_1.obs["pct_count
Bash Print cells 64-79 source
IN
cd "e:/RNA-seq-courses/scRNA-seq_python_course_2024" && python -c "
import json
nb=json.load(open('Demonstrations/02_iris_qc_normalisation_yuyao_mod_Feb2025_demo.ipynb',encoding='utf8'))
for i in range(64,80):
    c=nb['cells'][i]
    print('==== idx',i,'[',c.get('execution_count'),']',c['cell_type'])
    print(''.join(c['source'])[:600])
"

OUT
==== idx 64 [ 26 ] code
# create barcode metadata column indicating counts outliers
ETV6_RUNX1_1.obs["counts_outlier"] = is_outlier(ETV6_RUNX1_1, "log1p_total_counts", 5)

# visualise
sns.displot(ETV6_RUNX1_1.obs, x="log1p_total_counts", hue="counts_outlier", multiple="stack")
==== idx 65 [ None ] markdown
Note that we used the log-transformed counts, as its distribution is less skewed and therefore more suitable for the MAD-based filtering we are doing.  

We can repeat this for number of detected genes (also log-transformed) and the percentage of counts in the top 20 genes:
==== idx 66 [ 27 ] code
ETV6_RUNX1_1.obs["genes_outlier"] = is_outlier(ETV6_RUNX1_1, "log1p_n_genes_by_counts", 5)
sns.displot(ETV6_RUNX1_1.obs, x="log1p_n_genes_by_counts", hue="genes_outlier", multiple="stack")
==== idx 67 [ 28 ] code
ETV6_RUNX1_1.obs["topgenes_outlier"] = is_outlier(ETV6_RUNX1_1, "pct_counts_in_top_20_genes", 5)
sns.displot(ETV6_RUNX1_1.obs, x="pct_counts_in_top_20_genes", hue="topgenes_outlier", multiple="stack")
==== idx 68 [ None ] markdown
We also check for outliers with regards to percentage of mitochondrial counts, where we use more strict filters: 
==== idx 69 [ 29 ] code
ETV6_RUNX1_1.obs["mito_outlier"] = is_outlier(ETV6_RUNX1_1, "pct_counts_mt", 3) | (ETV6_RUNX1_1.obs["pct_counts_mt"] > 8)
sns.displot(ETV6_RUNX1_1.obs, x="pct_counts_mt", hue="mito_outlier", multiple="stack")
==== idx 70 [ None ] markdown
Finally, we create a variable which is the union of these conditions, i.e. if the barcode is determined an outlier of _any_ of our filters, then we consider it to be an outlier: 
==== idx 71 [ 30 ] code
ETV6_RUNX1_1.obs["outlier"] = ETV6_RUNX1_1.obs["genes_outlier"] | ETV6_RUNX1_1.obs["genes_outlier"] | ETV6_RUNX1_1.obs["topgenes_outlier"] | ETV6_RUNX1_1.obs["mito_outlier"]

ETV6_RUNX1_1.obs["outlier"].value_counts()
==== idx 72 [ None ] markdown
We can visualise our scatterplot of counts vs detected genes to see which barcodes will be removed: 
==== idx 73 [ 31 ] code
sns.scatterplot(ETV6_RUNX1_1.obs, 
                x = "total_counts", 
                y = "n_genes_by_counts",
                hue = "outlier")
==== idx 74 [ 32 ] code
ETV6_RUNX1_1 = ETV6_RUNX1_1[~ETV6_RUNX1_1.obs["outlier"], :].copy()
==== idx 75 [ None ] markdown
The barcodes we are left with after this filtering we will now consider to be cells. 
==== idx 76 [ None ] markdown
<a id="filtering_genes"></a>
### 4.2 Filtering Genes
==== idx 77 [ None ] markdown
In the same way that we explored several metrics for barcodes, we can also explore them for genes. However, as we will see, downstream analysis can focus on variable genes and will mostly ignore genes for which there is very little data. Therefore, we don't perform as strict filtering on genes as we do on barcodes.    

Still, it is useful to remove undetected genes, i.e. those with zero total counts:
==== idx 78 [ 33 ] code
# count of genes with zero counts
ETV6_RUNX1_1.var["total_counts"].eq(0).value_counts()
==== idx 79 [ None ] markdown
We can use the `sc.pp.filter_genes()` function to do this:
Grep "obs\[\\"outlier\\"\] = " (in e:\RNA-seq-courses\scRNA-seq_python_course_2024\Demonstrations)
4 lines of output
Cell [25] itself draws nothing — it only defines is_outlier. But it's the engine for the graphs that follow, cells [26]–[29] and [31]. That block is what I'll walk through, since it's one continuous piece of reasoning.

The four flagging graphs
Each one follows the identical two-line pattern: flag, then look.


adata.obs["<name>_outlier"] = is_outlier(adata, "<metric>", n)   # decide
sns.displot(adata.obs, x="<metric>", hue="<name>_outlier", multiple="stack")   # check
Cell	Metric	nmads	Catches
[26]	log1p_total_counts	5	droplets with too little or too much RNA
[27]	log1p_n_genes_by_counts	5	too few genes (debris) or too many (doublets)
[28]	pct_counts_in_top_20_genes	5	cells dominated by a handful of genes (low complexity)
[29]	pct_counts_mt	3 + >8%	dying cells
Each graph is a histogram of one metric, with bars stacked in two colours — blue for kept, orange for flagged. multiple="stack" piles them rather than overlapping, so total bar height stays the true cell count.

What each contributes: it shows you where the cut landed. You should see a large blue body and a thin orange fringe out at the tail. If a graph shows orange chewing into the main peak, your threshold is too aggressive. If a visibly separate junk population is still blue, it's too loose. This is the check that stops you deleting real cells on autopilot.

Why log1p_ on the first two
Raw counts are heavily right-skewed — most cells at a few thousand UMIs, a long tail stretching to tens of thousands. MAD assumes a roughly symmetric spread, so on raw counts the flagging would be lopsided.


raw counts:    ▁▇█▅▂▁▁▁▁_______________  ← skewed, long right tail
log1p:         ▁▂▅█▇▅▂▁                   ← roughly symmetric, MAD behaves
log1p = log(x + 1). The +1 exists so that a cell with 0 counts gives 0 rather than log(0) = negative infinity.

Note pct_counts_mt and pct_counts_in_top_20_genes are not logged — they're already percentages bounded 0–100, so they don't have that skew problem.

The combining step [30]

adata.obs["outlier"] = genes_outlier | genes_outlier | topgenes_outlier | mito_outlier
adata.obs["outlier"].value_counts()
| = OR, so a cell is condemned if any single test flagged it. value_counts() prints the tally — how many True vs False. That's your first hard number: how many cells am I about to lose?

⚠️ There's a genuine bug in this line in your notebook. genes_outlier appears twice and counts_outlier is missing — so cell [26]'s work is computed, plotted, then silently thrown away. The other notebook in the same folder, 02_iris_qc_normalisation_demo.ipynb:711, has the correct version (that was commit a4bc3e4b, "Correct omission of counts_outlier"). The fix here is one word:


ETV6_RUNX1_1.obs["outlier"] = (
    ETV6_RUNX1_1.obs["counts_outlier"]     # ← was genes_outlier
    | ETV6_RUNX1_1.obs["genes_outlier"]
    | ETV6_RUNX1_1.obs["topgenes_outlier"]
    | ETV6_RUNX1_1.obs["mito_outlier"]
)
The verdict graph [31]

sns.scatterplot(adata.obs, x="total_counts", y="n_genes_by_counts", hue="outlier")
This is the payoff. Every dot is a cell; colour shows the final decision.


   n_genes │        ░░░░░░░
            │      ░░░░░░░░░       ░ = keep
            │    ░░░░░░░░           ● = discard
            │  ●●░░░░
            │ ●●●
            └──────────────── total_counts
             ↑ low counts, low genes = junk corner
Healthy cells form a tight diagonal band — more RNA means more genes detected. The flagged cells should sit off that band: bottom-left (empty droplets), or far right at low genes (doublets and debris).

What it contributes: it's the sanity check across all metrics at once. The four histograms each show one dimension in isolation; this shows whether the combined verdict carved out a sensible region. If orange dots are scattered randomly through the middle of the healthy band, something's wrong.

Then [32] does the deleting

ETV6_RUNX1_1 = ETV6_RUNX1_1[~ETV6_RUNX1_1.obs["outlier"], :].copy()
~ means NOT — keep every cell where outlier is False. This is the first destructive step; everything before it was labelling and looking.

The block's logic in one line
Define a rule [25] → apply it to four metrics and eyeball each [26–29] → combine into one verdict and count [30] → confirm on a combined plot [31] → only then delete [32].

Want me to apply the counts_outlier fix to the notebook?

e:\RNA-seq-courses\scRNA-seq_python_course_2024\Demonstrations\02_iris_qc_normalisation_yuyao_mod_Feb2025_demo.ipynb
the log part is getting very complicated e3xplain in simpler terms and aalso what is MAD?
Let me do them separately, starting from scratch.

Part 1: What log does
Log squishes big numbers. That's the whole idea.

Watch what happens to a spread of counts:


before log:   100    1,000    10,000    100,000
after log:      2        3         4          5
Each time the number gets 10× bigger, the log goes up by just 1. Huge numbers get pulled in close; small numbers barely move.

Why that helps here
Cell counts are lopsided. Most cells sit around 2,000 UMIs, but a few reach 50,000. Drawn as a histogram:


BEFORE log — squashed left, long tail right

  │█
  │██
  │███
  │████▄▄▂▁▁___________________
  └────────────────────────────
   0                      50,000
Almost everything is crushed into the left edge, and one enormous tail stretches right.


AFTER log — nicely balanced

  │      ▄█▄
  │    ▄███▄
  │  ▄█████▄
  │▁▄███████▄▁
  └────────────
Now it's a sensible hump with roughly equal amounts on each side. Same data — just re-spaced so the picture is readable.

Why the flagging needs that
The outlier rule draws its cutoffs the same distance above and below the middle. That only makes sense if the data is roughly symmetric.

On skewed data it goes wrong:


skewed:  the lower cutoff falls below zero (meaningless — counts can't be negative),
         while the upper cutoff sits so far out it catches almost nothing
Log makes the data symmetric, so equal-sized cutoffs on both sides are fair.

What log1p means
Just log(x + 1).

The +1 is a safety patch. Some cells have 0 counts, and log(0) is undefined — it breaks. Adding 1 first means 0 becomes log(1) = 0. Clean.


log1p(0)   = 0      ← safe
log1p(9)   = 2.3
log1p(999) = 6.9
The 1p literally reads "1 plus."

Which metrics get logged
Metric	Logged?	Why
total_counts	✅ yes	ranges 0 → 50,000, badly skewed
n_genes_by_counts	✅ yes	same problem
pct_counts_mt	❌ no	already a % (0–100), not skewed
pct_counts_in_top_20_genes	❌ no	same
Percentages are already on a tame scale, so they don't need squishing.

Part 2: MAD
MAD = the typical distance from the middle.

Worked through with tiny numbers
Seven cells, mitochondrial %:


4, 5, 6, 6, 7, 8, 95
Six normal cells, one dying cell.

Step 1 — find the middle value (the median).


4, 5, 6, [6], 7, 8, 95
             ↑  middle = 6
Step 2 — how far is each one from 6? (ignore minus signs)


value:      4   5   6   6   7   8   95
distance:   2   1   0   0   1   2   89
Step 3 — find the middle of those distances.


0, 0, 1, [1], 2, 2, 89
           ↑  MAD = 1
MAD = 1. In words: a typical cell sits about 1 unit from the middle.

The name explained
Median Absolute Deviation, read backwards:

Deviation = how far from the middle
Absolute = ignore the minus sign
Median = take the middle one of those
You use the median twice — once to find the centre, once to find the typical distance from it. That's the trick to remembering it.

Why not the ordinary average?
Because one extreme value wrecks an average.


average distance = (2+1+0+0+1+2+89) ÷ 7 = 13.6   ← the 95 dragged it up
MAD              =                        1      ← unmoved
You're using this number to catch the weird cell. If the weird cell inflates your yardstick, it makes itself look normal and escapes. MAD refuses to be pushed around, which is why it's called robust.

And nmads=5 means

median = 6,  MAD = 1,  5 MADs = 5

lower cutoff = 6 − 5 = 1
upper cutoff = 6 + 5 = 11
Anything outside 1 to 11 gets flagged. The cell at 95 is caught. ✓

How the two fit together

raw counts  ──log──►  symmetric  ──MAD──►  fair cutoffs
 (skewed)              (balanced)          both sides
Log makes the data a sensible shape. MAD then measures spread in a way that outliers can't corrupt. Using MAD on unlogged counts would still technically run, but the cutoffs would be lopsided and mostly useless.
