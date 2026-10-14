Interview Prep: Differential Expression (Notebook 06)
Question this notebook answers: Within the same cell type, which genes are expressed differently in leukaemia (ETV6-RUNX1) compared with healthy bone marrow (PBMMC)?

Analogies alternate between Harry Potter (HP) and Avatar: The Last Airbender (ATLA). ⭐ = non-negotiable.

🔥 READ THIS FIRST: the 10 things to remember
⭐ Pseudobulk. Sum the raw counts of all cells of one cell type in one sample, so you get one "bulk-like" profile per sample × cell type. Then use bulk RNA-seq tools (edgeR, DESeq2, limma).
⭐ Why pseudobulk? Cells from the same person are not independent replicates. The replicate is the sample (patient), not the cell. Testing cells individually (e.g. Wilcoxon on 5,000 cells) gives thousands of false positives. This is called pseudoreplication (Squair et al., Nat Commun 2021).
⭐ Use raw counts, not log-normalised values. edgeR models counts with a negative binomial distribution and does its own normalisation (TMM).
⭐ Drop pseudobulk samples with too few cells (here min_cells=20). Profiles built from very few cells are noisy.
⭐ Filter lowly expressed genes (filterByExpr): 31,801 → 14,543 genes. Lowly expressed genes have no statistical power and only add to the multiple-testing burden.
⭐ Design formula ~SampleGroup tells the model what explains expression. Add other factors if needed (e.g. ~batch + SampleGroup, or ~patient + condition for paired samples).
⭐ edgeR quasi-likelihood F-test (glmQLFit / glmQLFTest). It borrows information across genes (empirical Bayes) because with only 3–4 replicates each gene's variability can't be estimated reliably on its own. robust=True protects against outlier genes.
⭐ Use the adjusted p-value (FDR, Benjamini–Hochberg), not the raw p-value. With 14,543 tests, about 727 genes would reach p < 0.05 by chance alone.
⭐ Know the direction of the fold change. baseline="ETV6-RUNX1", group_to_compare="PBMMC" means positive log_fc = higher in PBMMC (healthy). See the ⚠️ warning in Step 9.
⭐ DE ≠ marker genes. Marker genes compare clusters within the data, which is exploratory. DE compares conditions across biological replicates, which is real statistics.
🗺️ The notebook in one picture
30,461 cells (7 samples: 4 ETV6-RUNX1 + 3 PBMMC), already clustered & labelled
   │
   ▼  Step 3: Pseudobulk: sum counts per (sample × cell type), drop if < 20 cells
71 pseudobulk "samples" × 31,801 genes
   │
   ▼  Step 4: Pick one cell type: B (c1) → 7 pseudobulk samples (4 vs 3)
   ▼  Step 5: filterByExpr → 14,543 genes
   ▼  Step 6: edgeR model ~SampleGroup, fit (NB GLM + quasi-likelihood)
   ▼  Step 7: QC: PCA of the 7 pseudobulk samples
   ▼  Step 8: Test contrast PBMMC vs ETV6-RUNX1 → results table
   ▼  Step 9: Volcano plot
   ▼  Step 10: Gene set enrichment (ORA with MSigDB Hallmarks)
   ▼  Steps 11–13: Loop over ALL cell types, count up/down genes, plot one gene everywhere

Step 1: Load the libraries
import decoupler as dc      # pseudobulk + gene set analysis
import pertpy as pt         # Python wrapper around edgeR (R)
import scanpy as sc
import numpy, pandas, matplotlib, seaborn

Library	Its job here
decoupler (Saez-Rodriguez lab)	get_pseudobulk, MSigDB gene sets, enrichment (ORA)
pertpy (Theis lab)	Runs edgeR from Python (pt.tl.EdgeR)
rpy2 (used inside a helper function)	Lets Python call R directly
seaborn	Nicer statistical plots
📌 Point to make: edgeR and DESeq2 are R packages. Python reaches them through rpy2, or uses native Python alternatives such as PyDESeq2.

Step 2: Load the clustered data and look at it
adata = sc.read_h5ad('.../Caron_clustered.PBMMCandETV6RUNX1.h5ad')
sc.pl.embedding(adata, "X_umap_corrected", color=["label"], legend_loc="on data")
adata.obs[["label", "SampleName"]].value_counts()

What's in the object:

30,461 cells × 33,102 genes: the full dataset (not the 500-cells-per-sample subset), only 2 groups: ETV6-RUNX1 (4 patients) + PBMMC (3 healthy donors).
obs["label"]: cell-type labels from clustering + annotation, e.g. B (c1), T (c4), Mono (c10), Erythrocytes (c5).
obsm["X_umap_corrected"]: UMAP from batch-corrected data (from the integration step).
value_counts() counts cells per cell type × sample:

B (c1) in ETV6-RUNX1_2 = 4,939 cells ← the leukaemic B-cell blasts
Some combinations have just 1 cell, which is far too few to say anything.
📌 Expert insight: B (c1) is a huge cluster in the leukaemia samples (the cancer cells) but contains only 40–187 cells in the healthy samples. So the "same cell type" comparison is really malignant blasts vs a small set of normal B cells. Say this, it shows you think about the biology.

Step 3: Create the pseudobulk ⭐⭐⭐ (the most important concept)
pdata = dc.get_pseudobulk(adata, sample_col="SampleName", groups_col="label",
                          mode="sum", min_cells=20)

🧠 Simple explanation
For every sample × cell type combination, add up the raw counts of all its cells for each gene. Each combination becomes one row, like a mini bulk RNA-seq sample.

Argument	Meaning
sample_col="SampleName"	the biological replicate (patient/donor)
groups_col="label"	the cell type
mode="sum"	sum counts (keeps them as counts, which edgeR needs)
min_cells=20	throw away combinations with fewer than 20 cells
Result: pdata = 71 pseudobulk samples × 31,801 genes (genes with no counts at all were dropped). New obs columns: psbulk_n_cells (how many cells were summed), psbulk_counts (total counts). New layer psbulk_props: the fraction of cells expressing each gene. pdata.X shows whole numbers (510, 249, 2...): still raw counts ✅

🪄 Analogy (HP)
You want to know whether Hogwarts students are better at Charms than Durmstrang students.

❌ Wrong: test 1,000 Hogwarts students vs 1,000 Durmstrang students as 2,000 independent data points. They share the same school, teachers and castle, so they're not independent. You'd "prove" a difference that's just about one school's quirks.
✅ Right: each school-year is one data point. Summarise each school (sample), then compare schools. You need several schools per side to draw a conclusion.
That's pseudobulk: the patient is the unit of replication, not the cell.

❓ Expert questions
Q: Why not just run a Wilcoxon test on single cells between conditions? A: Cells from the same donor are correlated. Treating them as independent inflates the sample size from 7 to thousands, so p-values become tiny and you get many false positives. Squair et al. 2021 showed pseudobulk methods are much more accurate. Single-cell tests effectively test "these cells differ", not "these conditions differ across individuals".

Q: Sum or mean? A: Sum, so the data stay as counts that a negative binomial model can handle. Differences in library size are taken care of by edgeR's normalisation.

Q: Why min_cells=20? A: A pseudobulk built from 3 cells is very noisy and dominated by chance. A minimum (commonly 10–30) keeps profiles reliable. The trade-off is that you may lose replicates for rare cell types.

Q: What are alternatives to pseudobulk? A: Mixed models with a random effect for the donor (e.g. MAST with random effects, NEBULA, muscat's mixed models). They're more complex and slower; pseudobulk usually performs as well or better.

Q: What if a cell type exists in only one condition? A: You can't do DE for it, because there's nothing to compare. That's a differential abundance question (notebook 07).

Step 4: Pick one cell type: B (c1)
b_c1_pdata = pdata[pdata.obs['label'] == "B (c1)"].copy()

Keep only rows where the label is B (c1) → 7 pseudobulk samples (4 ETV6-RUNX1 + 3 PBMMC).
.copy() makes an independent object (not a "view" that points back to pdata).
Look at psbulk_n_cells: ETV6-RUNX1 = 1,761–4,939 cells each, PBMMC = 40, 69, 187 cells. The healthy pseudobulks are built from far fewer cells, so they're noisier.

Step 5: Filter lowly expressed genes (filterByExpr) ⭐
def filter_by_expression(adata, group):
    ... edger.filterByExpr(adata.X.T, ro.FactorVector(adata.obs[group]))
b_c1_pdata = filter_by_expression(b_c1_pdata, 'SampleGroup')

🧠 Simple explanation
The helper function calls edgeR's filterByExpr in R (via rpy2):

adata.X.T → transpose (.T), because R/edgeR wants genes × samples, while AnnData is samples × genes.
FactorVector(adata.obs[group]) → tells edgeR which group each sample belongs to.
filterByExpr keeps a gene only if it has enough counts (roughly ≥10 counts, converted to CPM) in at least as many samples as the smallest group (here 3).
Returns True/False per gene → keep the True ones.
Result: 31,801 → 14,543 genes.

🪄 Analogy (ATLA)
Before a council meeting in Ba Sing Se, the Dai Li remove everyone who only whispers. You can't judge whether someone's opinion changed if you can barely hear them. Removing whisperers also means fewer people to cross-check, so less multiple-testing penalty.

❓ Expert questions
Q: Why filter genes before DE? A: (1) Genes with very low counts can't reach significance anyway. (2) Fewer tests means a smaller multiple-testing correction, so more power for the rest. (3) They distort the mean–variance (dispersion) trend that edgeR estimates.

Q: Why filter based on the smallest group size? A: So a gene expressed in only one group (e.g. switched on just in leukaemia) is still kept. That's exactly the kind of gene you want to find.

Step 6: Build and fit the edgeR model ⭐
b_c1_edgr = pt.tl.EdgeR(b_c1_pdata, design="~SampleGroup")
b_c1_edgr.fit(robust=True)

🧠 Simple explanation
Design ~SampleGroup: "expression of each gene depends on whether the sample is ETV6-RUNX1 or PBMMC." The ~ means "is modelled by".

What fit() does inside edgeR, in order:

TMM normalisation (calcNormFactors): corrects for library size and composition bias (when a few very highly expressed genes take up a big share of reads in one sample, all the other genes look lower there).
Estimate dispersion (estimateDisp): how much each gene varies between replicates, beyond Poisson noise. It shrinks each gene's estimate toward the trend across all genes (empirical Bayes).
Quasi-likelihood GLM fit (glmQLFit): fits a negative binomial generalised linear model for each gene, with an extra layer of gene-specific variability. robust=True means outlier genes don't distort the shared estimate.
🪄 Analogy (ATLA)
With only 3–4 samples per group, judging each gene's variability on its own is like judging a bender's skill from 3 moves: unreliable. The Order of the White Lotus approach: each member's judgement is combined with the collective wisdom of the whole order. Each gene's variability estimate is pulled toward what similar genes show. That's empirical Bayes shrinkage.

❓ Expert questions
Q: Why a negative binomial distribution? A: RNA-seq counts are overdispersed: the variance is larger than the mean, because biological variability adds to sampling (Poisson) noise. The negative binomial has an extra dispersion parameter to capture that.

Q: What is TMM normalisation? A: Trimmed Mean of M-values. It assumes most genes are not DE. It compares each sample to a reference, trims extreme genes, and uses the weighted mean of the remaining log-ratios as a scaling factor. That corrects for composition bias, which simple CPM doesn't.

Q: Why the quasi-likelihood F-test instead of the likelihood ratio test? A: QL accounts for uncertainty in the dispersion estimate, so it controls false positives better, especially with few replicates. The F column in the results shows this is the QL F-test.

Q: How would you handle a batch effect or paired samples in DE? A: Put it in the design: ~batch + SampleGroup or ~patient + condition. Don't use batch-corrected expression values for DE. Correct in the model instead, using raw counts.

Q: edgeR vs DESeq2 vs limma-voom? A: All are well-validated bulk methods that work well on pseudobulk. edgeR and DESeq2 use a negative binomial GLM with shrinkage; limma-voom transforms to log-CPM with precision weights and uses linear models. Benchmarks (Squair 2021, Crowell et al. 2020 muscat) show all three perform well on pseudobulk.

Step 7: QC: PCA of the pseudobulk samples
b_c1_pdata.layers['counts'] = b_c1_pdata.X.copy()   # back up raw counts
sc.pp.normalize_total(b_c1_pdata); sc.pp.log1p(b_c1_pdata)
sc.pp.scale(b_c1_pdata, max_value=10)               # each gene mean 0, sd 1, clip at 10
sc.tl.pca(b_c1_pdata)
sc.pl.pca(b_c1_pdata, color='SampleGroup', ..., size=100)

🧠 Simple explanation
Plot the 7 pseudobulk samples on PCA to check do samples group by condition, and are there outliers?

What the plot shows:

The 4 ETV6-RUNX1 samples sit tightly together. The 3 PBMMC samples are separated from them along PC1 but spread out from each other.
→ There's a strong condition effect (good for DE). PBMMC is more variable, which makes sense because those pseudobulks come from few cells (40–187) and different healthy donors.
📌 Points:

Saving layers['counts'] before normalising is good practice. The model was already fitted on raw counts.
The warning ImplicitModificationWarning: ... view appears because filter_by_expression returned a view. Harmless, but adding .copy() avoids it.
Ideally do this QC before fitting the model, to catch outlier samples or sample swaps.
Q: What if one sample sat with the other group in this PCA? A: Check for a sample swap, a failed library, or real biology (e.g. a patient in remission). Also check its cell number and counts. Consider removing it or adding a covariate.

Step 8: Test the contrast → the results table ⭐
b_c1_res_df = b_c1_edgr.test_contrasts(
    b_c1_edgr.contrast(column="SampleGroup", baseline="ETV6-RUNX1", group_to_compare="PBMMC"))

The contrast: "compare PBMMC against the baseline ETV6-RUNX1", so log fold change = log2(PBMMC / ETV6-RUNX1).

Reading the results table (know every column):

Column	Meaning	Example (HTR1F)
log_fc	log2 fold change (group_to_compare vs baseline)	9.93 → about 2^9.9 ≈ 975-fold difference
logCPM	average expression (log2 counts per million): how highly expressed the gene is overall	6.80
F	quasi-likelihood F-statistic: strength of evidence	1005
p_value	raw p-value	3.2e-12
adj_p_value	FDR (Benjamini–Hochberg) adjusted p-value: use this one	4.7e-08
The bottom of the table (RAG1, CLIC4: log_fc ≈ 0, adj p ≈ 1) shows no change genes.

🪄 Analogy (HP): multiple testing
Professor Trelawney makes 14,543 predictions. Even if she has no real powers, about 5% (~727) will "come true" by luck. Benjamini–Hochberg asks: of all the predictions I'm calling true, what fraction are probably luck? An FDR of 0.05 means you accept that about 5% of your significant genes are false discoveries.

❓ Expert questions
Q: What's the difference between p-value and FDR? A: The p-value is the chance of seeing this result if there's no real effect, for one test. FDR controls the expected proportion of false positives among all genes called significant. With thousands of tests you must use FDR.

Q: Should you also filter on fold change? A: Often yes (e.g. |log2FC| > 1) for biological relevance. A better way is to test against a threshold (edgeR glmTreat) rather than filter after the fact.

Q: What does logCPM tell you? A: How highly the gene is expressed overall. Large fold changes for genes with very low logCPM are less reliable.

Step 9: Volcano plot ⭐
b_c1_edgr.plot_volcano(b_c1_res_df, pval_thresh=0.05, log2fc_thresh=0)

🧠 How to read it
x-axis = log2 fold change (right = higher in PBMMC, left = higher in ETV6-RUNX1, according to the contrast).
y-axis = −log10(adjusted p): higher = more significant. The dashed horizontal line = FDR 0.05 (−log10(0.05) ≈ 1.3).
Red "Up" / blue "Down" = significant genes on each side. Grey = not DE.
Top genes are labelled: HTR1F, FXYD2, L3MBTL4 (right), OVCH2, NRN1, IGF2BP1, ECM1 (left).
🪄 Analogy (ATLA)
A volcano plot is like the Fire Nation's wanted posters: the further someone is from the centre (x = how much they changed) and the higher up (y = how sure we are), the more "wanted" they are. The top corners hold the most important genes.

⚠️ CRITICAL: the direction contradiction in this notebook
The code sets baseline="ETV6-RUNX1", group_to_compare="PBMMC". In pertpy, the contrast is group_to_compare minus baseline, so positive log_fc (red, "Up") = higher in PBMMC (healthy). The design also agrees: by default the reference level is the first alphabetically, and "ETV6-RUNX1" comes before "PBMMC".
But the notebook's text later says "HTR1F is highly upregulated… in leukemic cells". That contradicts the code. One of the two is wrong.
Before the interview, check it yourself (run this in the course environment, after Step 7):
cpm = b_c1_pdata[:, "HTR1F"].layers["counts"].flatten() / b_c1_pdata.obs["psbulk_counts"] * 1e6
print(pd.Series(cpm.values, index=b_c1_pdata.obs["SampleName"]))

If the PBMMC rows are higher, the text is wrong and HTR1F is higher in healthy cells.
What to say if asked: "I always check the contrast direction: which group is the reference. Here log fold change is PBMMC relative to ETV6-RUNX1, so positive values mean higher in healthy samples." Experts love this, because it's a very common mistake.
Step 10: Gene set analysis (ORA with MSigDB Hallmarks)
msigdb = dc.get_resource('MSigDB')                     # ~4 million gene–set pairs
msigdb = msigdb[msigdb['collection']=='hallmark']      # keep the 50 Hallmark sets
msigdb = msigdb[~msigdb.duplicated(['geneset', 'genesymbol'])]   # remove duplicates
msigdb['geneset'] = [n.split('HALLMARK_')[1] for n in ...]       # shorten names
enr_pvals = dc.get_ora_df(df=<significant genes with log_fc > 0>, net=msigdb,
                          source='geneset', target='genesymbol')
dc.plot_dotplot(...)                                   # top 15 enriched sets

🧠 Simple explanation
A list of 500 DE genes is hard to interpret. Gene set analysis asks: do my DE genes cluster in known biological pathways?

MSigDB = a big library of gene sets. Hallmark = 50 curated, non-redundant sets (e.g. INTERFERON_GAMMA_RESPONSE, MYC_TARGETS_V1, OXIDATIVE_PHOSPHORYLATION).
ORA (Over-Representation Analysis) = a Fisher's exact / hypergeometric test: are there more of my DE genes in this set than expected by chance?
Input here: genes with adj p < 0.05 and log_fc > 0 (higher in PBMMC, given the contrast direction).
Dot plot: x = combined score, dot size = odds ratio, colour = FDR.
🪄 Analogy (HP)
If 1% of all wizards are Death Eaters, but 30% of the people at a suspicious meeting are Death Eaters, that's enrichment. Fisher's test tells you how unlikely that is by chance.

⚠️ In this notebook, this step is BROKEN
The filtering cell has output text accidentally pasted into the code (...=='hallmark']AnnData object with n_obs × ...), which causes a SyntaxError. So the duplicate-removal line never ran, and ORA failed with ValueError: net contains repeated edges. The dot plot and all the cells after it have no outputs (they never ran).
Fix: delete the pasted text so the line reads msigdb = msigdb[msigdb['collection']=='hallmark'].
💡 If they ask whether you ran it: "The ORA cell failed because of a corrupted line, so duplicates weren't removed. I spotted the cause and how to fix it." That's a good sign of debugging skill.
❓ Expert questions
Q: ORA vs GSEA? A: ORA uses a cut-off list (significant genes only) and ignores how strongly each gene changes. GSEA uses the whole ranked list (e.g. by t-statistic or signed −log p) with no arbitrary cut-off, so it can catch small, coordinated shifts across many genes. decoupler also offers ULM/MLM on the statistics.

Q: Why does the background (universe) matter in ORA? A: The background should be the genes that could have been detected, i.e. the 14,543 tested genes, not the whole genome. Using the whole genome inflates enrichment.

Q: Why run up- and down-regulated genes separately? A: Pathways can be activated in one direction. Mixing them can hide or blur the signal.

Steps 11–13: Repeat for every cell type and summarise
def pseudo_bulk_dge(adata, label, design, group, column, baseline, group_to_compare): ...
de_results = pseudo_bulk_dge(pdata, 'label', '~SampleGroup', 'SampleGroup',
                             'SampleGroup', 'ETV6-RUNX1', 'PBMMC')
is_de = decide_tests_per_label(de_results, 0.05)
is_de.apply(lambda x: x.value_counts()).T

🧠 Simple explanation
pseudo_bulk_dge wraps Steps 4–8 in a loop over every cell type: subset → filterByExpr → EdgeR → fit → test. The results go into a dictionary {cell type: results}.
try / except: if a cell type can't be tested (e.g. it has pseudobulks in only one condition after min_cells=20), print a message and skip it instead of crashing.
@dataclass dgeResults is a small container holding the edgeR object and the results table together, plus a plot_volcano shortcut.
decide_tests_per_label labels each gene as +1 (up), −1 (down), 0 (not significant) at FDR < 0.05. This is like decideTests in limma.
.value_counts() → a table of how many genes go up/down in each cell type. It shows which cell types change most between conditions.
de_results["B (c1)"] and its volcano plot reproduce Steps 8–9 (a sanity check that the loop gives the same answer).
Then: look at one gene across all cell types
pdata.layers["counts"] = pdata.X.copy()
sc.pp.normalize_total(pdata); sc.pp.log1p(pdata)
plot_data["logcounts"] = pdata[:, "HTR1F"].X.flatten().tolist()
sns.catplot(data=plot_data, x="SampleGroup", y="logcounts", col="label", col_wrap=5, ...)

Normalise and log the pseudobulk data (for plotting only, never for edgeR).
Plot HTR1F per pseudobulk sample (each dot = one patient), with one panel per cell type.
This checks whether a change is specific to one cell type or shared across several. The notebook says it's also changed in other B-cell clusters.
📌 Expert point: always look at the expression of your top hits per sample. One extreme sample can create a "significant" result.

🧨 Critical review: issues to raise yourself
Fold-change direction vs the text (Step 9): positive log_fc = higher in PBMMC, which contradicts the notebook's statement about leukaemia.
Cell type confounded with disease: B (c1) in leukaemia = malignant blasts (thousands of cells). In healthy samples it's a few normal B cells (40–187). The DE mixes "cancer vs normal" with "different B-cell maturation states".
Very few replicates (4 vs 3), and the PBMMC pseudobulks are built from few cells → limited power and noisy estimates.
ORA cell broken by pasted output text, so the later cells never ran.
Not checked that adata.X holds raw counts before summing. Always verify (integers, large maximum).
No covariates (e.g. sex, age, processing batch) in the design.
filterByExpr returns a view → harmless warning; add .copy().
❓ Rapid-fire expert questions (cover the answer, practise)
What's the unit of replication in scRNA-seq DE? → The biological sample (patient/donor/mouse), not the cell.
What is pseudoreplication? → Treating non-independent measurements (cells from one donor) as independent replicates, which inflates significance.
When is a single-cell Wilcoxon test OK? → For marker gene discovery between clusters within a dataset (sc.tl.rank_genes_groups). Even then the p-values are inflated (cells were clustered on the same data, "double dipping"), so use them for ranking, not as real statistics.
Why can't you use batch-corrected/integrated values for DE? → Most integration methods correct the embedding, not counts. Corrected values break the count model and can create false signals. Use raw counts and put batch in the design.
What's overdispersion? → Variance larger than the mean (more than Poisson). Handled by the negative binomial dispersion.
What does a log2FC of −2 mean? → 4-fold lower in group_to_compare than in the baseline.
What if a gene is significant but with a tiny fold change? → Probably statistically real but biologically minor. Use a fold-change threshold or glmTreat.
How many replicates do you need? → At least 3 per group; more is much better. Cells don't replace replicates.
DE vs differential abundance (DA)? → DE: does a gene's expression change within a cell type? DA: does the proportion of a cell type change? (Notebook 07: milo, scCODA.)
What is composition bias? → If a few genes dominate the reads in one sample, other genes look artificially lower. TMM (edgeR) and median-of-ratios (DESeq2) correct this.
Paired design? → If you have before/after samples from the same patient, use ~patient + condition to compare within each patient. This gives more power.
What does robust=True do? → Robust empirical Bayes: genes with unusually high or low variability don't distort the shared prior used for all genes.
🎤 The 90-second explanation (say it out loud)
"For differential expression between conditions, I use a pseudobulk approach. Cells from the same patient aren't independent, so the real replicates are the patients. I sum the raw counts per sample and cell type, dropping combinations with fewer than 20 cells. That gives a bulk-like count matrix with one row per sample per cell type.

For each cell type I filter lowly expressed genes with edgeR's filterByExpr, then fit a negative binomial GLM with design ~SampleGroup. edgeR applies TMM normalisation, estimates dispersions with empirical Bayes shrinkage, which matters with only 3–4 replicates, and I test with the quasi-likelihood F-test. I use FDR-adjusted p-values because we test about 14,000 genes.

I check sample-level PCA for outliers, look at volcano plots, and always confirm the contrast direction. Here the baseline is ETV6-RUNX1, so positive fold changes mean higher in healthy samples. Then I interpret the gene lists with gene set analysis, either ORA on Hallmark sets or, preferably, GSEA on the full ranked list. Finally I repeat across all cell types to see which ones change most."
