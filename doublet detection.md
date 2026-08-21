What a doublet actually is
In 10x Genomics droplet-based scRNA-seq, cells are suspended in liquid and pushed through a microfluidic chip so that — ideally — each droplet catches one cell and one barcoded bead. That barcode is what makes the whole method work: every mRNA molecule from that cell gets stamped with the same barcode, so afterwards you can say "all these reads came from one cell."

But loading is random. Sometimes two cells end up in the same droplet. Both get stamped with the same barcode, so downstream you have no way to tell them apart — the software sees one barcode, and calls it one cell. That's a doublet (or multiplet for 3+).

This isn't rare. Doublet rate scales with how many cells you load: roughly ~0.8% per 1,000 cells recovered. Load for 10,000 cells and about 8% of your "cells" are actually two cells glued together. That's a big contaminating fraction.

Why it wrecks downstream analysis
This is the part that matters, and it's why the notebook says doublets "can be misclassified and thus confound downstream analysis."

1. They invent cell types that don't exist. Suppose your sample has T cells (expressing CD3E) and monocytes (expressing LYZ). A T-cell + monocyte doublet expresses both markers. When you cluster, that hybrid doesn't sit neatly in either group — and if you have enough of them, they form their own cluster. A naive interpretation: "I've discovered a novel CD3E+/LYZ+ intermediate population!" You haven't. You've discovered your pipetting.

This is a genuinely famous failure mode — several published "novel transitional cell states" have turned out to be doublet artefacts.

2. They fake developmental trajectories. If you run pseudotime/trajectory analysis, doublets sit between two real cell types in expression space and look exactly like cells transitioning from one to the other. You get a smooth, convincing-looking differentiation path that is pure artefact.

3. They corrupt differential expression. A doublet contributes monocyte genes to your T-cell cluster's average. Your "T cell markers" get polluted with genes T cells don't express.

4. In your specific dataset it's worse. This is ETV6_RUNX1 — a childhood B-cell acute lymphoblastic leukaemia sample. The whole biological question is about identifying malignant blast populations and how they differ from normal cells. A leukaemic-cell + normal-cell doublet is exactly the kind of thing that could be mistaken for a rare aberrant subclone. High stakes for a false positive.

Why QC filtering doesn't already catch them
Reasonable question — section 4.1 already removed outliers on log1p_total_counts and log1p_n_genes_by_counts, and a doublet has roughly double the RNA. Doesn't that filter handle it?

Only partly. Two reasons it isn't enough:

Cells naturally vary hugely in size and RNA content. A large activated cell legitimately has more RNA than a small resting lymphocyte. A doublet of two small cells can have a perfectly normal total count. There's no clean count threshold that separates doublets from big cells.
Homotypic vs heterotypic. Two T cells in one droplet (homotypic) look like... a T cell with a bit more RNA. Mostly harmless and nearly undetectable. The dangerous ones are heterotypic doublets — two different cell types — and those are detectable by their expression profile, not their count total.
That's the gap Scrublet fills: it works on the transcriptome pattern, not just the totals.

How Scrublet does it (the clever bit)
sc.pp.scrublet(ETV6_RUNX1_1) in cell 86 runs this:

Build fake doublets. Take your real cells, pick random pairs, and add their count vectors together. That's a simulated doublet — and crucially, you know it's a doublet, because you made it. Do this thousands of times.
Mix simulated doublets into the real data and project everything into shared expression space (PCA + a nearest-neighbour graph).
Score each real cell by asking: of my nearest neighbours, what fraction are simulated doublets? A real singlet sits among other singlets → low score. A real doublet sits in the same region as the simulated ones → high score.
Threshold the score distribution (Scrublet picks the cut automatically, at the dip between the two humps) to call each barcode True or False.
It adds two columns to .obs:

doublet_score — the continuous fraction (0 to 1)
predicted_doublet — the boolean call after thresholding
Note it needs no external reference and no marker genes — it generates its own ground truth from your data. That's why the method is popular.

The filtering step

ETV6_RUNX1_1 = ETV6_RUNX1_1[~ETV6_RUNX1_1.obs["predicted_doublet"], :].copy()
ETV6_RUNX1_1.obs["predicted_doublet"] — a column of True/False, one per cell
~ — the NOT operator, flipping every True→False and vice versa. So this becomes "cells that are not doublets"
[ rows , : ] — AnnData indexing is [cells, genes]; the bare : means "keep all genes"
.copy() — the slice above is only a view (a window onto the original object). .copy() makes it a real standalone object. Without it, scanpy will throw ImplicitModificationWarning when you later try to write to it.
One thing worth noticing about the ordering
Cell 92 mentions there's a choice here: filter doublets now, or wait until after clustering and drop whole clusters with high average doublet scores. The notebook takes the simple route for teaching purposes.

The post-clustering approach is often considered more robust, because a single cell's score is noisy, but a whole cluster with a uniformly elevated score is strong evidence. It also lets you inspect what you're removing rather than trusting a threshold blind. Worth remembering when you run this on your own data.

Also note Scrublet is run per sample, not on merged data — doublets form inside one droplet emulsion, so simulating them across samples would be meaningless. In the sc_preprocess() function at cell 103, sc.pp.scrublet(adata) sits inside the per-sample function, which is correct.
