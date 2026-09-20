Written for: you — interview prep, not for the notebook.

Here's the whole notebook broken down, with PCA and UMAP as the main event.

---

# The big picture

Your data is a table: **5,500 cells × 29,346 genes**. That's 29,346 numbers describing each cell — a point in 29,346-dimensional space. You cannot plot that, cluster it reliably, or even measure distances in it sensibly.

The entire notebook does one thing: **shrink 29,346 numbers per cell down to ~50, and then to 2 for a picture** — while throwing away noise, not biology.

The chain is:

```
normalise → log → pick 2,934 variable genes → PCA (→50 dims) → neighbour graph → UMAP (→2 dims)
```

Each arrow discards information on purpose. The skill is discarding the *right* information.

---

# Step 1 — Normalise and log (cell 12)

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

**Why it's here:** cell 8 checked `adata.X.sum(axis=1)` and got 4557, 7745, 7181… — different totals per cell. Cell 9 showed a max of 8900. Both say: *this is raw counts*. The R→Python conversion dropped the normalised layer, so they redo it.

- `normalize_total` → every cell is rescaled to the same total (10,000 counts). Removes **sequencing depth** differences. Without it, a deeply-sequenced cell looks "different" from a shallow one for purely technical reasons, and PC1 would just be library size.
- `log1p` → `log(x+1)`. Gene expression is wildly skewed; a handful of genes have counts in the thousands while most are 0–5. Log compresses that, so a change from 10→20 counts counts the same as 100→200 (a 2-fold change either way). It also stops a few loud genes from dominating the variance.

> **Order matters and interviewers check this:** HVG selection with `flavor="seurat"` *expects* log-normalised data. `flavor="seurat_v3"` expects raw counts. Getting this backwards silently gives you garbage genes.

---

# Step 2 — Feature selection: highly variable genes (cells 13–16)

```python
sc.pp.highly_variable_genes(adata, n_top_genes=int(0.1 * adata.n_vars), flavor="seurat")
# → 2,934 genes kept out of 29,346
```

**The idea in one line:** a gene that's expressed at the same level in every cell tells you nothing about how cells differ. Only genes that *vary* carry the signal that separates cell types.

Most of the 29,346 genes are either always-off (dropouts, unexpressed) or always-on at a flat level (housekeeping genes like ribosomal/actin). Keeping them adds noise and compute cost, not information.

**The catch — and this is the clever part.** In count data, variance is coupled to mean: highly expressed genes are *automatically* more variable (Poisson noise scales with the mean). If you just ranked genes by raw variance, you'd get a list of the most-expressed genes, not the most-informative ones.

The **`seurat` flavor** fixes this:
1. For each gene compute mean expression and **dispersion** (roughly variance/mean).
2. **Bin genes by their mean expression.**
3. Within each bin, z-score the dispersions (subtract the bin's mean dispersion, divide by its SD).
4. Take the top genes by this *normalised* dispersion.

So you're asking: *"is this gene more variable than other genes expressed at a similar level?"* That's a fair comparison. You end up picking variable genes from across the whole expression range — low, medium and high.

**That's exactly what the plot in cell 14 is for.** `sc.pl.highly_variable_genes` plots dispersion vs mean, with selected genes in black. The sanity check the comment describes: the black points should be **spread along the x-axis**, not all bunched at the high-expression end. If they're all clumped at high means, your mean–variance correction didn't work.

Cell 16's violin plots are the biological confirmation: each HVG is high in *some* cells, near zero in most others. That bimodal shape is exactly what a cell-type marker looks like.

---

# Step 3 — PCA (cells 18–21) ⭐

## What PCA actually is

Imagine a 3D cloud of points shaped like a flat, tilted pancake. It lives in 3D, but really it's 2D — the thickness is negligible. PCA finds that: it identifies the directions in which the data actually spreads out, and tells you how much spread lives in each.

Formally:
- **PC1** is the single direction (axis) through the data along which the points are most spread out — the direction of **maximum variance**.
- **PC2** is the direction of most remaining variance, **constrained to be perpendicular (orthogonal) to PC1**.
- PC3 is perpendicular to both. And so on.

Each PC is a **weighted sum of your 2,934 genes**:

```
PC1 score for a cell = w₁·(gene A) + w₂·(gene B) + w₃·(gene C) + ...
```

Those weights (`w`) are the **loadings**, stored in `adata.varm['PCs']`. Because the weights are fixed for all cells, PCA is a **linear** method — that matters later when we compare it to UMAP.

**The biological reading:** genes don't vary independently. All the T-cell genes go up and down together; all the B-cell genes go up and down together. PCA discovers these **co-varying gene modules** automatically. So a PC is usually not "one gene" — it's "a coordinated programme of hundreds of genes," which in practice often corresponds to a cell type, a cell cycle phase, or a batch.

## What you actually get

```python
sc.tl.pca(adata, use_highly_variable=True, n_comps=100)
print(adata.obsm["X_pca"][:10, :5])
```

`adata.obsm["X_pca"]` is a **5,500 × 100** matrix. Each cell now has 100 coordinates instead of 29,346. The printed numbers (−5.36, 0.188, −2.77…) are just *"where this cell sits along PC1, PC2, PC3…"*. Negative/positive is arbitrary — only the relative positions mean anything.

Note `use_highly_variable=True`: PCA runs on the 2,934 HVGs only, not all genes.

## Choosing how many PCs — the scree plot (cell 20)

```python
sc.pl.pca_variance_ratio(adata, n_pcs=70, log=True)
```

PCs come out ranked: PC1 explains the most variance, PC2 less, PC3 less still. Plot that decay and you get a curve that drops steeply then flattens into a plateau — the **elbow**.

- **Before the elbow:** real, structured biological signal.
- **After the elbow:** each PC explains a tiny sliver, and it's mostly technical noise.

The comment says the plateau is around **50–60 PCs**, which is why cell 24 uses `n_pcs=50`. The notebook computed 100 deliberately to *over-compute* so you can see where the curve flattens — you can't find the elbow if you stop at 50.

The trade-off, stated plainly:
- **Too few PCs** → you merge genuinely distinct cell types, losing rare populations.
- **Too many PCs** → you feed noise into the neighbour graph, blurring structure.
- In practice the result is fairly robust in the 30–50 range. 50 is scanpy's default and a defensible choice.

## Why PCA is essential, not optional

This is the answer to *"why not just run UMAP on the genes directly?"*

1. **Denoising.** Averaging hundreds of co-varying genes into one PC cancels out random per-gene dropout noise. PCA's output is *cleaner* than its input.
2. **Distances become meaningful.** In 29,346 dimensions, almost all pairs of points are roughly equidistant — the "curse of dimensionality." Euclidean distance stops discriminating. In 50 dimensions it works again, and everything downstream (neighbours, clustering, UMAP) is built on distances.
3. **Speed.** 50 columns instead of 2,934.
4. **It's the substrate for everything after.** Clustering, UMAP, neighbour graphs, integration/batch correction (Harmony works *in PC space*) — they all consume `X_pca`, not the gene matrix.

## Reading the PCA plots (cell 21)

```python
sc.pl.pca(adata, color=["SampleGroup", "SampleGroup",
                        "subsets_Mito_percent", "subsets_Mito_percent"],
          dimensions=[(0,1), (2,3), (0,1), (2,3)], ncols=2, size=4)
```

This is a **diagnostic**, not a result. You colour the PC plot by known variables and ask: *what is driving the main axes of variation?*

- Coloured by **`SampleGroup`** → the groups separate on the PCs. Biology *or* batch effect — here the notebook concludes batch effect, which is why module 05 is data integration.
- Coloured by **`subsets_Mito_percent`** → no structure. **This is good news.** If mito% had driven PC1, your top axis of variation would be dying cells, i.e. a QC failure contaminating your biology, and you'd go back and tighten filtering or regress it out.

The general habit: *every time you see unexpected structure, check whether it's technical before calling it biological.*

---

# Step 4 — The neighbour graph (cell 24)

```python
sc.pp.neighbors(adata, n_pcs=50)
```

Easy to skim past, but it's the bridge between PCA and everything visual.

For each cell, find its **k nearest neighbours** (default `n_neighbors=15`) in the 50-dimensional PC space, using Euclidean distance. Draw an edge to each one, weighted by closeness. Do that for all 5,500 cells and you get a **graph**: cells are nodes, "transcriptomic similarity" are edges.

This graph is the shared foundation for:
- **UMAP** (module 03 — this notebook)
- **Leiden/Louvain clustering** (module 04)

Which is why, if you change `n_pcs` or `n_neighbors`, **both your UMAP and your clusters change.** They're not independent results.

---

# Step 5 — UMAP (cells 25–27) ⭐

## What it does

UMAP takes that neighbour graph and **lays it out flat in 2D**, trying to place cells that are connected in the graph close together on the page, and unconnected cells far apart.

The physical intuition: treat edges as **springs** pulling connected cells together, with a repulsive force pushing everything else apart. Let the system settle. Where it lands is your UMAP.

```python
sc.tl.umap(adata)   # writes adata.obsm['X_umap'], shape 5500 × 2
sc.pl.umap(adata, color="SampleName", size=3)
```

## PCA vs UMAP — the comparison to have ready

| | **PCA** | **UMAP** |
|---|---|---|
| Type | **Linear** — each PC is a fixed weighted sum of genes | **Non-linear** — no formula from genes to coordinates |
| Preserves | Global variance structure; distances are real | Local neighbourhoods only |
| Axes mean | "Expression of gene module X" — interpretable | **Nothing.** UMAP1/UMAP2 have no units and no meaning |
| Deterministic | Yes (up to sign) | No — stochastic, changes with random seed |
| Output dims | Choose ~50 | 2 (for the picture) |
| Use for | **Downstream analysis:** clustering, integration, trajectories | **Visualisation only** |
| Reversible | Yes — you can project back to genes | No |

Cell 23 in the notebook states the rule directly: *"The results of a UMAP projection should be used for visualisation only and not for downstream analysis (such as cell clustering)."*

**Why that rule exists:** you cannot cram genuine 50-dimensional structure into 2 dimensions without distortion. Something must be sacrificed. UMAP chooses to preserve *local* structure (who your neighbours are) and sacrifice *global* structure (how far apart the blobs are). So clustering on the 2D UMAP coordinates would cluster on the distortion. Correct practice: **cluster on the neighbour graph / PC space, then colour the UMAP by those cluster labels.** The UMAP is the display, not the analysis.

## What the plots show

- Coloured by **`SampleName`** → each sample forms its own island instead of samples mixing within shared cell-type islands.
- Coloured by **`SampleGroup`** → same pattern at group level.

The conclusion in cell 28: **a large batch effect.** If the biology were driving the picture, you'd see one T-cell blob containing cells from *all* samples. Instead each sample sits apart, so sample identity dominates cell identity. That's the motivation for module 05 (integration/batch correction).

---

# The two exercise questions — these are interview questions

### Q1: `n_neighbors` = 5 vs 50 vs 500

`n_neighbors` sets **how local UMAP's view is** — literally, how many neighbours each cell is connected to in the graph.

- **n_neighbors = 5** → very local view. Many small, fragmented islands. UMAP preserves fine detail but can shatter one real population into pieces and exaggerate tiny differences.
- **n_neighbors = 50** → balanced. Sensible-sized clusters, reasonable global arrangement. (Default is 15.)
- **n_neighbors = 500** → very global view. Everything gets smoothed into a few large continuous blobs; fine structure and rare cell types disappear.

The one-liner: **`n_neighbors` is a local ↔ global dial.** Low = detail, risk of fragmentation. High = broad structure, risk of over-merging.

*(Note: the notebook's cells 33–40 actually vary `n_pcs` (5 vs 50), not `n_neighbors` — that's a separate but related dial. `n_pcs=5` throws away most of the signal, so structure collapses. If you want to actually answer Q1, it's `sc.pp.neighbors(adata, n_neighbors=5, n_pcs=50)` etc.)*

### Q2: if cluster A's centroid is nearer B than C on the UMAP, is A transcriptomically more similar to B than C?

**No.** This is the single most common misreading of a UMAP, and it's a great question to nail.

Reasons:
1. **UMAP does not preserve between-cluster distances.** It optimises local neighbourhoods. The gap between two separated blobs is essentially arbitrary — a by-product of the layout, not a measured quantity.
2. **It's stochastic.** Re-run with a different seed and blobs can land in different relative positions.
3. **Blob size and density are also meaningless** — a big sparse blob isn't "more variable" than a small tight one.

**What you'd do instead:** measure distance in **PC space** — e.g. Euclidean distance between cluster centroids over the 50 PCs, or build a correlation/dendrogram on mean expression profiles (`sc.tl.dendrogram`). That's a real distance in a space that preserves it.

The safe formulation: *"On a UMAP, I trust that points sitting together are genuinely similar. I don't trust the distance between separate clusters, their sizes, or their orientation."*

---

# Things an interviewer may push on

**"Why did you not scale the data before PCA?"**
This notebook goes straight from log-normalised HVGs to `sc.tl.pca`. Scanpy's PCA **centres** the data but does not scale genes to unit variance (Seurat's standard workflow includes `ScaleData` for this). Scaling gives every gene equal weight; not scaling lets highly expressed genes contribute more. Both are defensible — but know which you did. Log-transformation already does much of the variance stabilisation.

**"Why 10% of genes as HVGs?"**
`int(0.1 * adata.n_vars)` = 2,934. It's a convention, typically 2,000–3,000. The result is fairly insensitive in that range. The important thing is that it's a deliberate choice, not a default you can't explain.

**"Is PCA lossy?"**
Yes — keeping 50 of 100 PCs discards variance. But the discarded portion is mostly noise, which is why the output is often *better* for downstream use than the input. Losing information and losing signal are not the same thing.

**"Why does the batch effect show up in PCA too?"**
Because a batch effect *is* real variance — a coordinated shift across many genes. PCA can't tell "coordinated because same cell type" from "coordinated because same 10x run." That's precisely why you need explicit integration methods (Harmony, scVI, MNN) in module 05.

**"What's stored where?"**

| Slot | Contents |
|---|---|
| `adata.X` | log-normalised expression, 5500 × 29346 |
| `adata.var['highly_variable']` | True/False per gene (2,934 True) |
| `adata.obsm['X_pca']` | 5500 × 100 PC coordinates |
| `adata.varm['PCs']` | gene loadings — the weights defining each PC |
| `adata.uns['pca']['variance_ratio']` | the scree plot values |
| `adata.obsp['distances']`, `['connectivities']` | the neighbour graph |
| `adata.obsm['X_umap']` | 5500 × 2 UMAP coordinates |

---

## The 30-second version, if they ask you to summarise

> "Each cell starts as ~29,000 gene measurements, which is too high-dimensional to analyse and mostly noise. I keep the ~2,900 genes that vary more than expected for their expression level, since flat genes can't distinguish cell types. PCA then compresses those into ~50 components — each one a co-varying gene module — chosen from a scree-plot elbow. That's my denoised working representation, and everything downstream uses it. From it I build a k-nearest-neighbour graph in PC space, and UMAP lays that graph out in 2D purely so I can look at it. Colouring the UMAP by sample showed samples separating rather than mixing, which flagged a batch effect and motivated integration as the next step."
