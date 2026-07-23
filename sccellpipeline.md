Scanpy's whole pipeline is really one object being progressively annotated. Understanding **AnnData** first makes everything else click.

## The AnnData mental model

```
adata.X        # matrix: cells × genes (the "current" expression values)
adata.obs      # DataFrame: one row per cell   (QC metrics, clusters, labels)
adata.var      # DataFrame: one row per gene   (mito flag, HVG flag, means)
adata.layers   # alternative matrices, same shape as X ('counts', 'lognorm')
adata.obsm     # per-cell arrays: X_pca, X_umap, X_harmony
adata.obsp     # cell × cell: 'distances', 'connectivities' (the kNN graph)
adata.uns      # unstructured: parameters, DE results, color maps
```

Almost every `sc.*` function mutates `adata` in place and writes its result into one of those slots. The pipeline is a chain where each step reads what the previous one deposited.

Naming convention: `pp` = preprocessing (changes X or adds graph), `tl` = tools (computes and stores results), `pl` = plotting (reads and draws).

---

## 1. Load

```python
adata = sc.read_10x_mtx("filtered_feature_bc_matrix/", var_names="gene_symbols")
adata.var_names_make_unique()
```

X is a sparse matrix of raw UMI counts. Rows = barcodes, not yet real cells — droplets that passed CellRanger's cutoff.

## 2. Quality control

```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=False)
```

This computes per-cell `n_genes_by_counts`, `total_counts`, `pct_counts_mt` into `.obs`. Nothing is filtered yet — you look at the distributions first (violin plots, and a scatter of total_counts vs pct_counts_mt).

The biology behind the three metrics:
- **Low gene count** → empty droplet or dying cell with degraded RNA.
- **Very high counts/genes** → likely a doublet (two cells, one barcode).
- **High mitochondrial %** → the cell membrane ruptured, cytoplasmic mRNA leaked out, mito transcripts stayed trapped in the organelle. So a high ratio means a stressed/lysed cell.

```python
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.pct_counts_mt < 15].copy()
sc.pp.scrublet(adata)          # doublet scores → .obs
```

Thresholds are dataset-dependent, not universal. 5% mito is fine for PBMCs, wrong for cardiomyocytes or hepatocytes, which are genuinely mito-rich. Ideally you set cutoffs per-sample from the distribution rather than hardcoding.

## 3. Normalization

```python
adata.layers["counts"] = adata.X.copy()   # preserve raw counts!
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata                         # snapshot of log-norm, all genes
```

`normalize_total` divides each cell by its total counts and rescales — this removes sequencing-depth differences, since a cell with 20k UMIs isn't twice as transcriptionally active as one with 10k, it was just captured better.

`log1p` does `log(x+1)`. Counts are roughly log-normal with a huge dynamic range; logging makes distances behave, stops a few enormously expressed genes from dominating, and turns fold-changes into additive differences.

The `layers["counts"]` copy matters — some methods (scVI, `seurat_v3` HVG flavor, some DE tests) require raw integer counts and will silently misbehave on transformed data.

## 4. Feature selection

```python
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", batch_key="sample")
```

Most of your ~20-30k genes are either off or uniformly expressed housekeeping genes — they contribute noise, not signal. Scanpy bins genes by mean expression and picks the ones with the highest variance *relative to other genes at similar expression level* (because variance scales with mean; without the binning you'd just select highly-expressed genes). Result is a boolean `.var["highly_variable"]`.

`batch_key` selects genes variable within each batch and takes the consensus, which avoids picking up genes that just track batch.

## 5. Scaling

```python
sc.pp.scale(adata, max_value=10)
```

Z-scores each gene to mean 0, variance 1, so PCA weights a rare transcription factor as heavily as a ribosomal gene. `max_value` clips extreme outliers. Note this makes X dense and produces negative values — which is why you keep `.raw`.

`sc.pp.regress_out(["total_counts", "pct_counts_mt"])` exists but is slow and often does more harm than good; most current workflows skip it.

## 6. PCA

```python
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
```

Writes `.obsm["X_pca"]`. This does two things at once: denoises (each PC is a weighted average over many genes, averaging away dropout noise) and drops dimensionality from 2000 genes to ~30-50 components. Everything downstream operates in PC space, not gene space.

Choose `n_pcs` from the elbow in the variance plot. Too few loses rare populations; too many just adds noise, though the effect is fairly forgiving.

## 7. Integration (only if multiple samples)

```python
sce.pp.harmony_integrate(adata, key="sample")   # → .obsm['X_pca_harmony']
```

Harmony iteratively nudges PC coordinates so that clusters are batch-balanced. Alternatives: BBKNN (corrects the graph instead), scVI (deep generative model, uses raw counts, generally strongest for hard cases), Combat (linear, weakest).

Crucial caveat: correct for *technical* batch, not biology. If your "batches" are treatment vs control, integrating can erase the effect you're studying.

## 8. Neighbor graph — the pivot point

```python
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)   # or use_rep="X_pca_harmony"
```

For each cell, find its k nearest neighbors in PC space, then convert distances into edge weights using UMAP's fuzzy-simplicial-set procedure — each cell's distances get scaled by its own local density, so the graph works in both dense and sparse regions. Stored in `.obsp["connectivities"]`.

This is the most consequential step in the pipeline. **Both UMAP and Leiden read from this same graph.** That's why your clusters look coherent on the UMAP — it's not independent validation, it's two views of one object. If the graph is wrong, both are wrong together.

## 9. Embedding

```python
sc.tl.umap(adata)     # → .obsm['X_umap']
```

UMAP optimizes a 2D layout so that graph neighbors stay close and non-neighbors are pushed apart. It preserves local structure reasonably well. It does **not** preserve global distances — cluster sizes and the gaps between distant clusters are essentially meaningless. Read it as a topology sketch, not a map.

## 10. Clustering

```python
sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2)
```

Leiden partitions the graph to maximize modularity — more within-community edges than expected at random — with a guarantee (that Louvain lacked) that every community is internally connected. Labels land in `.obs["leiden"]`.

`resolution` is the knob: higher → more, smaller clusters. There is no correct value. Standard practice is to sweep 0.2–2.0 and pick the level where clusters correspond to markers you can actually defend. Clustering algorithms will always return clusters, including from noise.

## 11. Marker genes

```python
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)
```

One-vs-rest test per cluster per gene, results into `.uns`. Wilcoxon is the sensible default — t-tests assume normality that scRNA-seq violates, and `logreg` is a different (multivariate) question.

Run this on log-normalized values, not scaled ones. If you scaled `.X`, pass `use_raw=True` (default when `.raw` is set) so fold-changes stay interpretable.

The statistical caveat everyone glosses over: you clustered and tested on the same data, so the p-values are anti-conservative by construction. Treat them as a ranking, not as inference.

## 12. Annotation

```python
adata.obs["cell_type"] = adata.obs["leiden"].map({
    "0": "CD4 T", "1": "CD14 Monocyte", ...
}).astype("category")
```

Match markers against literature or a reference atlas, or automate with CellTypist/scArches. Expect iteration: you'll often find a cluster that's clearly two cell types and go back to raise the resolution, or subcluster that group alone.

---

## Where things actually go wrong

**Losing track of what's in `.X`.** After the pipeline above, `.X` holds scaled z-scores — negative values, HVGs only if you subsetted. Plotting expression from that, or feeding it to a counts-based method, produces confident nonsense. Always know which representation you're on. Many people now keep `.X` as log-normalized permanently and pass scaled data to PCA only.

**Treating the UMAP as the result.** It's a lossy projection of the graph. Two clusters sitting apart may be a rendering artifact; validate with markers and, if it matters, PAGA.

**Filtering before looking.** Applying someone else's thresholds can silently delete a real population — plasma cells legitimately have low gene counts, some cell types legitimately have 30% mito.

Minimal version of the whole thing:

```python
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
sc.pp.filter_cells(adata, min_genes=200); sc.pp.filter_genes(adata, min_cells=3)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata); sc.tl.leiden(adata, resolution=1.0)
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
```

If you tell me what your data looks like — species, tissue, number of samples, 10x or something else — I can point out which of these defaults you'd want to change.
