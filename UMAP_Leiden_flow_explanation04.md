Really good question, and it's the single most common conceptual knot in the whole field. The short answer: **a UMAP shows you that groups exist. It doesn't tell you which group any given cell is in.**

Those are two different deliverables.

---

## They produce completely different things

| | `sc.tl.umap` (notebook 03) | `sc.tl.leiden` (notebook 04) |
|---|---|---|
| **Output** | 2 numbers per cell — *coordinates* | 1 label per cell — *a group ID* |
| **Stored in** | `ad.obsm['X_umap']` (3500 × 2) | `ad.obs['leiden_res0.3']` (3500 × 1) |
| **Type** | continuous floats | categorical: `'0'`, `'1'`, `'2'`… |
| **Answers** | "what does the landscape look like?" | "which cells belong together?" |
| **For** | your **eyes** | your **code** |

That last row is the whole thing. A UMAP is a picture. **You cannot do arithmetic on a picture.**

Look at what notebook 04 does with the Leiden labels — none of it is possible with UMAP coordinates:

```python
ad.obs['leiden_res0.3'].value_counts()                    # how many cells per group?
sc.tl.rank_genes_groups(ad1, groupby="label")              # markers per group
ad1.obs['cellType'] = ad1.obs['label'].map(annotation)     # name each group
pd.crosstab(ad1.obs['SampleName'], ad1.obs['cellType'])    # composition per patient
```

Every one of those needs a **discrete label column**. `groupby=` requires something to group *by*. You cannot write `groupby='X_umap'`.

> **Analogy:** notebook 03 handed you the Marauder's Map. You can see four dense knots of students in the castle and think *"clearly there are houses here."* But no name on that map says **Gryffindor**. Leiden is the Sorting Hat walking through and stamping a house onto every single student. Only then can you ask "how many Gryffindors are there?", "what do Gryffindors have in common?", or "did Gryffindor grow this year?"
>
> **The map shows structure. The Hat assigns membership.** You need both, and they are not the same act.

---

## The part that makes it click: they're siblings from the same parent

Here's what notebook 04 actually does, and why the two stay consistent with each other:

```
                X_corrected  (the batch-corrected 40-dim space)
                             │
                  sc.pp.neighbors(use_rep="X_corrected")     ← cell 14
                             │
                   obsp['connectivities']
                    ⬅ THE KNN GRAPH — the real object of interest
                             │
                    ┌────────┴────────┐
                    ▼                 ▼
             sc.tl.leiden        sc.tl.umap
              cells 16–18           cell 20
                    │                 │
                    ▼                 ▼
          obs['leiden_res*']     obsm['X_umap']
             the LABELS           the PICTURE
```

**Both read the same KNN graph.** Leiden partitions it; UMAP lays it out in 2D. They're two different renderings of one underlying structure — which is exactly why the Leiden clusters land as contiguous patches on the UMAP rather than confetti.

> **Analogy:** one Marauder's Map, two uses. Leiden reads it and writes a list of who's in which gang. UMAP reads it and draws the floor plan. Same source, different products.

### ⚠️ And the sharp corollary — an interview favourite

**Because they share the graph, your clusters looking clean on the UMAP is *not* independent validation of your clustering.** It's circular. You built both from the same connectivities, so of course they agree. Real validation is **marker genes** — cells in cluster 11 all expressing `LYZ`, `CST3`, `TYROBP` is external evidence. Pretty blobs are not.

This trips up a lot of people who present a tidy UMAP as proof the clustering is right.

---

## Why can't you just eyeball the blobs?

Three reasons, all of which bite in practice:

1. **Boundaries are ambiguous.** Where exactly does one blob end? Two people drawing lassos get different answers, and neither is reproducible.
2. **UMAP lies about structure.** It can **split one cell type into two visual islands** purely from initialisation, and it can **smear two distinct types into one blob** because 2D isn't enough room. Leiden works in 40 dimensions, where the separation actually exists.
3. **It doesn't scale.** Hand-drawing 23 clusters (your resolution 1.5 result) across 3,500 cells, then re-doing it at four resolutions? No.

And the deep version: **Leiden sees things UMAP cannot show.** Two clusters may be genuinely distinct along dimension 27 of `X_corrected` and sit right on top of each other in the 2D projection. Leiden splits them correctly; your eye never could.

---

## The bit specific to your course — notebook 03's UMAP was thrown away

This is the detail that probably drove the question, and it's worth seeing clearly.

**Notebook 03's UMAP was a diagnostic, not a product.** Look at what it did:

```python
sc.pp.neighbors(adata, n_pcs=50)   # ← no use_rep → defaults to UNCORRECTED X_pca
sc.tl.umap(adata)
```

Then the very next markdown cell says:

> *"We can see that the different samples have quite a large batch effect"*
> *"So lets proceed to batch correction and data integration"*

**That UMAP existed to reveal a problem.** Its entire purpose was to show you that cells were grouping by *patient* instead of by *cell type*. It was never meant to be clustered on — clustering it would have given you 7 clusters named after 7 patients.

And notebook 04 doesn't reuse it. **Cell 20 computes a brand-new UMAP** from the corrected graph, overwriting `X_umap` entirely (the silent overwrite I flagged earlier).

So the real sequence is:

```
03:  uncorrected PCA → graph → UMAP    → 👁 "there's a batch effect"   → DISCARD
05:  batch correction                   → X_corrected
04:  X_corrected → graph → ├→ LEIDEN   → labels      ⬅ the actual product
                           └→ UMAP     → new picture ⬅ to display the labels on
```

**Notebook 04 doesn't replace notebook 03's UMAP with Leiden. It builds both again, on data that's finally trustworthy.**

> **Analogy:** notebook 03's map was drawn while everyone was still under Polyjuice. Useful — it's how you *discovered* everyone was Polyjuiced. But you don't sort students from that map. You wait for the potion to wear off (integration), redraw the map, *and then* bring in the Hat.

---

## The one-sentence answer

**UMAP is a diagnostic and a display surface; Leiden is the measurement.** Notebook 03 used UMAP to find the batch effect. Notebook 04 clusters the corrected space to produce the labels that every downstream analysis — markers, annotation, DE in 06, abundance in 07 — literally cannot run without.

If you deleted `obsm['X_umap']` from notebook 04, every result would be identical and you'd just be working blind. Delete `obs['leiden_res0.3']` and notebooks 04, 06 and 07 all stop dead.
