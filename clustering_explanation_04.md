No errors in your run, and the outputs are all there — so I can walk you through it with the real numbers. Here's the whole notebook, step by step.

---

# The story of this notebook in one breath

You are given a pile of **3,500 cells** from children with leukaemia plus healthy controls. Someone has already cleaned them, normalised them, reduced them to a handful of dimensions, and removed batch effects. Your job in this notebook is two things:

1. **Sorting** — split those 3,500 cells into groups of "cells that look alike" (clustering).
2. **Naming** — work out *what each group actually is* (annotation).

That's it. This is **the Sorting Hat chapter** of scRNA-seq. Every cell walks up, the Hat looks at its transcriptome, and shouts a house. And just like at Hogwarts, the hard part isn't the sorting — it's that nobody tells you what "Hufflepuff" means until you go look at what those students have in common.

The dataset is the **Caron et al. childhood acute lymphoblastic leukaemia (B-ALL)** data. Hold that thought — it explains a genuinely confusing result at the end.

---

# PART 0 — Setup (cells 1–6)

## Cell 1 — the imports

```python
import scanpy as sc
import pandas as pd
import os
```

`scanpy` is your wand. `pandas` is your ledger for tables. `os` is for touching the filesystem.

> **Analogy:** `import scanpy as sc` is Ollivander handing you the wand. From here on, `sc.` is you saying *"wand, do the thing."* The `sc.pp.`, `sc.tl.`, `sc.pl.` prefixes are three different kinds of spell, and knowing which is which is half of learning scanpy:
>
> | Prefix | Means | Does what | Analogy |
> |---|---|---|---|
> | `sc.pp.` | **p**re**p**rocessing | *Changes your data* | Potion-making — the cauldron contents are permanently different |
> | `sc.tl.` | **t**oo**l**s | *Computes and stores a result* | Casting a spell that adds a new object to the room |
> | `sc.pl.` | **pl**otting | *Draws a picture, changes nothing* | Looking in the Mirror of Erised. Pretty. Harmless. |

**Tip:** when you're lost about whether a line is going to mutate your data, look at the middle letters. `pp` = danger, `pl` = safe.

## Cell 2 — the output folder

```python
prefix_output = "/mnt/e/RNA-seq-courses/.../Data/results/04"
#os.mkdirs(prefix_output, exist_ok=True)
```

This just stores a path in a variable. The `mkdirs` line is commented out — **and it's misspelled anyway**. The real function is `os.makedirs`, not `os.mkdirs`. If you ever uncomment it, it will crash.

**Two things an experienced person spots instantly here:**
- `prefix_output` is defined and then **never used anywhere in the notebook.** Nothing is saved. You run 20 minutes of clustering and it all dies when you close the kernel. This is the single most common junior mistake in scRNA-seq — see the tips section for the fix.
- The path is `/mnt/e/...`, which is a **WSL/Linux** path for your Windows `E:\` drive. So this notebook only runs inside your WSL or your `.sif` container, not in native Windows Python. Worth knowing before you spend an hour debugging "file not found."

## Cells 3–4 — cosmetics

```python
pd.set_option('display.max_rows', 6)
pd.set_option('display.max_columns', 10)
plt.rcParams['figure.dpi'] = 80
```

"When you print a table, show me 6 rows and 10 columns, not 3,500." And "make plots this crisp." Zero effect on the analysis. Purely so your notebook doesn't turn into a wall of numbers.

**Tip:** bump `figure.dpi` to 120–150 when you're actually reading UMAPs. At 80, small clusters genuinely disappear and you will mis-annotate. Keep it low while exploring (fast), raise it for the plots you'll stare at.

## Cell 6 — loading the data

```python
ad = sc.read_h5ad(".../Caron_batch_corrected.500.h5ad")
```

`.h5ad` is the standard single-cell file on disk. `ad` is now your **AnnData** object — the one thing you must truly understand in this whole course.

The `.500` in the filename means **500 cells were randomly sampled per sample** to keep the demo fast. 3,500 cells total ÷ 500 = seven samples' worth.

**Interviewers love this:** *"This is a downsampled dataset. What does downsampling to 500 cells per sample cost you?"* Answer: you lose **rare cell types**. A population at 0.5% frequency is ~2–3 cells in 500 — below any clustering algorithm's ability to form a community. Downsampling is fine for pipeline demos and catastrophic for rare-cell discovery.

---

# PART 1 — Meet the AnnData object (cells 7–10)

## Cell 8 — `ad`

```
AnnData object with n_obs × n_vars = 3500 × 33102
    obs:  'Sample', 'Barcode', 'SampleName', 'SampleGroup', 'sum', 'detected',
          'subsets_Mito_sum', 'subsets_Mito_detected', 'subsets_Mito_percent',
          'total', 'sizeFactor'
    var:  'ID', 'Symbol', 'Chromosome'
    obsm: 'X_corrected', 'X_pca', 'X_tsne', 'X_tsne_corrected', 'X_umap', 'X_umap_corrected'
```

**This printout is the most important thing on the screen.** Read it every single time you load or modify data. Let me decode it.

> **Analogy — AnnData is the Room of Requirement.** One room, many compartments, and it grows new compartments as you need them:
>
> | Slot | What lives there | Shape | Room of Requirement version |
> |---|---|---|---|
> | `.X` | The expression matrix — the actual numbers | cells × genes | The room itself |
> | `.obs` | One row per **cell**: metadata | 3500 rows | The register of every student who walked in |
> | `.var` | One row per **gene**: metadata | 33102 rows | The catalogue of every object on the shelves |
> | `.obsm` | Per-cell **multi-dimensional** things: PCA, UMAP | 3500 × k | The Pensieve — compressed views of the same memories |
> | `.obsp` | Cell × cell **relationships** | 3500 × 3500 | The Marauder's Map — who is standing next to whom |
> | `.uns` | Unstructured leftovers: settings, colours, results | anything | The box of junk in the corner that turns out to matter |
> | `.layers` | Alternative versions of `.X` (counts, lognorm, scaled) | cells × genes | Time-Turner copies of the room at different moments |
> | `.raw` | A frozen snapshot of `.X` from earlier | cells × genes | A Horcrux — a preserved earlier self |

**Memorise the mnemonic:** `obs` = **obs**ervations = **cells** (rows). `var` = **var**iables = **genes** (columns). Every single beginner mixes these up at least once. The `m` in `obsm` means "matrix", the `p` in `obsp` means "pairwise".

Notice what is **missing** from that printout: no `layers`, no `raw`, no `uns`. Which tells you `.X` is the only copy of the expression data you have. Remember this — it bites in cell 27.

Notice what **is** there: `X_pca`, `X_tsne`, `X_umap` *and* `X_corrected`, `X_tsne_corrected`, `X_umap_corrected`. Two parallel sets — before and after batch correction. The `_corrected` versions came from an R workflow (fastMNN, most likely) and were converted to `.h5ad`.

## Cell 9 — `ad.obs.columns`

Just listing the per-cell metadata. Quick translation of what each means:

| Column | Meaning |
|---|---|
| `Sample`, `SampleName` | Which patient/library this cell came from |
| `SampleGroup` | The biological condition (ALL subtype vs healthy control) |
| `sum` / `total` | Total UMI counts in this cell — its "library size" |
| `detected` | How many distinct genes were detected in this cell |
| `subsets_Mito_percent` | % of reads from mitochondrial genes — the classic dying-cell flag |
| `sizeFactor` | The normalisation scaling factor computed upstream in R (scran's pooled deconvolution) |

`sum`, `detected`, `sizeFactor` are R/Bioconductor naming (`scater`/`scran`). Native scanpy would call these `total_counts`, `n_genes_by_counts`. **The mixed naming is your tell that this object was born in R.** Interviewers notice when you notice.

## Cell 10 — `ad.var.head()`

```
                              ID           Symbol Chromosome
MIR1302-2HG      ENSG00000243485      MIR1302-2HG       chr1
ENSG00000238009  ENSG00000238009  ENSG00000238009       chr1
ENSG00000239945  ENSG00000239945  ENSG00000239945       chr1
```

**⚠️ This is one of the most important quiet details in the notebook and it's easy to skim past.**

Row 1's index is `MIR1302-2HG` — a human-readable gene symbol. Rows 2–4's index is `ENSG00000238009` — a raw Ensembl ID. Why the inconsistency? Because **those genes have no assigned symbol**, so the pipeline fell back to the Ensembl ID.

Why you must care: later, every marker gene lookup is done **by index name** — `color='CD79A'`, `'LYZ'`, `'CST3'`. That works because those genes happen to have symbols. But:

- If a gene is only stored as `ENSG...`, `color='YOURGENE'` → **KeyError**, and you conclude the gene isn't expressed when actually it's just named differently.
- Gene symbols are **not unique**. Several Ensembl IDs can map to the same symbol. Ensembl IDs are unique; symbols are not. If duplicates existed here, scanpy would have appended suffixes or thrown a warning.

> **Analogy:** Ensembl IDs are the gene's *legal name on the Ministry of Magic registry*. Symbols are its *nickname*. "Padfoot" is friendlier than "Sirius Orion Black", but two different people can be called Padfoot and some people have no nickname at all. Index on the legal name, display the nickname.

**Tip — the defensive one-liner, use it before every marker plot:**
```python
genes = ["CD79A", "LYZ", "CST3", "CD3D"]
missing = [g for g in genes if g not in ad.var_names]
print("MISSING:", missing)
```
Cell 49 at the very end of this notebook does exactly this. It should have been at the top.

---

# PART 2 — Look before you leap (cells 11–13)

```python
sc.pl.embedding(ad, basis='X_pca',  color='SampleGroup', title='PCA')
sc.pl.embedding(ad, basis='X_corrected', color='SampleGroup', title='PCA')
# ... same pattern for tsne, then umap
```

Six plots, all the same question asked three ways: **"has the batch effect been dealt with?"**

Each pair is *uncorrected on the left, corrected on the right*, coloured by `SampleGroup`. What you want to see:

- **Left (uncorrected):** colours sitting in separate blobs. Each sample forming its own island. That's a batch effect — cells are grouping by *where they were processed*, not by *what they are*.
- **Right (corrected):** colours intermixed within each cluster. Same cell type from different patients now sitting together.

> **Analogy — Polyjuice Potion.** The batch effect is everyone in the room having drunk Polyjuice. Two T cells from two different patients look like totally different creatures, purely because of a technical artifact. Batch correction is the potion wearing off. You need to *watch it wear off* before you trust anything downstream.
>
> **ACOTAR version:** uncorrected data is the seven Courts, sealed behind walls, each speaking its own dialect. Batch correction is tearing the walls down so you can finally see that the Winter Court healer and the Summer Court healer are the same *kind* of person.

## The three embeddings, and why all three exist

| Method | What it is | Trustworthy for | Analogy |
|---|---|---|---|
| **PCA** | Linear. Finds the directions of greatest variance. | **Distances — real and meaningful.** The workhorse. | An honest but boring map. Scale is correct everywhere. |
| **t-SNE** | Non-linear. Preserves *local* neighbourhoods. | Seeing that groups exist | A map that magnifies the village you're in and shrinks everything else |
| **UMAP** | Non-linear. Local + a bit of global structure. Fast. | Same, plus rough between-cluster layout | The Marauder's Map — brilliant for "who's near whom," lying about "how far apart" |

**🔴 The single biggest conceptual trap in the entire field, and interviewers ask it constantly:**

> *"In this notebook, clustering uses `X_corrected` and the plots use `X_umap_corrected`. Why not just cluster on the UMAP coordinates?"*

**Answer:** UMAP is a **2D cartoon for human eyes only**. It compresses 40+ dimensions into 2, and in doing so it distorts distance non-linearly, invents empty space between clusters, and splits continuous gradients into fake islands. Clustering, distance metrics, and neighbour graphs must all run in the **high-dimensional corrected space** (`X_corrected`). UMAP is the *display layer*, never the *compute layer*.

> **Analogy:** the Marauder's Map shows Snape three corridors from Harry. It does not tell you the true walking distance — there are moving staircases, secret passages, and the geometry is a lie. You use the map to *see* that Snape is nearby. You'd never use it to *calculate* how long the walk takes.

Corollary you should internalise: **two clusters sitting far apart on a UMAP are not necessarily more different than two sitting adjacent.** And cluster *sizes* on a UMAP mean essentially nothing.

---

# PART 3 — Building the neighbour graph (cells 14–15)

```python
sc.pp.neighbors(ad, use_rep="X_corrected", n_neighbors=10, n_pcs=40)
```

**This one line is the foundation the rest of the notebook stands on.** If you only deeply understand one cell, make it this one.

What it does: for every one of the 3,500 cells, find its 10 nearest neighbours in the corrected 40-dimensional space, and record who-is-next-to-whom. The result is a **k-nearest-neighbour (KNN) graph** — cells are dots, "you're one of my 10 closest" is a line.

> **Analogy — the Marauder's Map, built from scratch.** You're not drawing a picture yet. You're writing down, for every person in the castle, the ten people standing closest to them. Once you have that list for everybody, the social groups — Gryffindor Quidditch team, the Slug Club, Dumbledore's Army — emerge from the connection pattern alone. Nobody had to tell you the groups existed. You inferred them from who hangs out with whom.

Breaking down the arguments:

- **`use_rep="X_corrected"`** — "measure distances in the *batch-corrected* space." ⭐ **This is the most important argument in the cell.** Leave it out and scanpy defaults to `X_pca` (uncorrected), and you'd cluster patients instead of cell types. Every single downstream result would be wrong, and the plots would look perfectly plausible.
- **`n_neighbors=10`** — how many neighbours. This is a **resolution knob in disguise**:
  - Low (5–10) → tight, local, sensitive graph → more, smaller clusters, more noise-driven splits
  - High (30–50) → smoother, globally-connected graph → fewer, broader clusters
  - Rough rule: scale with cell count. 10 is fine for 3,500 cells; for 100,000 cells you'd want 15–30.
- **`n_pcs=40`** — use only the first 40 dimensions of `X_corrected`. Later dimensions mostly encode noise.

### ⚠️ Three things experienced people check on this line

1. **`n_pcs` with `use_rep` is a known confusion point.** When you pass `use_rep`, scanpy slices that representation to its first `n_pcs` columns. So `n_pcs=40` here means "the first 40 columns of `X_corrected`" — *not* "40 principal components." If `X_corrected` had only 30 columns, `n_pcs=40` would be meaningless or error. Always check `ad.obsm['X_corrected'].shape` before setting this. Nobody in the notebook does, which is exactly why you should.
2. **`n_pcs=40` is a number someone chose, not a number someone computed.** The honest way is to look at an elbow plot (`sc.pl.pca_variance_ratio`) upstream and pick where the variance flattens. 40 is a very common convention-by-inertia. If asked "why 40?", the correct answer is "it was inherited from the upstream notebook, and I'd verify it against the variance ratio."
3. **The TensorFlow message in the output is noise.** `oneDNN custom operations are on...` — that's some library in the environment initialising TensorFlow. Nothing to do with your clustering. Learning to recognise irrelevant log spam is an underrated skill; beginners lose hours debugging messages that were never errors.

## Cell 15 — proof it worked

```
    uns:  'SampleGroup_colors', 'neighbors'
    obsp: 'distances', 'connectivities'
```

Two new compartments appeared in the Room of Requirement:

- **`obsp['distances']`** — the raw KNN distances. Sparse: only ~10 entries per row out of 3,500.
- **`obsp['connectivities']`** — the **weighted, fuzzified** graph. This is what Leiden actually consumes. Rather than a hard "neighbour: yes/no," each edge carries a strength between 0 and 1.
- **`uns['neighbors']`** — the receipt: which parameters were used.

> **Analogy:** `distances` is "Ron is 4 metres from Harry." `connectivities` is "Ron and Harry are *close*, weight 0.9" — a judgement, not a measurement. Leiden cares about relationships, not metres.

**Tip:** `uns['neighbors']` is a paper trail. Six months later, when you can't remember whether you clustered on corrected or uncorrected space, `ad.uns['neighbors']['params']` tells you. This is why you save the object rather than the plots.

---

# PART 4 — The Sorting Hat: Leiden clustering (cells 16–19)

## Cell 16 — one clustering

```python
sc.tl.leiden(ad, resolution=1, key_added="leiden_res1")
```

**Leiden** is a community-detection algorithm. Give it the graph, and it partitions the cells into groups that are densely connected internally and sparsely connected to each other.

> **Analogy — the Sorting Hat, but it never read *Hogwarts: A History*.** The Hat doesn't know the houses exist. It has no list of "brave," "clever," "loyal." All it has is the Marauder's Map of who-stands-near-whom. It looks for knots of people who cluster tightly together, and declares each knot a house. It then hands you **numbered** houses: 0, 1, 2, 3... It has absolutely no idea that group 3 is "Gryffindor." **Naming them is your job, and that's Part 8 of this notebook.**

### `resolution` — the one dial that changes everything

| Resolution | Behaviour | Sorting Hat analogy |
|---|---|---|
| 0.1–0.3 | Few, broad clusters | "Four houses. Done." |
| 0.5–1.0 | Moderate | "Four houses, split by year group." |
| 1.5–3.0 | Many, fine clusters | "Every dormitory is its own house. Also Neville is his own house." |

Higher resolution → more clusters. Always. It's the algorithm's willingness to split.

### The FutureWarning in the output — read it, it's telling you something real

```
FutureWarning: In the future, the default backend for leiden will be igraph instead of leidenalg.
To achieve the future defaults please pass: flavor="igraph" and n_iterations=2.
directed must also be False to work with igraph's implementation.
```

**Don't ignore this.** It means your scanpy version's default backend is changing. If you re-run this notebook in a year on a newer scanpy, **you will get different cluster numbers** — and if a figure in your thesis says "cluster 7 is the NK T cells," your figure quietly becomes wrong.

**The fix — use this form from now on:**
```python
sc.tl.leiden(
    ad,
    resolution=1.0,
    key_added="leiden_res1.0",
    flavor="igraph",      # the future default — pin it explicitly
    n_iterations=2,       # required companion to flavor="igraph"
    directed=False,       # required for the igraph implementation
    random_state=0,       # reproducibility
)
```

### 🎯 Interview favourite: Leiden vs Louvain

*"Why does everyone use Leiden now instead of Louvain?"*

**Answer:** Louvain has a documented failure mode — it can produce communities that are **internally disconnected**. You get a "cluster" whose members aren't actually all linked to each other, which is biologically nonsensical. Leiden (Traag, Waltman & van Eck, 2019) adds a refinement step that **guarantees well-connected communities**, converges to a better partition, and runs faster. Louvain is legacy; Leiden is the default in both Scanpy and modern Seurat.

## Cell 18 — clustering at four resolutions

```python
resolutions = [0.3, 0.5, 1.0, 1.5]
for res in resolutions:
    key = f"leiden_res{res}"
    sc.tl.leiden(ad, resolution=res, key_added=key)
```

**This is the right instinct, and it's the most valuable habit in this whole notebook.** There is no single correct resolution. So you compute several, then compare and choose with your eyes and your biology.

The `f"leiden_res{res}"` is an **f-string** — Python builds the name on the fly. `res=0.3` → `"leiden_res0.3"`. Clean, and it means the resolution is baked into the column name so you never lose track.

> **Analogy:** instead of asking the Sorting Hat once and living with it, you ask it four times at four sensitivity settings, write down all four answers, and *then* decide which sorting actually carves reality at the joints.

### ⚠️ The duplicate-column trap, right here in the notebook

Cell 16 made `leiden_res1` (resolution=1). Cell 18 makes `leiden_res1.0` (resolution=1.0). **These are the same clustering, computed twice, stored under two names.** Look at cell 19's output — both columns are sitting there.

```
obs: ..., 'leiden_res1', 'leiden_res0.3', 'leiden_res0.5', 'leiden_res1.0', 'leiden_res1.5'
```

Harmless here — but this is exactly how people end up with 40 near-identical columns, plotting the wrong one and not noticing.

> **Analogy:** Fred and George. Functionally identical, differently labelled, and you will absolutely mix them up at the worst possible moment.

**Tip:** always build keys from a single format string (`f"leiden_res{res:.1f}"` → always `1.0`, never `1`), and never hand-type a cluster key twice.

---

# PART 5 — Redrawing the map (cells 20–22)

## Cell 20 — `sc.tl.umap(ad)`

One short line. **And it silently does something that catches out even experienced users.**

```python
sc.tl.umap(ad)
```

This computes a **fresh UMAP** from the neighbour graph you built in cell 14 (the corrected one). Good — that's what you want, a UMAP consistent with the graph the clusters came from.

### 🔴 But: it writes to `ad.obsm['X_umap']` — overwriting the one that came from the file.

Look at cell 8's output again. `X_umap` was already there, imported from R. After cell 20, **it's gone**, replaced by a scanpy UMAP built on corrected data.

Why this matters concretely: **cell 24** then plots `X_umap` beside `X_umap_corrected` as though it were an uncorrected-vs-corrected comparison. It isn't anymore. Both are now built on corrected data — one by scanpy, one by R. The comparison in cell 24 is no longer the comparison it appears to be.

> **Analogy — the Vanishing Cabinet.** You put the original map in the cabinet, cast a spell, and a *different* map comes back out with the same label on it. Nothing errors. Nothing warns. The label still says `X_umap`. You will trust it.

**This is the class of bug that never announces itself.** Silent in-place overwrites of `obsm`, `obs`, and `.X` are responsible for a large share of irreproducible single-cell figures.

**The fix — in recent scanpy, name your output:**
```python
sc.tl.umap(ad, key_added="X_umap_scanpy")   # leaves X_umap untouched
```
Or, if your version doesn't support `key_added`, defend it manually before you cast:
```python
ad.obsm["X_umap_original"] = ad.obsm["X_umap"].copy()   # keep the imported one
sc.tl.umap(ad)
```

**Tip:** make a habit of printing `ad` before *and* after any `sc.tl.` / `sc.pp.` call while you're learning. Diff the two printouts in your head. You'll catch overwrites and typos within seconds instead of two notebooks later.

## Cell 21 — the resolution comparison plot

```python
sc.pl.umap(ad, color=['leiden_res0.3','leiden_res0.5','leiden_res1.0','leiden_res1.5'],
           wspace=0.2, frameon=False)
```

Four panels, **the same UMAP coordinates every time**, coloured four different ways. The dots never move; only the colouring changes. That's the whole point — you're seeing four different opinions about the same landscape.

What to look for as resolution climbs: do new clusters carve out **genuinely separate territory on the map**, or do they just slice an existing blob into arbitrary pie wedges? Territory = probably real. Pie wedges = over-clustering.

- `wspace=0.2` — horizontal gap between panels
- `frameon=False` — drop the axis box. **This is more than cosmetic:** UMAP axis values are meaningless, so displaying them invites people to read distances off them. Hiding the frame is honest.

## Cell 22 — the sanity-check plot (do not skip this one)

```python
sc.pl.umap(ad, color=['leiden_res0.3','SampleGroup','SampleName'], wspace=0.2, frameon=False)
```

Three panels: clusters, then condition, then individual sample. You are asking: **"are my clusters biology, or are they just patients?"**

- Cluster that is **one colour** in the `SampleName` panel → that cluster is **one patient**. Suspicious.
- Cluster that is a **mix** of `SampleName` colours → cell type shared across individuals. Reassuring.

> **Analogy:** you've sorted the students and are now checking whether "Gryffindor" secretly just means "everyone who took the 9 a.m. train." If your houses map perfectly onto train times, you haven't discovered houses — you've rediscovered the train timetable.

### 🔴 The plot twist specific to this dataset

In this leukaemia data, **some clusters genuinely are single-patient — and that's real biology, not a batch effect.** Each child's leukaemic B-cell blasts are a *clonal* population carrying that child's own mutations. They legitimately don't resemble any other patient's blasts. Meanwhile the healthy T cells, NK cells and monocytes cluster beautifully across everyone.

This is why the final annotation in cell 50 has **six separate B-cell clusters** (c1, c2, c3, c4, c5, c12) but only **one** T cluster and **one** monocyte cluster.

> **Analogy:** the healthy immune cells are the Hogwarts houses — the same four exist at every school, and you can compare across schools. The leukaemic blasts are **Horcruxes**: each one is unique to the individual who made it, and expecting them to cluster together across patients is a category error.

**🎯 Interview gold.** *"You see a cluster containing cells from only one patient. Batch effect or biology?"* The answer that gets you hired: **"It depends on the biology, and here's how I'd distinguish them."** Batch effect suspects — the cluster is driven by technical covariates (library size, mito %, ribosomal content), it persists after correction across *all* cell types not just one, and its markers are housekeeping genes. Real biology suspects — it has a coherent, patient-specific marker programme; the patient has a known clonal disease; it survives correction while other populations mix fine. In cancer, tumour cells clustering per-patient is the **expected** result. Over-correcting them away destroys your signal.

---

# PART 6 — Choosing a resolution (cells 23–25)

## Cell 23 — silhouette score

```python
from sklearn.metrics import silhouette_score
for res in [0.3, 0.5, 1.0, 1.5]:
    key = f"leiden_res{res}"
    score = silhouette_score(ad.obsm['X_corrected'], ad.obs[key])
    print(f"Silhouette Score for resolution {res}: {score:.4f}")
```

Your actual output:

```
Silhouette Score for resolution 0.3: 0.2940
Silhouette Score for resolution 0.5: 0.2647
Silhouette Score for resolution 1.0: 0.1826
Silhouette Score for resolution 1.5: 0.1715
```

**What silhouette measures, per cell:** *"Am I closer to my own clustermates than to the nearest other cluster?"* Formally, for each cell: `a` = mean distance to its own cluster, `b` = mean distance to the nearest *other* cluster, and `s = (b − a) / max(a, b)`. Average over all cells.

- **+1** — perfectly tight, perfectly separated
- **0** — sitting right on a boundary
- **−1** — in the wrong cluster entirely

> **Analogy:** at the Yule Ball, everyone checks whether the people at their own table are more their kind than the people at the next table over. Score near 1 = the tables are real friendship groups. Score near 0 = the seating was arbitrary and everyone's mingling anyway.

**Note this well:** it's computed on `ad.obsm['X_corrected']`, the **high-dimensional** space — *not* the UMAP. That's correct and deliberate. Computing silhouette on UMAP coordinates gives you beautiful, meaningless numbers, because UMAP was *designed* to make blobs look separate. It's a common and serious error.

### 🔴 The big caveat nobody mentions when they show you this metric

Look at the pattern: the score falls monotonically as resolution rises. **That is not evidence that 0.3 is biologically correct.** Silhouette is *structurally biased* toward fewer clusters. Two facts make it nearly useless as a standalone resolution-picker on scRNA-seq:

1. **More clusters → clusters are closer together → `b` shrinks → score drops.** Almost mechanically. If you tested resolution 0.05 you'd get two clusters and an even better score. By silhouette's logic, the ideal analysis has two clusters — which is obviously wrong.
2. **Silhouette assumes round, well-separated, Euclidean blobs.** Biology gives you continuous differentiation gradients, branching trajectories, and elongated manifolds. A perfectly correct clustering of a developmental continuum scores *badly* on silhouette.

Also worth knowing: it's O(n²) in memory and time. Fine on 3,500 cells. On 200,000 cells it will kill your kernel. Use `sample_size=5000` on large data.

**🎯 Interview trap, and it's a good one:** *"Silhouette says resolution 0.3 is best. Do you go with 11 clusters?"*

The answer that separates candidates: **"No — I'd use silhouette as one weak input among several, and let biology decide."** Then list your real toolkit:

| Approach | What it actually tells you |
|---|---|
| **Marker interpretability** ⭐ | Can I name every cluster with a known, coherent marker programme? If two clusters share identical markers, merge them. **This is the real criterion.** |
| **Cluster stability** | Bootstrap-resample the cells, re-cluster, measure agreement. Clusters that survive resampling are real; ones that dissolve are noise. |
| **Clustree / Sankey across resolutions** | Watch clusters split as resolution climbs. Clean bifurcations = real substructure. Cells sloshing between clusters = over-clustering. |
| **ARI / AMI between resolutions** | Quantifies how much the partition actually changed. |
| **Minimum viable size** | A "cluster" of 8 cells out of 3,500 is usually a doublet nest or an ambient-RNA artifact. |
| **Domain knowledge** | If you're studying B-cell subsets, you *need* resolution high enough to split them, silhouette be damned. |

And the pragmatic honest answer: **deliberately over-cluster slightly, then merge by marker identity.** Splitting a cluster you failed to split is impossible after the fact; merging two you over-split is trivial. Asymmetric risk → bias toward splitting.

## Cell 24 — comparing embeddings again

```python
sc.pl.embedding(ad, basis='X_umap', color=['leiden_res0.3','SampleName'], title='UMAP')
sc.pl.embedding(ad, basis='X_umap_corrected', color=['leiden_res0.3','SampleName'], title='UMAP_corrected')
```

Note the warning in the output:

```
WARNING: The title list is shorter than the number of panels. Using 'color' value instead for some plots.
```

Two colours were requested but only one `title` given, so scanpy titled the second panel itself. Harmless — but it's a good illustration of scanpy's style: it warns rather than crashing, and **those warnings are worth reading.** Silent recoveries are how wrong figures get made.

And remember from cell 20: `X_umap` is no longer the uncorrected embedding. This comparison isn't what it looks like.

## Cell 25 — how many clusters did we get

```python
for res in [0.3, 0.5, 1.0, 1.5]:
    num_clusters = ad.obs[f"leiden_res{res}"].nunique()
    print(f"Number of clusters for resolution {res}: {num_clusters}")
```

```
resolution 0.3 → 11 clusters
resolution 0.5 → 13 clusters
resolution 1.0 → 18 clusters
resolution 1.5 → 23 clusters
```

The dial in action: **resolution 0.3 → 11 clusters, resolution 1.5 → 23.** Same cells, same graph, same algorithm. Just a different willingness to split. Put these numbers next to the silhouette scores and the trade-off is right there in front of you.

**Tip — a more informative version, because cluster *count* hides the thing that matters:**
```python
for res in [0.3, 0.5, 1.0, 1.5]:
    vc = ad.obs[f"leiden_res{res}"].value_counts()
    print(f"res {res}: {len(vc)} clusters | smallest={vc.min()} | largest={vc.max()}")
```
A jump from 18 to 23 clusters where all five new clusters have <15 cells is a very different story from five new clusters of 150 cells each. **Always look at the size distribution, not just the count.**

---

# PART 7 — First look at markers, and a normalisation landmine (cells 26–30)

## Cell 26 — plotting two marker genes

```python
sc.pl.umap(ad, color=['leiden_res0.3','CD79A',"LYZ","SampleName"], wspace=0.2, frameon=False)
```

The first attempt at *naming* the clusters. Two classic markers:

- **`CD79A`** — part of the B-cell receptor complex. **B cells.**
- **`LYZ`** — lysozyme, an antibacterial enzyme. **Monocytes / myeloid cells.**

You're looking for: which numbered cluster lights up for which gene. That's how a number becomes a name.

> **Analogy:** the Sorting Hat handed you houses 0 through 10 with no names. Now you walk into house 7 and everyone's wearing a red-and-gold scarf. Aha — *that's* Gryffindor. `CD79A` is the scarf. **Marker genes are the house colours of cell biology.**

### ⚠️ But this plot is showing you raw counts

Look at the very next cell:

```python
sc.pp.normalize_total(ad, target_sum=None)
sc.pp.log1p(ad)
```

**Normalisation happens *after* cell 26.** So the `CD79A` and `LYZ` colouring in cell 26 is raw UMI counts, unnormalised, un-logged. Which means:

- A cell with a big library (20,000 UMIs) looks "high" for every gene, just because it was sequenced more deeply.
- Without log transformation, a handful of extreme cells dominate the colour scale and everything else looks uniformly dark. You lose all the mid-range signal.

Cell 28 re-plots the same genes after normalisation. **Compare cells 26 and 28 side by side — that visual difference is the entire argument for normalisation, for free, in a demo you already ran.** Genuinely one of the most instructive accidents in this notebook.

## Cell 27 — normalisation and log transform

```python
sc.pp.normalize_total(ad, target_sum=None)
sc.pp.log1p(ad)
```

**`normalize_total`** — divide each cell's counts by that cell's total, then rescale. This removes **library-size / sequencing-depth** differences, so you're comparing *proportions* of a cell's transcriptome rather than absolute counts.

> **Analogy:** two students write essays, one 500 words and one 5,000. You can't compare raw counts of the word "however." You compare *rate per thousand words*. Normalisation converts counts to rates.

**`target_sum=None`** is the interesting choice. Your options:

| Value | Every cell rescaled to | Note |
|---|---|---|
| `1e4` | 10,000 counts | "CP10K". The scanpy tutorial default, most common in the literature |
| `1e6` | 1,000,000 | "CPM". Familiar from bulk RNA-seq |
| `None` | **the median total count across cells** | Data-driven; avoids inflating a shallow dataset onto an artificially large scale |

`None` is a defensible, arguably better choice — it keeps numbers in the range your data actually occupies. Just know that **your values are no longer comparable to someone else's CP10K-normalised values.**

**`log1p`** — compute `log(x + 1)` for every value.

Why the `+1`? Because `log(0)` is negative infinity, and your matrix is mostly zeros. The `+1` maps 0 → 0 cleanly. Hence "log1p."

Why log at all? Three real reasons:
1. **Variance stabilisation.** In count data, highly expressed genes have much larger variance. Log shrinks that so a few loud genes don't dominate every distance calculation.
2. **Makes fold-changes additive.** A difference of 1 on a log₂ scale = a 2-fold change, everywhere on the scale. Biology thinks multiplicatively; log lets your linear methods keep up.
3. **Kills the long tail.** Without it, one cell with 8,000 copies of a mitochondrial transcript bends your entire analysis around itself.

> **Analogy:** the Richter scale. Without log, a magnitude-9 earthquake is 100,000,000× a magnitude-1, and you can't plot them on the same axis. With log, they're 9 and 1. Same information, a scale a human can reason about.

### 🔴 Three landmines in these two innocent lines

**1. Normalising twice is a silent disaster.** If `.X` had already been normalised upstream, running `normalize_total` again produces garbage that looks completely plausible. `log1p` at least tries to protect you — it stamps `ad.uns['log1p']` and newer scanpy warns if you try again — but `normalize_total` has **no such guard**. Nothing stops you.

**How to check before you cast:**
```python
import numpy as np
x = ad.X[:20].toarray() if hasattr(ad.X, "toarray") else ad.X[:20]
print("max value:", x.max())
print("all integers?", np.allclose(x, np.round(x)))
# integers + max in the hundreds/thousands → raw counts, safe to normalise
# floats + max around 5-10             → already lognormalised, STOP
```

**2. The raw counts are now destroyed.** Remember cell 8: no `layers`, no `raw`. `.X` was the only copy of the counts, and these two lines overwrote it in place. Counts are gone unless you re-read the file.

This matters because several methods **require raw counts** and will give wrong answers on lognormalised input — notably many differential-expression models (negative-binomial based: DESeq2, edgeR, glmGamPoi), `scVI`/`totalVI`, and `sc.pp.highly_variable_genes(flavor="seurat_v3")`.

**The habit that will save you — do this immediately after loading, every time, forever:**
```python
ad.layers["counts"] = ad.X.copy()   # ⭐ before ANY sc.pp. call
sc.pp.normalize_total(ad, target_sum=None)
sc.pp.log1p(ad)
ad.raw = ad                          # freeze the lognormalised full gene set too
```
Costs you one line and a little memory. Saves you re-running the whole pipeline.

> **Analogy:** a Horcrux, for once used responsibly. Before you go into the fight that will change you, you split off a piece of yourself and hide it somewhere safe. `layers["counts"]` is your Horcrux.
>
> **ACOTAR version:** Feyre Under the Mountain doesn't get to go back to who she was before. Your `.X` doesn't either. Copy it while you still can.

**3. Order of operations.** Normalise, *then* log. Always. `log1p` then `normalize_total` is meaningless — you'd be dividing logged values by a sum of logged values, which corresponds to no statistical model at all.

## Cells 29–30 — violin plots

```python
sc.pl.violin(ad, 'LYZ',   groupby='leiden_res0.3', color='leiden_res0.3', use_raw=False)
sc.pl.violin(ad, 'CD79A', groupby='leiden_res0.3', color='leiden_res0.3', use_raw=False)
```

One violin per cluster, showing the **full distribution** of that gene's expression across the cells in that cluster. Width at any height = how many cells sit at that expression level.

> **Analogy:** the UMAP colour plot tells you "house 7 looks reddish on average." The violin tells you **whether all 200 students wear a red scarf, or whether 15 wear very red scarves and 185 wear none.** That distinction decides whether the cluster is truly a B-cell cluster or a mixed cluster with some B cells hiding in it.

**Why violins matter more than UMAP colouring for annotation:** a mean can be created by a minority. The shape can't lie to you the same way.

- **Bimodal violin** (blob at zero, blob up high) → your cluster contains two populations. Consider raising resolution.
- **Unimodal and high** → coherent, confidently expressing population. Good marker, good cluster.
- **Everything at zero with a few dots up top** → not a marker for this cluster. Ambient RNA or dropout noise.

### `use_raw=False` — small argument, important meaning

Many scanpy plotting and testing functions **default to using `.raw` if it exists**. `use_raw=False` says "no, use `.X`."

Here `.raw` doesn't exist (cell 8 confirmed it), so the argument is redundant — *but writing it explicitly is good practice*, because on an object where `.raw` *does* exist, forgetting it means you silently plot a different version of your data than you think.

**🎯 Interview question:** *"What is `.raw` and when does it burn people?"* Answer: `.raw` is a frozen snapshot of the expression matrix, conventionally stored **after** normalisation and log-transform but **before** subsetting to highly-variable genes and before scaling. It exists so you can still plot and test genes that were filtered out of `.X`. It burns people because **the defaults are inconsistent between functions** — `sc.pl.violin`, `sc.pl.dotplot` and `sc.tl.rank_genes_groups` all consult `.raw` by default, so if you stored *raw counts* in `.raw` instead of lognormalised values, you will silently run statistical tests on unnormalised data and never see a warning.

> **Analogy:** the Pensieve. `.raw` holds a memory of your data as it was at a specific moment. Extremely useful — and dangerous if you can't remember *which* moment you bottled.

A note on `color='leiden_res0.3'` in these two calls: `sc.pl.violin` has no `color` parameter, so this falls through to the underlying seaborn call. It ran without complaint, but it isn't doing what the name suggests. **The real argument for colouring violins is `palette=`.** Minor, but worth knowing so you don't copy the pattern and then wonder why it does nothing.

---

# PART 8 — The annotation half begins (cells 31–38)

## Cell 32 — load the pre-clustered object

```python
ad1 = sc.read_h5ad(".../Caron_clustered.500.h5ad")
```

A **different object**, `ad1`. Same cells, but clustering was already done upstream (in R), and the answer lives in `ad1.obs['label']`.

Why switch? Partly teaching convenience — everyone in the room now has identical cluster numbers, so the course's annotations match. **But be aware:** `ad1`'s clusters are *not* the ones you computed in cells 16–18. Different tool, different parameters, different numbering. Don't carry insights about "cluster 6" between `ad` and `ad1`.

## Cell 33 — `ad1.X.todense()`

```
matrix([[0., 0., 0., ..., 0., 0., 0.],
        [0., 0., 0., ..., 0., 0., 0.],
        ...
```

Two lessons in one line of output.

**Lesson 1 — scRNA-seq data is mostly zeros.** Typically 90–95% of entries. Partly biology (most genes are off in most cells), partly technical (**dropout** — the transcript was there, but at ~20% capture efficiency the sequencer never saw it). This is why `.X` is stored as a **sparse matrix**: it only records the non-zero values and their coordinates.

> **Analogy:** Hermione's beaded bag. From the outside it's tiny. It holds a tent, a library, and Grimmauld Place's worth of belongings — because it only stores what's actually in it, not the empty space. `.todense()` is dumping the entire bag onto the floor of the Great Hall.

**Lesson 2 — ⚠️ never do this on real data.** Let's do the arithmetic: 3,500 × 33,102 = **115,857,000 values × 8 bytes = ~930 MB**. Nearly a gigabyte, from one keystroke, to look at a corner of a matrix that's almost entirely zeros. Scale to a realistic 100,000-cell dataset and you're asking for **26 GB**. Your kernel dies instantly with no useful error message.

**Tip — inspect sparse matrices safely:**
```python
print(ad1.X.shape, ad1.X.dtype, type(ad1.X))
print("sparsity:", 1 - ad1.X.nnz / (ad1.X.shape[0] * ad1.X.shape[1]))
print(ad1.X[:5, :5].toarray())     # densify a 5×5 corner, not the universe
```
Also: `.todense()` returns a `np.matrix` (a deprecated, quirky subclass). Prefer **`.toarray()`**, which gives you a proper `ndarray`.

## Cell 34 — `ad1.raw`

Output: **nothing.** `.raw` is `None`.

A deliberate check, and exactly the right one to make. It's *why* `use_raw=False` appears throughout, and why `rank_genes_groups` in cell 39 will operate on `.X`. **Check `.raw` on every object you didn't create yourself.**

## Cell 35 — `ad1.obs["label"]`

```
1_CGACTTCGTCCAGTTA-1      1
11_ACGGAGAGTTAAGATG-1    11
Name: label, Length: 3500, dtype: category
Categories (12, object): ['1', '2', '3', ..., '9', '10', '11', '12']
```

Three details here, all of which matter:

1. **12 clusters, labelled `'1'` through `'12'`.** Note they start at **1** — the R convention. Scanpy's own `leiden` starts at **0**. Mixing R-derived and scanpy-derived labels in one project is a reliable off-by-one bug factory.
2. **They are strings, not integers.** `'1'`, not `1`. So `ad1.obs['label'] == 1` silently matches nothing, while `== '1'` works. This is why cell 50's annotation dictionary uses string keys.
3. **`dtype: category`** — pandas `Categorical`, not plain strings. Memory-efficient, and it lets you fix the plotting order. But categoricals have sharp edges: assigning a value that isn't in `.categories` fails, and `.map()` with an incomplete dictionary produces `NaN` rather than an error. Cell 50 is standing right next to that hole.
4. **The barcode index is informative.** `1_CGACTTCGTCCAGTTA-1` — the prefix before the underscore is the sample, the middle is the cell barcode, the `-1` is the 10x GEM well suffix. The prefix exists because barcodes are **not unique across samples** — the same 16-mer appears in different libraries as different cells. Merging samples without prefixing barcodes silently collapses distinct cells together.

## Cell 36 — normalise `ad1` too

```python
sc.pp.normalize_total(ad1, target_sum=None)
sc.pp.log1p(ad1)
```

Same as cell 27, applied to the new object. **Necessary** — everything that follows (marker plots, Wilcoxon tests, violins, dotplots) assumes lognormalised data. Skip this and your fold-changes and p-values are all computed on raw counts.

And the same landmines apply: no counts layer saved, no `.raw` set. This object's raw counts are also now gone.

## Cells 37–38 — first look at the labelled map

```python
sc.pl.embedding(ad1, basis='X_umap_corrected', color=['label','CD79A',"LYZ"], legend_loc="on data")
sc.pl.embedding(ad1, basis='X_umap_corrected', color=['CST3'])
```

`legend_loc="on data"` prints cluster numbers **directly on the clusters** instead of in a side legend. **Use this constantly during annotation.** With 12 clusters and a colour legend off to the side, you spend your whole life eye-matching pastel shades. With numbers on the blobs, you read the answer straight off the plot.

`CST3` (cystatin C) is another **monocyte/dendritic marker**, checked as a second opinion alongside `LYZ`.

> **Analogy:** you don't identify someone from one detail. Scarf colour *and* the badge *and* which common room they walk toward. **One marker is a hypothesis. Three concordant markers is an identification.** Interviewers care enormously about this habit — "I'd confirm with a panel, not a single gene" is the answer of someone who has actually done this work.

---

# PART 9 — Finding markers computationally (cells 39–43)

## Cell 39 — the statistical test

```python
sc.tl.rank_genes_groups(ad1, groupby="label", method="wilcoxon",
                        key_added="rank_genes_groups_wilcoxon")
```

So far you've been *guessing* markers from prior knowledge. Now you flip it: **let the data tell you which genes define each cluster.**

For each of the 12 clusters, this tests every one of the 33,102 genes: *"is this gene higher in this cluster than in all the other cells combined?"* Then it ranks them. That's **one-vs-rest** differential expression, 12 times over.

> **Analogy:** instead of walking into house 7 looking for red scarves, you line up all 12 houses and ask *"what does house 7 own that the other eleven don't?"* The answer comes back as a ranked list, and **you never had to know what to look for in advance.** This is how you find cell types you weren't expecting — and how genuinely novel populations get discovered.

### The arguments

- **`groupby="label"`** — the cluster column to compare across
- **`method="wilcoxon"`** — the **Wilcoxon rank-sum test** (= Mann–Whitney U)
- **`key_added="..."`** — store results under a named key in `.uns`

### Why Wilcoxon, and not t-test?

Wilcoxon is **rank-based**: it throws away the actual values and keeps only the ordering. That makes it robust to exactly what scRNA-seq data is — wildly non-normal, zero-inflated, heavy-tailed. A t-test assumes approximate normality, which single-cell expression violates enthusiastically.

> **Analogy:** judging a Quidditch season. The t-test cares about *how many points* each team scored, so one freak 400-point match distorts the whole table. Wilcoxon only cares *who finished above whom*. The outlier match can't dominate.

Your method options:

| Method | Character | When |
|---|---|---|
| `'wilcoxon'` | Rank-based, robust | **Default recommendation.** Scanpy's own docs prefer it |
| `'t-test'` | Fast, assumes normality | Quick look only |
| `'t-test_overestim_var'` | t-test with a conservative variance | Slightly safer t-test |
| `'logreg'` | Multivariate logistic regression | Interesting **ranking** — considers genes jointly, not one at a time. No p-values. |

**Tip:** with sparse data there are enormous numbers of tied values (all those zeros). Passing `tie_correct=True` gives a more accurate Wilcoxon statistic, at some speed cost. Worth it for results going into a paper.

### 🔴 The deepest conceptual issue in the whole notebook — "double dipping"

This is the question that separates people who ran a tutorial from people who understand what they ran.

**🎯 Interview question:** *"Your top marker gene has an adjusted p-value of 1e-121. How much do you believe it?"*

**Answer: as a p-value, not at all. As a ranking, quite a lot.**

Here's the circularity. You used the expression data to *define* the clusters. Now you're using **the same data** to test whether those clusters differ in expression. But the clusters were *constructed* to differ in expression — that's what clustering means. So the null hypothesis ("no difference between these groups") was already violated by how the groups were built. The p-values are **anti-conservative by orders of magnitude** and are not valid inference. This is known as **double dipping**, **data snooping**, or the **circular analysis problem**.

> **Analogy:** you sort students into houses *specifically by* how brave they are. Then you run a statistical test showing that Gryffindors are significantly braver than Hufflepuffs, p < 10⁻¹²¹, and publish it as a discovery. You didn't discover anything. You measured your own sorting criterion and called it a finding.

**What to do about it:**
- Use the p-values as a **ranking device** to surface candidate genes. That's legitimate and useful.
- Never report them as evidence that the clusters are real.
- For genuine inference on the *conditions* (ALL vs healthy), do **pseudobulk** DE: sum counts per sample per cell type, then run DESeq2/edgeR on those with the **sample** as the unit of replication, not the cell. This is what notebook 06 is for, and it's the single biggest step up in statistical rigour in the whole course.
- If you truly need valid per-cluster p-values, the advanced answer is **count-splitting** or data-thinning (split each cell's counts into independent train/test halves — cluster on one, test on the other), or methods built for the problem such as ClusterDE or the TN test. Mentioning any of these by name in an interview is a strong signal.

Also note: **cells are not independent replicates.** 3,500 cells from 7 patients is n=7, not n=3,500. Treating cells as replicates inflates significance enormously — the "pseudoreplication" problem, and a well-documented source of false positives in the single-cell literature.

## Cell 40 — the dotplot

```python
sc.pl.rank_genes_groups_dotplot(ad1, groupby="label", standard_scale="var",
                                n_genes=5, key="rank_genes_groups_wilcoxon")
```

**The single most information-dense plot in single-cell biology.** Learn to read it and you can annotate a dataset in minutes.

Each dot encodes **two** numbers at once:

| Visual channel | Encodes |
|---|---|
| **Colour** (intensity) | *Mean expression* of that gene, among the cells that express it |
| **Size** (diameter) | *Fraction of cells in the cluster* that express it at all |

That second channel is what makes it powerful. Consider:

- **Big and dark** → nearly every cell in the cluster expresses it, strongly. **This is a real marker.** ⭐
- **Small and dark** → a handful of cells express it very highly. Could be a rare subpopulation, could be doublets, could be ambient contamination. **Not a cluster marker.**
- **Big and pale** → everybody has a little. Housekeeping-ish. Useless for discrimination.
- **Absent/tiny** → off in this cluster.

> **Analogy:** you're assessing which house owns "having a pet owl." Colour = *how impressive the owls are*. Size = *how many students actually have one*. One student with a magnificent eagle owl (small, dark) tells you nothing about the house. Eighty percent of students each with a decent owl (big, dark) tells you everything.

**`standard_scale="var"`** — scale each gene's expression to the 0–1 range **across clusters**. Without it, a single enormously-expressed gene (a ribosomal protein, say) saturates the colour scale and every other gene looks blank. With it, you see *relative* pattern per gene, which is what you want for annotation. The cost: you lose absolute magnitude, so you can no longer tell a highly-expressed gene from a barely-expressed one.

**`n_genes=5`** — top 5 genes per cluster. 12 clusters × 5 = 60 columns. Push it to 10 and the plot becomes unreadable. 3–5 is the sweet spot.

Note the warning in your output:
```
WARNING: dendrogram data not found (using key=dendrogram_label).
Running `sc.tl.dendrogram` with default parameters.
```
Scanpy auto-clustered the *clusters themselves* into a tree (so similar clusters sit next to each other on the plot) and is telling you it used defaults. **That dendrogram is genuinely useful** — clusters that branch together are transcriptionally similar and are your **merge candidates** if you suspect over-clustering. The warning suggests running `sc.tl.dendrogram(ad1, groupby='label')` yourself so you control the parameters. Good advice; take it.

## Cell 41 — cluster 11's top genes

```python
sc.get.rank_genes_groups_df(ad1, group="11", key="rank_genes_groups_wilcoxon").head(5)
```

```
    names     scores  logfoldchanges          pvals      pvals_adj
0  LGALS1  23.814713        7.195841  2.351423e-125  7.783679e-121
1     LYZ  23.706751        9.338234  3.071662e-124  5.083908e-120
2  TYROBP  23.641224        7.701289  1.453029e-123  1.603272e-119
3    CST3  23.495899        5.844524  4.492093e-122  3.717432e-118
4  S100A6  23.450804        5.011444  1.297160e-121  8.587719e-118
```

`sc.get.*` pulls results out of `.uns` as a tidy **pandas DataFrame** — filterable, sortable, exportable. This is how you work with markers programmatically instead of squinting at plots.

**Reading the columns:**

| Column | Meaning |
|---|---|
| `names` | Gene symbol |
| `scores` | The Wilcoxon z-score. **This is what the ranking is by.** Higher = more cluster-specific |
| `logfoldchanges` | **log₂** fold change, cluster vs rest. `9.34` means ~2⁹·³ ≈ **650-fold** higher |
| `pvals` | Raw p-value |
| `pvals_adj` | **Benjamini–Hochberg FDR-adjusted** p-value. Always read this one, never `pvals` |

**And now the biology — this is a textbook-perfect result:**

`LYZ` (lysozyme), `TYROBP` (a myeloid signalling adaptor), `CST3` (cystatin C), `LGALS1` (galectin-1), `S100A6` — every one of these is a **monocyte/myeloid** gene. Four independent markers pointing the same direction, with log₂FC of 5–9 (30× to 650× enrichment). That's not a marginal call. **Cluster 11 is monocytes**, and cell 50 labels it exactly that.

> **Analogy:** this is a cluster wearing the scarf, the badge, the tie, *and* carrying the house-specific textbook. No ambiguity. When annotation looks like this, trust it and move on.

**⚠️ One caution on `logfoldchanges`:** scanpy computes it as log₂ of the ratio of mean expressions, with a tiny pseudocount (1e-9). When a gene is essentially **zero in the reference group**, that ratio explodes and you get spectacular but unstable fold-changes. A log₂FC of 15 usually means "off in the rest" rather than "1000× more interesting than log₂FC of 5." Always sanity-check enormous fold-changes against `pct_nz_group`/`pct_nz_reference` (fraction of cells expressing) — scanpy can add these with `pts=True`:

```python
sc.tl.rank_genes_groups(ad1, groupby="label", method="wilcoxon", pts=True,
                        key_added="rank_genes_groups_wilcoxon")
```
**Turn `pts=True` on by default.** It gives you the dotplot's size channel as a number you can filter on, and it's how you separate a true marker from a handful of outlier cells.

## Cells 42–43 — confirming visually

```python
sc.pl.embedding(ad1, basis='X_umap_corrected', color=['label',"LYZ"], legend_loc="on data")
sc.pl.violin(ad1, 'LYZ', groupby='label', color='label', use_raw=False)
```

**The discipline worth copying:** the statistics said `LYZ` marks cluster 11, and rather than believing the table, you go look. The UMAP shows *where*; the violin shows *the distribution* — that `LYZ` is high in essentially all of cluster 11 and near-zero everywhere else, not just high on average.

> **Analogy:** Moody's advice. **Constant vigilance.** A p-value of 1e-124 is not a substitute for looking at the data. Plenty of published, beautifully-significant markers turn out to be driven by nine cells.

---

# PART 10 — The exercise, and the most instructive result in the notebook (cells 44–48)

## Cell 44 — the setup

```python
# CD3D suggests cluster 6 and 7 are T cells
sc.pl.embedding(ad1, basis='X_umap_corrected', color=['label',"CD3D"], legend_loc="on data")
sc.pl.violin(ad1, 'CD3D', groupby='label', color='label', use_raw=False)
```

`CD3D` is part of the T-cell receptor complex — **the** canonical T-cell marker. Clusters 6 and 7 both express it. The exercise: *why are they two clusters and not one?*

## Cell 45 — cluster 6's markers, and a genuine red flag

```python
sc.get.rank_genes_groups_df(ad1, group="6", key="rank_genes_groups_wilcoxon").head(5)
```

```
   names     scores  logfoldchanges          pvals      pvals_adj
0  RPS27  33.152874        1.672657  5.148473e-241  1.704247e-236
1  RPL21  32.862972        1.679237  7.434287e-237  1.230449e-232
2  RPS29  32.591293        1.629365  5.448890e-233  6.012305e-229
3  RPL34  31.237585        1.743909  3.291744e-214  2.724083e-210
4  RPL32  30.895584        1.514293  1.369251e-209  9.064992e-206
5  ...
```

**Stop and look at this properly. It is the single most educational output in the notebook.**

Every top-5 gene is `RPS*` or `RPL*` — **ribosomal protein** genes. And notice the numbers:

- p-values are **astronomically small** (1e-241!) — smaller than the monocyte cluster's
- but log₂FC is only **~1.5–1.7** — a mere 3-fold change

**That combination — insanely significant, barely different — is the classic signature of a technical artifact, not a cell type.** With 3,500 cells, a tiny but perfectly consistent shift becomes hyper-significant. Significance measures your confidence that a difference exists; **effect size** measures whether it matters. Here, confidence is enormous and the effect is small.

Compare cluster 11: log₂FC of 5–9. *That's* a cell type.

> **Analogy:** "Gryffindors are statistically significantly taller than the school average, p < 10⁻²⁴¹." Also, the difference is 4 millimetres. You have measured something absolutely real and completely uninformative. **Gryffindor is not the house of tall people.**

**Why ribosomal genes show up as "markers" so often — the standard suspects:**
1. **Cell size / total RNA content.** Bigger, more metabolically active cells have more ribosomes. Ribosomal fraction correlates with library size, which correlates with cell volume.
2. **Normalisation residue.** Ribosomal transcripts are a huge share of every cell's reads (often 20–40%), so anything that shifts library composition shifts them.
3. **Cell state, not cell type.** Proliferating and activated cells genuinely upregulate translation machinery. That's real biology — but "actively translating" is a *state* shared by many cell types, not an identity.
4. **Dissociation and ambient RNA** effects also push these around.

**The standard mitigation — filter housekeeping families out before ranking:**
```python
import numpy as np
mask = ~ad1.var_names.str.match(r'^(RP[SL]|MT-|MRP[SL]|HB[AB])')
sc.tl.rank_genes_groups(ad1[:, mask], groupby="label", method="wilcoxon",
                        key_added="rgg_filtered")
```
(`RPS`/`RPL` = ribosomal proteins, `MT-` = mitochondrial, `MRPS`/`MRPL` = mitochondrial ribosomal, `HBA`/`HBB` = haemoglobin. That last one you'd **keep** in this dataset — see cell 49.)

### 🌟 And now the beautiful twist, specific to this dataset

The Caron et al. paper that produced this data is, notably, **about ribosomal protein expression** — they reported it as a genuine source of intra-individual heterogeneity among leukaemic cells. So in *this* dataset, ribosomal variation is not purely noise.

**Which is exactly the lesson.** "Ribosomal genes = artifact" is a good heuristic, not a law. The right move is never to delete them reflexively — it's to **recognise them, ask whether they could be biology in your system, and check whether they correlate with technical covariates**:
```python
# does the "marker programme" just track library size?
print(ad1.obs.groupby('label')[['sum','detected','subsets_Mito_percent']].median())
```
If cluster 6 also has conspicuously higher `sum` and `detected`, you have your answer.

> **ACOTAR analogy:** a gene that looks like a technical artifact in every other dataset turning out to be the paper's central biological finding is very much a Tamlin-vs-Rhysand situation. The obvious reading is not always the correct one, and the only way to know is to look at the full context rather than the first impression.

**🎯 Interview question, very commonly asked:** *"Your top differentially-expressed genes are all ribosomal. What do you do?"* A strong answer covers: check **effect size, not just p-value**; check whether the pattern tracks technical covariates (total counts, genes detected, mito %); consider whether it reflects a proliferative *state* rather than an identity; re-rank with housekeeping families excluded to see what's underneath; and — the part most candidates miss — **don't assume it's an artifact if the biology of your system makes translation rate meaningful** (leukaemia, stem cells, plasma cells, anything highly proliferative).

## Cell 46 — cluster 7's markers

```
       names     scores  logfoldchanges         pvals     pvals_adj
0        B2M  17.629669        1.611062  1.458177e-69  4.826859e-65
1       SRGN  15.214087        3.370955  2.851717e-52  4.719877e-48
2  LINC-PINT  14.845518        3.029438  7.437826e-50  8.206897e-46
3     IFITM1  14.800387        3.123500  1.456405e-49  1.158922e-45
4       IL32  14.788010        4.372953  1.750532e-49  1.158922e-45
```

**A much more interpretable result, and notice the fold-changes are bigger (3–4.4) than cluster 6's.**

- **`B2M`** — beta-2-microglobulin, part of MHC class I. High in cytotoxic lymphocytes.
- **`SRGN`** — serglycin, the proteoglycan that packages **cytotoxic granules**. ⭐ This is the tell.
- **`IFITM1`** — interferon-induced. Activated lymphocyte.
- **`IL32`** — pro-inflammatory cytokine, characteristic of activated T/NK cells.
- **`LINC-PINT`** — a long non-coding RNA. Less informative; illustrates that marker lists contain genes nobody can interpret, and that's normal.

**Cytotoxic granules + interferon response + inflammatory cytokine = a cytotoxic lymphocyte.** Cluster 7 is labelled **NK T** in cell 50, and this supports it.

> **So the exercise's answer:** clusters 6 and 7 are both `CD3D`⁺ T-lineage cells, but 6 is the conventional T-cell pool (driven mostly by translational/metabolic differences) while **7 carries a genuine cytotoxic effector programme.** The Hat split them for a real reason.
>
> **Analogy:** both are Gryffindors. But cluster 7 is Dumbledore's Army — same house, and they're carrying wands drawn.

**Tip — confirm with a proper NK/cytotoxic panel:** `NKG7`, `GNLY`, `PRF1` (perforin), `GZMB`/`GZMA` (granzymes), `KLRD1`, `FCGR3A` (CD16). And to distinguish true NK from NK-T: real NK cells are `CD3D`-**negative**. Since these are `CD3D`⁺, "NK T" (or a cytotoxic/effector T subset) is the more careful call. Cell 49 checks `NKG7`, which is the right instinct.

## Cells 47–48 — plot the top genes automatically

```python
dc_cluster_genes = sc.get.rank_genes_groups_df(ad1, group="6",
                       key="rank_genes_groups_wilcoxon").head(10)["names"]
sc.pl.embedding(ad1, basis='X_umap_corrected',
    color=[*dc_cluster_genes, "label"],
    legend_loc="on data", frameon=False, ncols=3)
```

**A genuinely useful pattern — learn this idiom.** It goes from "the statistics say these genes" to "here's where they actually are on the map" with no hand-typing.

The mechanics worth understanding:
- `.head(10)["names"]` → a pandas Series of 10 gene names
- **`[*dc_cluster_genes, "label"]`** — the `*` is **unpacking**. It spreads the Series into individual elements inside a new list, then appends `"label"`. Result: `['RPS27', 'RPL21', ..., 'label']`. Without the `*` you'd get a nested list and scanpy would reject it.
- `ncols=3` → grid layout, 3 panels per row

**The always-include-`label` trick:** by putting `"label"` last in the colour list, you get the cluster map **in the same figure**, at the same coordinates, as the gene panels. No flipping between figures trying to remember where cluster 6 was. This is a small habit that makes annotation dramatically faster.

**Tip — wrap it in a function and stop retyping:**
```python
def peek_markers(adata, group, n=9, key="rank_genes_groups_wilcoxon"):
    genes = sc.get.rank_genes_groups_df(adata, group=group, key=key).head(n)["names"]
    sc.pl.embedding(adata, basis="X_umap_corrected",
                    color=[*genes, "label"], legend_loc="on data",
                    frameon=False, ncols=3)

peek_markers(ad1, "6")
```
Now annotating 12 clusters is 12 one-liners instead of 12 copy-paste-edit cycles where you'll eventually forget to change the group number and annotate cluster 6 twice. (People do this. Constantly.)

---

# PART 11 — Naming the houses (cells 49–51)

## Cell 49 — the known-marker panel

```python
known_genes = [
    "HBA1",    # erythrocytes
    "CST3",    # monocytes
    "CD3E",    # T cells
    "NKG7",    # NK T cells
    "CD79A",   # B cells
    "MS4A1"    # CD20 B cells
]
known_genes = [gene for gene in known_genes if gene in ad1.var_names]
sc.pl.violin(ad1, known_genes, groupby='label')
```

**This is the payoff cell — the moment 12 numbers become 12 cell types.** One curated marker per expected lineage, all plotted against all clusters at once.

The panel, decoded:

| Gene | Cell type | What it does |
|---|---|---|
| `HBA1` | **Erythrocytes** | Haemoglobin alpha. Red blood cells / erythroid precursors |
| `CST3` | **Monocytes** | Cystatin C. Myeloid |
| `CD3E` | **T cells** | TCR complex. Pan-T |
| `NKG7` | **NK / NK T** | Natural killer granule protein. Cytotoxic |
| `CD79A` | **B cells** | BCR complex. Pan-B |
| `MS4A1` | **CD20⁺ B cells** | This *is* CD20 — the rituximab target. Mature B |

Note the deliberate pairing: `CD79A` marks **all** B-lineage cells including immature blasts, while `MS4A1`/CD20 marks **mature** B cells. Two markers for one lineage, chosen to resolve a **maturation stage**. That's how you go from "B cells" to "which kind of B cells" — which, in a B-ALL dataset, is the entire clinical question.

> **Analogy:** `CD79A` tells you someone is a Gryffindor. `MS4A1` tells you they're a seventh-year Gryffindor rather than a first-year. In leukaemia, that distinction *is* the diagnosis — B-ALL is defined by B-lineage cells arrested at an immature stage.

### The defensive line — copy this habit

```python
known_genes = [gene for gene in known_genes if gene in ad1.var_names]
```

A **list comprehension** filtering out genes that aren't in the dataset. Without it, one missing gene → `KeyError` → the whole plot fails.

**⚠️ But note the trade-off: it fails *silently*.** If `HBA1` weren't present, it vanishes from the panel and you'd never know a lineage went unchecked. Better:

```python
missing = [g for g in known_genes if g not in ad1.var_names]
if missing:
    print(f"⚠️ NOT IN DATASET: {missing}")
known_genes = [g for g in known_genes if g in ad1.var_names]
```

**Why genes go missing, for real:** filtered out during QC for low expression; stored as an Ensembl ID because no symbol was assigned; or a **naming variant** — `HBA1` vs `HBA2`, `MS4A1` vs `CD20`, `FCGR3A` vs `CD16`, `PTPRC` vs `CD45`. That last category is the sneaky one: flow-cytometry names (CD-numbers) and gene symbols are frequently different words for the same protein, and you will look for `CD20` and conclude B cells are absent.

**Tip — search rather than guess:**
```python
print([g for g in ad1.var_names if g.startswith("HBA")])   # HBA1, HBA2
print([g for g in ad1.var_names if "MS4A" in g])
```

## Cell 50 — the annotation dictionary

```python
cell_annotation = {
    "1": "B (c1)",       "2": "B (c2)",
    "3": "B (c3)",       "4": "B (c4)",
    "5": "CD20+ B (c5)", "6": "T (c6)",
    "7": "NK T (c7)",    "8": "Erythrocytes (c8)",
    "9": "Erythrocytes (c9)",  "10": "Erythrocytes (c10)",
    "11": "Monocytes (c11)",   "12": "B (c12)"
}
ad1.obs['cellType'] = ad1.obs['label'].map(cell_annotation)
```

**The Sorting Hat finally learns the house names.** A plain Python dictionary from cluster number → biological identity, applied with `.map()`.

### Three things this dictionary teaches you about real annotation

**1. Multiple clusters, one cell type — and that's correct.** Six B clusters (1, 2, 3, 4, 5, 12). Three erythrocyte clusters (8, 9, 10). Why? Because in a **B-ALL** dataset, each patient's leukaemic B blasts are a clonal population with patient-specific expression. They *should* form separate clusters. Erythrocytes split by maturation stage and by a handful of genes so dominant (haemoglobins) that small shifts separate them.

> This closes the loop with cell 22: those single-patient clusters weren't a failure of batch correction. They were the disease.

**2. The naming convention is excellent — steal it.** `"B (c1)"` keeps **both** the biology *and* the original cluster ID. You can always trace a label back to the cluster that produced it, and when a reviewer asks "which cluster was that?", you already know. Never throw away the cluster number.

**3. Notice what's absent: no "unknown."** Every cluster got a confident name. In real work you will have clusters you cannot identify, and labelling one `"Unknown (c13)"` or `"Low-quality (c14)"` is **not failure — it's honesty.** Forcing a name onto an ambiguous cluster is how doublet clusters get published as novel intermediate cell types. (Which has happened. More than once.)

### ⚠️ The `.map()` trap, sitting right here

`ad1.obs['label']` is a **Categorical**. When you `.map()` a Categorical with a dictionary:

- Any category **missing from your dictionary** silently becomes `NaN`. No error, no warning.
- With 12 clusters and 12 keys this works. Add a 13th cluster at higher resolution, forget to add the key, and you get a `NaN` cell type that quietly propagates into every downstream plot and count.
- **And remember the string/int trap** from cell 35: keys must be `"1"`, not `1`. Integer keys would map *everything* to `NaN` — and you'd get a plot, just an all-grey one.

**The safe version:**
```python
ad1.obs['cellType'] = ad1.obs['label'].map(cell_annotation)
n_bad = ad1.obs['cellType'].isna().sum()
assert n_bad == 0, f"{n_bad} cells unmapped! Missing keys: " \
    f"{set(ad1.obs['label'].cat.categories) - set(cell_annotation)}"
ad1.obs['cellType'] = ad1.obs['cellType'].astype('category')
```

That `assert` is three seconds of typing and will one day save you a day of debugging. The final `astype('category')` matters too — `.map()` can hand back `object` dtype, and some scanpy plotting functions behave differently (continuous vs categorical colouring) depending on dtype.

**Tip — control the plotting order** so related types sit together in legends instead of appearing alphabetically:
```python
order = ["B (c1)","B (c2)","B (c3)","B (c4)","B (c12)","CD20+ B (c5)",
         "T (c6)","NK T (c7)","Monocytes (c11)",
         "Erythrocytes (c8)","Erythrocytes (c9)","Erythrocytes (c10)"]
ad1.obs['cellType'] = ad1.obs['cellType'].cat.reorder_categories(order)
```

## Cell 51 — the final figure

```python
sc.pl.embedding(ad1, basis='X_umap_corrected', color=['cellType',"label","SampleName"],
                legend_loc="on data")
```

Three panels: **named cell types**, the original cluster numbers, and the sample of origin. The complete picture, with the audit trail attached.

The `SampleName` panel is doing real work even now. It's your final check: the B clusters should be patient-dominated (leukaemic clones), and the T/NK/monocyte clusters should be shared across patients (normal immune cells). If the *monocytes* split by patient, you have an uncorrected batch effect and need to go back. **Keep a sample-of-origin panel in your final figure — always.** It's the difference between a result and a claim.

> **The chapter ends with the Hat's work explained.** Twelve numbers, seven patients, 33,102 genes, and a map where every region has a name and a reason.

---

---

# 🎯 Interview question bank

Grouped by how likely you are to be asked. I've marked the ones I'd bet money on.

## Clustering fundamentals

**⭐ Q: Walk me through clustering a single-cell dataset.**
Normalise → log-transform → select highly variable genes → scale → PCA → (batch-correct) → **KNN graph on the corrected reduced space** → **Leiden** at several resolutions → evaluate → annotate by markers. The two details that show you've done it: clustering happens on the *corrected, reduced* representation, and resolution is chosen by marker interpretability, not by a metric.

**⭐ Q: Why cluster on PCA/corrected space instead of the full gene matrix?**
Curse of dimensionality — in 33,000 dimensions Euclidean distances concentrate and become uninformative. Most of those dimensions are technical noise and dropout. Reducing to 30–50 dimensions keeps the biological signal, discards noise, and makes the computation tractable.

**⭐⭐ Q: Why not cluster on UMAP coordinates?**
UMAP is a 2D visualisation that non-linearly distorts distance, invents gaps, and fragments continuous structure. It's the display layer. Cluster in high dimensions; visualise in two.

**⭐ Q: Leiden vs Louvain?**
Louvain can yield internally-disconnected communities. Leiden's refinement step guarantees well-connected communities, converges better, and is faster. Leiden is the modern default.

**Q: What does `resolution` do, mechanically?**
Leiden optimises a quality function (modularity, or CPM). `resolution` weights the penalty for having many communities — higher resolution tolerates more, smaller communities. It's a granularity dial, not a "correctness" dial.

**Q: How does `n_neighbors` interact with `resolution`?**
Both control granularity, from different angles. `n_neighbors` shapes the graph (small = local/noisy, large = smooth/global); `resolution` shapes how that graph is partitioned. Change one at a time or you can't attribute what happened.

**⭐ Q: Is clustering deterministic?**
No. Leiden is stochastic. Set `random_state`. And note that the *cluster labels* — which group is called "3" — are arbitrary and change between runs even with identical partitions. Never hard-code a cluster number across re-runs.

## Choosing resolution

**⭐⭐ Q: How do you choose the number of clusters?**
Say up front there's no single correct answer, then give your toolkit: marker interpretability (primary), bootstrap stability, clustree across resolutions, ARI between resolutions, minimum cluster size, domain knowledge. Add the pragmatic point: **over-cluster slightly and merge by markers**, because merging is easy and un-splitting is impossible.

**⭐ Q: Silhouette score picks your lowest resolution. Convinced?**
No. Silhouette is structurally biased toward fewer clusters (more clusters → smaller `b` → lower score) and assumes convex Euclidean blobs, which continuous biological manifolds are not. It's one weak signal.

**Q: What if two clusters have essentially the same markers?**
Merge them, or check whether they're separated by a technical covariate (library size, mito %, cell cycle, patient) rather than identity. Confirm with the dotplot dendrogram — clusters branching together are merge candidates.

## Marker genes and DE

**⭐⭐ Q: Your top marker has adjusted p = 1e-121. Believe it?**
As a **ranking**, yes. As **inference**, no — that's double dipping. The clusters were defined by the same expression data you're now testing, so the null was violated by construction and p-values are wildly anti-conservative. Bonus points: mention count-splitting, ClusterDE, or the TN test as principled fixes, and pseudobulk for condition-level inference.

**⭐ Q: Why Wilcoxon over t-test?**
Rank-based, so robust to the non-normal, zero-inflated, heavy-tailed reality of scRNA-seq. A t-test's normality assumption doesn't hold.

**⭐⭐ Q: Top markers are all ribosomal. What now?**
Check **effect size vs significance** (huge p, tiny log₂FC = artifact signature). Check correlation with technical covariates. Consider a proliferative *state* rather than identity. Re-rank with `RP[SL]`/`MT-` excluded. **And** consider that in some systems (leukaemia, stem cells, plasma cells) translational output is genuine biology — as it happens to be in the Caron dataset this notebook uses.

**⭐ Q: One-vs-rest has a specific blind spot. What is it?**
If a cell type is split across several clusters — six B clusters here — then testing one B cluster against "the rest" includes the other five B clusters in the reference. The shared B-cell programme cancels out, and your top markers become whatever distinguishes that *patient's clone*, not what makes it a B cell. **This is exactly what happens in this notebook.** The fix is targeted pairwise comparisons (`groups=['1'], reference='6'`) or testing at a coarser resolution/merged label.

**Q: How would you do DE between conditions (ALL vs healthy) properly?**
**Pseudobulk.** Sum raw counts per sample per cell type, then DESeq2/edgeR with sample as the unit of replication. Cells within a patient are not independent; treating 3,500 cells as n=3,500 instead of n=7 patients inflates significance enormously (pseudoreplication).

**Q: What does `standard_scale='var'` do on a dotplot, and what does it cost?**
Rescales each gene 0–1 across clusters so one highly-expressed gene doesn't saturate the colour map. Cost: you lose absolute expression magnitude.

## Data handling

**⭐ Q: Why is scRNA-seq data stored sparse?**
90–95% zeros (biology + dropout at ~10–20% capture efficiency). Sparse storage keeps only non-zeros. Densifying this notebook's 3,500 × 33,102 matrix costs ~930 MB; at 100,000 cells it's ~26 GB and your kernel dies.

**⭐ Q: What's in `.X` vs `.layers` vs `.raw`?**
`.X` = the working matrix, whatever state your pipeline has it in. `.layers` = named alternative full-size matrices (`counts`, `lognorm`, `scaled`). `.raw` = a frozen snapshot, conventionally post-lognorm/pre-HVG-subsetting, so you can still plot filtered-out genes. The trap: several functions default to `.raw` if it exists, so a mislabelled `.raw` means silently testing the wrong data.

**⭐ Q: `target_sum=1e4` vs `None`?**
`1e4` = CP10K, the community default, comparable across datasets. `None` = median library size, data-driven, avoids inflating shallow data onto an artificial scale, but not comparable to others' CP10K values.

**Q: Why `log1p` and not `log`?**
`log(0) = -inf` and the matrix is mostly zeros. `log(x+1)` maps 0 → 0.

**Q: Gene symbols or Ensembl IDs as your index?**
Ensembl IDs are stable and unique; symbols are neither (they get renamed, and several IDs map to one symbol). Index on Ensembl, carry symbols in `.var['Symbol']` for display. If you must use symbols, call `.var_names_make_unique()` and know you've made an arbitrary choice.

## Batch correction

**⭐ Q: How do you know batch correction worked?**
Visually, batches intermix within clusters while distinct cell types stay distinct. Quantitatively, kBET, LISI/iLISI, or scIB's composite metrics. The key point: you want batches mixed **and** cell types preserved. Over-correction merges genuinely distinct populations and looks like success on a mixing metric.

**⭐⭐ Q: A cluster contains cells from one patient only. Batch effect or biology?**
"Depends on the biology, and here's how I'd tell." Batch: driven by technical covariates, affects all cell types, markers are housekeeping. Biology: coherent patient-specific programme, known clonal disease, other populations mix fine. **In cancer, per-patient tumour clusters are expected**, and correcting them away destroys your signal.

---

# 🔴 What experienced people watch for in a notebook like this

Ranked by how much damage they cause.

| # | Watch for | Why it's dangerous | Guard |
|---|---|---|---|
| 1 | **`use_rep` on `sc.pp.neighbors`** | Omit it and you silently cluster the uncorrected space. Every downstream result is wrong and everything still *looks* fine | Always pass it explicitly; verify with `ad.uns['neighbors']['params']` |
| 2 | **Silent in-place overwrites** — `sc.tl.umap` clobbering `X_umap` (cell 20) | No error, no warning, and later comparisons become meaningless | `key_added=`, or `.copy()` first. Print `ad` before and after every `tl`/`pp` call |
| 3 | **Double normalisation** | `normalize_total` has no guard. Result is garbage that looks plausible | Check `max()` and integer-ness of `.X` before normalising |
| 4 | **Raw counts destroyed** | Many methods *need* counts. Recovering them means re-running the pipeline | `ad.layers["counts"] = ad.X.copy()` before any `pp` call. Non-negotiable |
| 5 | **Believing marker p-values** | Double dipping makes them invalid. People build papers on this | Rank with them, never infer with them. Pseudobulk for real inference |
| 6 | **Significance without effect size** | p=1e-241 with log₂FC=1.6 (cell 45) is an artifact wearing a lab coat | Always read `logfoldchanges` and `pts` alongside `pvals_adj` |
| 7 | **One marker = one cell type** | Dropout means absence proves nothing; ambient RNA means presence proves little | Panels of 3–5 concordant markers, confirmed on violins |
| 8 | **`.raw` / `use_raw` confusion** | Different functions have different defaults. Silently tests the wrong matrix | Check `ad.raw`; pass `use_raw` explicitly, always |
| 9 | **`.map()` producing NaN** | Missing dict key → silent `NaN` cell types → grey dots you'll rationalise | `assert` no NaN after every mapping |
| 10 | **Hard-coded cluster numbers** | Leiden numbering is arbitrary and shifts between runs/versions | Annotate to a named column immediately, keep `(cN)` in the label, pin `random_state` |
| 11 | **Nothing ever saved** | `prefix_output` defined and unused. Kernel dies, work dies | `ad.write(...)` after every expensive step |
| 12 | **`.todense()` on a real matrix** | ~930 MB here; ~26 GB at realistic scale. Instant kernel death | `.toarray()` on a small slice only |
| 13 | **Ignoring `FutureWarning`** | Cell 16's warning means your cluster numbers will change on a scanpy upgrade | Pin the future behaviour explicitly, pin your environment |
| 14 | **Out-of-order execution** | Cell 1 shows `execution_count 8` — it was re-run later. Notebooks lie about what ran when | `Restart & Run All` before trusting any result you'll report |
| 15 | **No doublet detection anywhere** | Doublets form convincing fake "intermediate" clusters | `scrublet`/`scDblFinder` upstream; be suspicious of clusters co-expressing markers of two lineages |
| 16 | **Unlabelled ambiguous clusters** | Forcing names on everything turns doublet nests into novel cell types | `"Unknown (cN)"` is a legitimate, honest label |

---

# 💡 Tips, tricks, and snippets worth keeping

## The five-line safety preamble — paste at the top of every notebook

```python
import scanpy as sc
sc.settings.verbosity = 3            # tell me what you're actually doing
sc.settings.set_figure_params(dpi=100, facecolor='white')
sc.logging.print_header()            # ⭐ prints all package versions — reproducibility
```
`print_header()` output pasted into your methods section is most of what a reviewer needs.

## The load ritual — do this every time, without thinking

```python
ad = sc.read_h5ad(path)
print(ad)                               # read this properly, don't skim
print("raw:", ad.raw)                   # None, or something?
print("layers:", list(ad.layers))
print("obsm:", list(ad.obsm))
ad.layers["counts"] = ad.X.copy()       # ⭐ the Horcrux
```

## Is this matrix normalised? — settle it in three lines

```python
import numpy as np
x = ad.X[:50].toarray() if hasattr(ad.X, "toarray") else ad.X[:50]
print(f"max={x.max():.2f}  integers={np.allclose(x, np.round(x))}")
# integers, max in the hundreds/thousands → raw counts
# floats, max ~5-10                       → lognormalised. Don't touch.
```

## Save aggressively

```python
import os
os.makedirs(prefix_output, exist_ok=True)   # makedirs, not mkdirs
ad.write(f"{prefix_output}/clustered.h5ad")
```
Clustering a large dataset takes real time. Losing it to a kernel restart is entirely avoidable, and everyone learns this the hard way once.

## The right way to sweep resolutions

```python
for res in [0.2, 0.4, 0.6, 0.8, 1.0, 1.4]:
    key = f"leiden_res{res:.1f}"          # :.1f → always "1.0", never "1"
    sc.tl.leiden(ad, resolution=res, key_added=key,
                 flavor="igraph", n_iterations=2, directed=False, random_state=0)
    vc = ad.obs[key].value_counts()
    print(f"{key}: {len(vc)} clusters, smallest={vc.min()}, largest={vc.max()}")
```

## How much did the partition actually change?

```python
from sklearn.metrics import adjusted_rand_score as ari
keys = [f"leiden_res{r:.1f}" for r in [0.2,0.4,0.6,0.8,1.0,1.4]]
for a, b in zip(keys, keys[1:]):
    print(f"{a} → {b}: ARI = {ari(ad.obs[a], ad.obs[b]):.3f}")
```
**ARI near 1** = nothing meaningfully changed, you're just splitting hairs. **A sharp drop** = real restructuring; investigate what split. Far more informative than the cluster count alone.

## Marker table you can actually filter

```python
sc.tl.rank_genes_groups(ad1, groupby="label", method="wilcoxon",
                        pts=True, tie_correct=True,      # ⭐ both worth having
                        key_added="rgg")

df = sc.get.rank_genes_groups_df(ad1, group=None, key="rgg")   # group=None → ALL clusters
clean = df[(df.pvals_adj < 0.01) &
           (df.logfoldchanges > 1) &
           (df.pct_nz_group > 0.5) &                            # in >50% of cluster cells
           (~df.names.str.match(r'^(RP[SL]|MT-|MRP[SL])'))]     # drop housekeeping
clean.to_csv(f"{prefix_output}/markers.csv", index=False)
```
**`pct_nz_group > 0.5` is the filter that does the most work.** It removes every "marker" that's really a handful of outlier cells — the most common source of bogus annotations.

## Ask the specific question, not the generic one

One-vs-rest can't tell you why clusters 6 and 7 differ, because each is in the other's reference set. Ask directly:

```python
sc.tl.rank_genes_groups(ad1, groupby="label", groups=["7"], reference="6",
                        method="wilcoxon", key_added="c7_vs_c6")
sc.get.rank_genes_groups_df(ad1, group="7", key="c7_vs_c6").head(20)
```
**This is the single best answer to cell 44's exercise**, and it's more informative than either cell 45 or cell 46.

## Check whether your clusters are really technical

```python
qc = ['sum','detected','subsets_Mito_percent']
sc.pl.violin(ad1, qc, groupby='label', rotation=90)
print(ad1.obs.groupby('label')[qc].median())
```
Run this **before** annotating. A cluster with conspicuously low `detected` and high mito % is dying cells, not a cell type. This is the diagnostic that would have explained cluster 6's ribosomal signature.

## Cell composition per sample

```python
pd.crosstab(ad1.obs['SampleName'], ad1.obs['cellType'], normalize='index').round(3)
```
Instantly shows which clusters are patient-specific (leukaemic clones) and which are shared (normal immune cells). One line, and it's the table that explains this entire dataset. It's also the starting point for notebook 07 (differential abundance).

## Cluster relationships

```python
sc.tl.dendrogram(ad1, groupby='label', use_rep='X_corrected')
sc.pl.dendrogram(ad1, groupby='label')
```
Run it yourself rather than letting the dotplot silently do it with defaults (cell 40's warning). Clusters that branch together are your merge candidates.

## A marker panel worth memorising for PBMC/bone marrow

| Population | Markers |
|---|---|
| **T cells (pan)** | `CD3D`, `CD3E`, `CD2`, `TRAC` |
| CD4 T | `IL7R`, `CD4`, `CCR7` (naive) |
| CD8 T | `CD8A`, `CD8B` |
| **NK** | `NKG7`, `GNLY`, `KLRD1`, `NCAM1`, `FCGR3A` — and **`CD3D`-negative** |
| **B cells (pan)** | `CD79A`, `CD79B`, `CD19` |
| Mature B | `MS4A1` (=CD20) |
| Plasma | `JCHAIN`, `MZB1`, `SDC1` |
| **Monocytes** | `LYZ`, `CST3`, `S100A8/9`, `TYROBP` |
| CD14 mono | `CD14` |
| CD16 mono | `FCGR3A`, `MS4A7` |
| **Dendritic** | `FCER1A`, `CLEC9A`, `LILRA4` (pDC) |
| **Erythroid** | `HBA1`, `HBA2`, `HBB`, `GYPA`, `ALAS2` |
| **Platelets** | `PPBP`, `PF4` |
| **Progenitors** | `CD34`, `KIT`, `MPO` |
| **Cycling (state!)** | `MKI67`, `TOP2A`, `PCNA` — a *state*, not an identity |

**Reference databases when your panel isn't enough:** CellMarker, PanglaoDB, Azimuth, the Human Cell Atlas. And automated tools — `CellTypist`, `SingleR`, `scANVI` — which are excellent **first drafts** and should never be the final word. Always eyeball the markers behind an automated call.

## Two habits that outlast this course

**Print the object before and after every mutation.** It's free, and it catches overwrites, typos and double-normalisation within seconds rather than notebooks later.

**Restart & Run All before you report anything.** Cell 1 in this notebook has `execution_count 8` — it was re-run after later cells. Notebook state is a lie by default, and out-of-order execution is the cause of a genuinely large share of irreproducible results.

---

# The whole notebook in twelve lines

```python
ad = sc.read_h5ad("Caron_batch_corrected.500.h5ad")   # 3500 cells × 33102 genes
ad.layers["counts"] = ad.X.copy()                      # ⭐ the line the notebook forgot
sc.pp.neighbors(ad, use_rep="X_corrected",             # the KNN graph — the foundation
                n_neighbors=10, n_pcs=40)
sc.tl.leiden(ad, resolution=0.3, key_added="leiden_res0.3",
             flavor="igraph", n_iterations=2, directed=False, random_state=0)
sc.tl.umap(ad, key_added="X_umap_scanpy")              # ⭐ don't clobber X_umap
sc.pp.normalize_total(ad, target_sum=None); sc.pp.log1p(ad)
sc.tl.rank_genes_groups(ad, groupby="leiden_res0.3", method="wilcoxon", pts=True)
sc.pl.rank_genes_groups_dotplot(ad, n_genes=5, standard_scale="var")
ad.obs["cellType"] = ad.obs["leiden_res0.3"].map(cell_annotation)   # ⭐ assert no NaN
ad.write("results/04/annotated.h5ad")                  # ⭐ actually save it
```

**Graph → communities → markers → names.** Everything else is care and paranoia.

---

Three things in this notebook are worth fixing in your own copy, and I'd rank them: (1) `ad.layers["counts"] = ad.X.copy()` right after loading, (2) `sc.tl.umap(ad, key_added=...)` in cell 20 so the comparison in cell 24 means what it claims, and (3) an `assert` on the `.map()` in cell 50.

Want me to apply those three fixes to the notebook, or add these explanations as markdown cells inline so they sit next to the code? I can also put this guide on a page you can keep open next to the notebook rather than scrolling terminal output.
