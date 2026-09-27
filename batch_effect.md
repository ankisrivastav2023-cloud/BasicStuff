# Notebook 05: Batch correction and data integration, explained for your interview

## The idea in one line
You have two samples of cells (PBMMC_1a and PBMMC_1b, which are healthy paediatric bone marrow). They were processed separately, so there are technical differences between them (a **batch effect**). The notebook merges them and uses **Harmony** to remove those technical differences while keeping the real biology.

> 🧙 **HP analogy:** Hogwarts and Durmstrang students arrive for the Triwizard Tournament. A Seeker is a Seeker in either school, but the Durmstrang ones wear fur cloaks. If you sort students by how they look, you get "Hogwarts kids" and "Durmstrang kids" instead of "Seekers" and "Beaters". The cloak is the **batch effect**, and batch correction takes the cloaks off.
>
> 🌙 **ACOTAR analogy:** An Illyrian warrior is an Illyrian warrior whether they're in Velaris or visiting the Spring Court. Each court puts its own **glamour** on them, though. Integration means seeing through the glamour, the way Feyre sees Velaris for what it really is, to the true identity underneath.

---

## Step by step

### 1. Setup (cells 2–6)
- Silences Numba warnings, creates an output folder, and imports `scanpy`, `anndata`, `harmonypy`, and the usual tools.
- `np.random.seed(10)` makes the results **reproducible**.
- 💡 **Pro tip:** UMAP and Harmony both involve randomness. Always say you fix seeds so the results can be reproduced.

### 2. Load and combine the two samples (cells 9–19)
```python
ad.concat([adata_1a, adata_1b], axis=0, join='inner', merge='same',
          label='batch', keys=['batch_1','batch_2'], index_unique='_')
```
| Parameter | What it does | Analogy |
|---|---|---|
| `axis=0` | Stacks **cells** (rows) | Combines both schools' student registers into one |
| `join='inner'` | Keeps only **genes present in both** samples | Keeps only subjects both schools teach |
| `merge='same'` | Keeps gene annotations only where they match in both samples | |
| `label='batch'` | Adds a `.obs['batch']` column recording where each cell came from | A school badge on every student |
| `index_unique='_'` | Adds a suffix to cell barcodes (`AAAC..._batch_1`) | Stops two students named "Harry" from clashing |

- `.var.astype(str)` is a workaround so `concat` doesn't fail on mixed column types.
- `adata_combined_counts = adata_combined.copy()` keeps a backup of the **raw counts**.

### 3. Look at the data before correcting it (cells 22–27)
The standard pipeline with **no** batch awareness: store counts in a layer → `normalize_total` → `log1p` → HVGs → PCA → neighbours → UMAP, then colour the UMAP by batch.
- **Result:** the two batches separate on the UMAP, which means there's a batch effect.
- 💡 **Pro tip:** Always plot the uncorrected data first. If the batches already overlap, you may not need integration at all.

### 4. Batch-aware HVG selection (cells 29–32)
```python
sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000, flavor="cell_ranger", batch_key='batch')
```
- HVGs are computed **within each batch** and then combined. This stops genes that vary only because of the batch from being picked.
- It adds two columns:
  - `highly_variable_nbatches`: how many batches the gene is an HVG in
  - `highly_variable_intersection`: whether the gene is an HVG in **every** batch
- 🧙 It's like picking prefects from each House separately, so Slytherin doesn't take every spot just because it's the loudest House.

### 5. Harmony integration (cells 34–40)
```python
adata_combined.obsm['X_umapraw'] = adata_combined.obsm['X_umap']   # back up the old UMAP
sc.external.pp.harmony_integrate(adata_combined, key='batch', basis='X_pca')
sc.pp.neighbors(adata_combined, use_rep='X_pca_harmony', key_added='harmony')
sc.tl.umap(adata_combined, neighbors_key='harmony')
```

**How Harmony works (be ready to explain this):**
1. It starts from the **PCA embedding**, not the gene counts.
2. It runs **soft k-means clustering**, with a **diversity penalty (θ, theta)** that discourages clusters made mostly of one batch.
3. Within each cluster, it calculates a batch-specific correction and shifts the cells.
4. It repeats steps 2–3 until the result stops changing.
5. The output is a new embedding, **`X_pca_harmony`**. The count matrices are never touched.

> 🌙 Harmony works like the Night Court's Inner Circle. Rhys pulls people from different backgrounds (Illyrians, a Cursebreaker, an ancient being) into one group based on who they really are, not where they came from. θ is how strongly he insists on mixing.
>
> 🧙 Or think of it as the Sorting Hat being told: "Don't sort by which train carriage they arrived in."

Then neighbours and UMAP are recalculated from `X_pca_harmony`. The batches should now overlap, while the cell types stay separate.

### 6. Save (cell 41)
The result is written to `.h5ad`.

---

## Exercise answers (interviewers love these)

**Q1. What changes after each method?**

| Element | Harmony | BBKNN | MNN_correct |
|---|---|---|---|
| Raw counts | ❌ | ❌ | ❌ |
| Log-normalised matrix | ❌ | ❌ | ✅ (corrected) |
| Scaled matrix | ❌ | ❌ | ✅ |
| HVG list | ❌ | ❌ | ❌ (usually) |
| Original `X_pca` | ❌ (a **new** `X_pca_harmony` is added) | ❌ | ✅ (recomputed from the corrected matrix) |
| Embedding used for neighbours | ✅ | ❌ (same PCA) | ✅ |
| Neighbour graph | ✅ | ✅ (**this is all BBKNN changes**) | ✅ |
| UMAP | ✅ | ✅ | ✅ |

**Summary:** Harmony corrects the **embedding**, BBKNN corrects the **neighbour graph**, and MNN corrects the **expression matrix**.

**Q2. Can you use corrected pseudo-counts (MNN, Seurat CCA/RPCA) for differential expression?**
**No.** Here's why:
- They aren't real counts anymore, so they break the negative binomial / Poisson assumptions that DE models rely on.
- Each cell's corrected value borrows information from other cells, so the cells are no longer independent. That gives inflated p-values and false positives.
- Correction can remove real biology (over-correction) or create artificial signal.
- ✅ **What to do instead:** run DE on **raw counts** and put batch in the model as a covariate, e.g. pseudobulk + DESeq2/edgeR with `~ batch + condition`.

> 🧙 It's like grading a student on an essay that Polyjuice-Hermione wrote for them. The result looks convincing, but it isn't really theirs.

---

## 💡 Pro tips that will make you stand out

1. **Spot the subtle problem in this notebook:** PCA (cell 23) was run **before** the batch-aware HVG selection (cell 29) and was **never rerun**. Harmony therefore used PCs built from the old, batch-unaware HVGs. Ideally you rerun `sc.tl.pca` after the batch-aware HVG step. Mentioning this shows you understand the order of the pipeline.
2. **Inconsistent PCs:** the uncorrected neighbours use `n_pcs=12`, but the Harmony neighbours use `n_pcs=40`. For a fair comparison, use the same number.
3. **The Q1 question mentions a scaled matrix, but `sc.pp.scale` is never run** in this notebook.
4. **Over-correction is the real danger.** If the batches overlap perfectly *and* distinct cell types merge into one blob, you've corrected away biology. 🌙 That's like the Cauldron mixing everything into one: nobody is distinct anymore.
5. **Don't just trust how the UMAP looks. Measure it:** **iLISI** (are batches mixed?), **cLISI** (are cell types still separate?), **kBET**, **silhouette score**, and the **scIB** benchmark suite.
6. **Choosing a method:**
   - Harmony: fast, and good for simple batches.
   - BBKNN: very fast.
   - scVI/scANVI: deep learning, best for complex or large atlases.
   - Scanorama.
   - Seurat CCA/RPCA (R).
7. **Don't correct away the biology you're studying.** If batch is **confounded** with condition (for example, all disease samples in batch 1), no method can separate the two. That's an experimental design problem. 🧙 Even Dumbledore can't fix that.

---

## Likely interview questions (with short answers)

- **What is a batch effect?** Technical variation from being processed in different runs, days, chemistries or lanes, not from biology.
- **How do you detect one?** Colour the UMAP/PCA by batch. If cells group by batch instead of by cell type, there's a problem. You can also compute LISI or kBET.
- **Why use batch-aware HVGs?** So genes that vary only because of the batch don't drive the PCA.
- **What does Harmony output?** A corrected PCA embedding (`X_pca_harmony`). The counts stay the same.
- **Harmony vs BBKNN vs MNN?** Embedding vs graph vs expression matrix (see the table).
- **What does θ (theta) in Harmony do?** It sets the strength of the diversity penalty. Higher θ gives more mixing, with a higher risk of over-correction.
- **Can you run DE on integrated data?** Not on corrected values. Use raw counts with batch as a covariate.
- **What does `join='inner'` vs `'outer'` do?** Inner keeps only shared genes. Outer keeps all genes and fills the missing ones with zeros or NaN.
- **Why keep `layers["counts"]`?** Normalisation overwrites `.X`, and later steps (DE, scVI) need the raw counts.
- **How would you tell integration worked?** Batches mix within each cell type, known marker genes still separate the cell types, and LISI/ASW scores look good.
- **When should you *not* integrate?** When the batches already overlap, or when batch is confounded with the biology you care about.

Good luck tomorrow! 🪄🌙



# Notebook 05, explained from zero

For each step you'll get **what** it does in plain words, **why** it's done, the **code**, the **terms to say out loud**, and a one-line analogy.

---

## Part 0: Five basics to know first

**1. What is scRNA-seq data?**
Single-cell RNA sequencing measures which genes are switched on in **each individual cell**. The result is a big table:
- **Rows = cells** (thousands of them)
- **Columns = genes** (about 20,000–30,000)
- **Each value = a count**: how many RNA molecules of that gene were detected in that cell. These are called **UMI counts** (Unique Molecular Identifiers).

Most values are **0**, so the matrix is called **sparse**. This happens because each cell has only a small amount of RNA and capture is inefficient. The missing zeros are often called **dropouts**.

**2. What is AnnData?**
AnnData is the Python object that holds all of this. Think of it as Hermione's beaded bag: one container with labelled compartments.

| Slot | What's inside | Example |
|---|---|---|
| `.X` | The main matrix (cells × genes) | counts or normalised values |
| `.obs` | Table of **cell** information (one row per cell) | batch, cluster, QC metrics |
| `.var` | Table of **gene** information (one row per gene) | gene name, is it highly variable? |
| `.obsm` | Low-dimensional **embeddings** per cell | `X_pca`, `X_umap`, `X_pca_harmony` |
| `.obsp` | Cell-by-cell **graphs** | neighbour distances and connectivities |
| `.layers` | Alternative versions of `.X` | `counts`, `logcounts` |
| `.uns` | Unstructured extras | plot colours, parameters used |

`.h5ad` is the file format AnnData is saved in.

**3. What is a batch effect?**
A **batch** is a group of cells processed together, for example on the same day, in the same 10x Chromium run, or on the same sequencing lane. A **batch effect** is **technical (non-biological) variation** between batches. Causes include differences in sequencing depth, reagent lots, capture efficiency, or ambient RNA.
- **The problem:** identical cell types from two batches can look different, so the cells group by **batch** instead of by **cell type**.

**4. What is data integration / batch correction?**
Computational methods that **remove the technical batch differences while preserving the real biological differences**. There are always two goals, and interviewers want to hear both:
- ✅ **Batch mixing:** cells of the same type from different batches overlap.
- ✅ **Biological conservation:** different cell types stay separate.

**5. The dataset**
`PBMMC_1a` and `PBMMC_1b` are **Paediatric Bone Marrow Mononuclear Cells**: healthy control samples from the Caron et al. 2020 childhood leukaemia (ALL) dataset used throughout this course. `_dimRed` in the filename means each file has already gone through QC, normalisation and dimensionality reduction in an earlier notebook.

---

## Step 1: Setup (cells 2–6)

```python
from numba.core.errors import NumbaDeprecationWarning, ...
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
```
- **What:** hides warning messages from **Numba**, a library scanpy uses to speed up code. These warnings are harmless noise.

```python
prefix_output = ".../05 Batch Correction"
os.makedirs(prefix_output, exist_ok=True)
```
- **What:** creates a folder to save results in. `exist_ok=True` means "don't crash if the folder already exists".

```python
import scanpy as sc      # main single-cell analysis toolkit
import anndata as ad     # the data container
import numpy as np       # maths on arrays
import pandas as pd      # tables
import matplotlib.pyplot as plt  # plotting
import harmonypy         # the Harmony integration algorithm
```

```python
pd.set_option('display.max_columns', 40)   # show up to 40 table columns when printing
seed = 10
np.random.seed(seed)
```
- **Why:** PCA solvers, UMAP and Harmony all use **random numbers** internally. Fixing the **random seed** means you get the same result every time you run the notebook. This is **reproducibility**.
- 🧙 It's like a Time-Turner: you can replay the exact same events.

**🗣️ Say:** "I set a random seed for reproducibility, because UMAP and Harmony are stochastic."

---

## Step 2: Load the two samples (cells 9–12)

```python
adata_1a = sc.read_h5ad(".../PBMMC_1a_dimRed.h5ad")
adata_1b = sc.read_h5ad(".../PBMMC_1b_dimRed.h5ad")
adata_1a.shape   # (number of cells, number of genes)
```
- **What:** loads each sample as its own AnnData object and checks its size.
- `.shape` returns `(n_obs, n_vars)`, i.e. (cells, genes).

---

## Step 3: Merge them into one object (cells 13–19)

```python
adata_1a.var = adata_1a.var.astype(str)
adata_1b.var = adata_1b.var.astype(str)
```
- **What:** converts every column of the gene table to text.
- **Why:** it's a **technical workaround**. `concat` can fail when the same column has different data types in the two objects. It isn't a scientific step.

```python
adata_combined = ad.concat([adata_1a, adata_1b],
                           axis=0, join='inner', merge='same',
                           label='batch', keys=['batch_1','batch_2'],
                           index_unique='_')
```
Parameter by parameter:
- `axis=0`: stack along the **cell** axis, so the cells from 1b go underneath the cells from 1a.
- `join='inner'`: keep only the **genes found in both** objects (the intersection). `'outer'` would keep all genes and fill the gaps.
- `merge='same'`: in `.var`, keep only the annotation values that are **identical** in both objects.
- `label='batch'`: add a new `.obs` column called `batch` recording which sample each cell came from. **This is the column integration will use.**
- `keys=['batch_1','batch_2']`: the values written into that column.
- `index_unique='_'`: add the key to each cell barcode, e.g. `AAACCTGA-1_batch_1`.
  - **Why:** the same barcode sequence can appear in both samples, and every cell name must be **unique**.

```python
adata_1a.n_obs + adata_1b.n_obs   # sanity check: combined cell count = sum of both
adata_combined_counts = adata_combined.copy()   # backup of raw counts
```
- 🧙 Two House registers are merged into one Hogwarts register. Every student gets a House badge (`batch`), and name clashes are fixed by adding the House name.

**🗣️ Terms:** *concatenation, inner join, cell barcode, obs/var, sanity check.*

---

## Step 4: Analyse without correction first, as the baseline (cells 22–27)

This runs the **standard scanpy pipeline** while ignoring batch, so you can see whether a batch effect exists.

### 4a. Keep the raw counts
```python
adata_combined.layers["counts"] = adata_combined.X.copy()
```
- **Why:** the next step overwrites `.X`. Raw counts are needed later for **differential expression** and for methods like **scVI** that model counts directly.

### 4b. Normalisation
```python
sc.pp.normalize_total(adata_combined)
```
- **The problem it fixes:** cells are sequenced to different depths. One cell might have 2,000 total counts and another 10,000. The second cell's genes all look "higher" simply because more of its RNA was sequenced, not because of biology.
- **What it does:** **library-size normalisation**. It scales each cell so that all cells have the same total count. With no `target_sum` given, scanpy uses the **median total count** across cells.
- **Terms:** *library size, sequencing depth, size factor.*

### 4c. Log transform
```python
sc.pp.log1p(adata_combined)
adata_combined.layers["logcounts"] = adata_combined.X.copy()
```
- **What:** replaces each value x with **log(x + 1)**. The "+1" avoids log(0), which is undefined.
- **Why:**
  - Gene expression is highly **skewed**: a few genes have huge counts. The log compresses that range.
  - It **stabilises variance**, so highly expressed genes don't dominate the analysis.
  - Differences become **fold changes**, which is how biologists think about expression.
- **Terms:** *log-normalised expression, variance stabilisation.*

### 4d. Highly Variable Genes (HVGs)
```python
sc.pp.highly_variable_genes(adata_combined, flavor='seurat')
```
- **What:** finds the genes whose expression **varies most between cells**, relative to their average expression.
- **Why:** most genes are either constant (e.g. housekeeping genes) or pure noise. Genes that vary between cells are the ones that **distinguish cell types**. Keeping only them (usually 1,000–3,000) reduces noise and speeds up computation.
- **How (seurat flavor):** for each gene it computes the **dispersion** (variance ÷ mean). Genes are grouped into bins by mean expression, dispersion is **z-scored within each bin**, and the top genes are marked with `.var['highly_variable'] = True`.
- **Term:** *feature selection.*
- 🧙 Choosing Quidditch team members by what makes them *different* (speed, agility), not by what everyone shares (they all have two legs).

### 4e. PCA (Principal Component Analysis)
```python
sc.tl.pca(adata_combined)
sc.pl.pca_variance_ratio(adata_combined, n_pcs=30)
```
- **What:** a **linear dimensionality reduction**. It compresses about 2,000 genes into around 50 new axes called **principal components (PCs)**. PC1 captures the most variation, PC2 the next most, and so on. By default scanpy uses only the HVGs.
- **Why:** removes noise and redundancy (many genes move together), and makes the next steps much faster.
- **Stored in:** `.obsm['X_pca']` (cells × 50).
- **Variance ratio plot ("elbow plot"):** shows how much variance each PC explains. Choose the number of PCs where the curve **flattens (the elbow)**. PCs after that point are mostly noise.
- 🌙 Like the Suriel: it takes an overwhelming amount of information and tells you only the few things that matter most.

### 4f. Neighbour graph
```python
sc.pp.neighbors(adata_combined, n_pcs=12)
```
- **What:** for each cell, finds its **k nearest neighbours** (default k = 15) in PCA space, using the first 12 PCs chosen from the elbow plot. The result is a **kNN graph**: cells are dots connected to their most similar cells.
- **Stored in:** `.obsp['distances']` and `.obsp['connectivities']`.
- **Why:** this graph is the foundation for **clustering** (Leiden/Louvain) and for **UMAP**.

### 4g. UMAP
```python
sc.tl.umap(adata_combined)
sc.pl.umap(adata_combined, color=['batch'])
```
- **What:** **UMAP (Uniform Manifold Approximation and Projection)** is a **non-linear** method that turns the neighbour graph into a **2D picture** for visualisation.
- ⚠️ **Scientific caution:** UMAP preserves **local** structure well (which cells are near each other). **Distances between clusters and cluster sizes are not meaningful.** Use UMAP only for visualisation, never for quantitative conclusions.
- **Stored in:** `.obsm['X_umap']`.

### 4h. Result
When the UMAP is coloured by batch, **the two batches separate**, so there is a **batch effect**.
- 🧙 On the Marauder's Map, Gryffindors cluster on one side and Slytherins on the other, even when they're taking the same class.

**🗣️ Say:** "I always visualise the unintegrated data first to confirm a batch effect exists before correcting, because unnecessary correction can remove real biology."

---

## Step 5: Batch-aware HVG selection (cells 29–32)

```python
sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000,
                            flavor="cell_ranger", batch_key='batch')
```
- **The problem:** when HVGs are selected on the merged data, genes that differ **between batches** (technical differences) look "highly variable" and get selected. They then drive the PCA, which makes the batch effect worse.
- **The fix, `batch_key='batch'`:**
  1. HVGs are calculated **separately within each batch**.
  2. Genes are ranked first by **how many batches they're highly variable in**, then by their average normalised dispersion.
  3. The top 2,000 are selected.
- **`flavor="cell_ranger"`:** a variant of the dispersion method that normalises dispersion using the median and **median absolute deviation** within expression bins. It expects **log-normalised** data.
- **New columns in `.var`:**
  - `highly_variable_nbatches`: the number of batches in which the gene is an HVG (0, 1 or 2 here).
  - `highly_variable_intersection`: `True` if the gene is an HVG in **all** batches.
- The **bar plot** shows how many genes are HVGs in 0, 1 or 2 batches. Genes that are HVGs in both batches are **robust**: they're variable because of biology, not batch.
- 🌙 Picking the Inner Circle based on qualities valued in **every** court, not just the Night Court.

**🗣️ Terms:** *batch-aware feature selection, robust HVGs.*

---

## Step 6: Integration with Harmony (cells 34–40)

### 6a. Back up the old UMAP
```python
adata_combined.obsm['X_umapraw'] = adata_combined.obsm['X_umap']
```
- **Why:** the next `sc.tl.umap` call will overwrite `X_umap`. Saving it under a new name lets you compare before and after.

### 6b. Run Harmony
```python
sc.external.pp.harmony_integrate(adata_combined, key='batch', basis='X_pca')
```
- `sc.external` holds **third-party tools wrapped by scanpy**. This one calls the `harmonypy` package.
- `key='batch'`: which `.obs` column defines the batches.
- `basis='X_pca'`: Harmony works on the **PCA embedding**, not on gene counts.

**How Harmony works (Korsunsky et al., *Nature Methods*, 2019):**
1. It starts from the cells' PCA coordinates.
2. It groups cells with **soft k-means clustering**. "Soft" means each cell has a *probability* of belonging to each cluster rather than a single assignment.
3. A **diversity penalty (parameter θ, theta)** discourages clusters made mostly of one batch, pushing clusters to contain cells from all batches.
4. Inside each cluster, it computes a **batch-specific linear correction** (a mixture-of-experts linear model) and shifts cells to remove the batch offset.
5. Steps 2–4 **repeat** until the result **converges** (stops changing).
6. The output is **`.obsm['X_pca_harmony']`**, a corrected embedding.

**Key facts to remember:**
- ✅ Harmony changes **only the embedding**. **The expression matrix (`.X`, layers) is not modified.**
- ✅ The original `X_pca` is kept, and a new `X_pca_harmony` is added.
- It is fast and scales to millions of cells.

🌙 Harmony is like Rhysand's Inner Circle meetings. Members come from different origins (Illyrian, Cursebreaker, Feyre from the human lands), but they're repeatedly grouped by *role and nature*, and nobody's origin gets to define their group. θ is how firmly he insists that every group is mixed.

### 6c. Rebuild the neighbours and UMAP on the corrected embedding
```python
sc.pp.neighbors(adata_combined, use_rep='X_pca_harmony', key_added='harmony',
                n_neighbors=15, n_pcs=40)
sc.tl.umap(adata_combined, neighbors_key='harmony')
```
- `use_rep='X_pca_harmony'`: build the kNN graph from the **corrected** coordinates.
- `key_added='harmony'`: store this graph under a **separate name** (`.obsp['harmony_connectivities']`, etc.) so the original graph isn't overwritten.
- `neighbors_key='harmony'`: tells UMAP which graph to use.

### 6d. Plot and save
```python
sc.pl.umap(adata_combined, color=['batch'])
adata_combined.obsm['X_umapharmony'] = adata_combined.obsm['X_umap']
```
- The batches should now **overlap**, meaning they're well mixed.
- The Harmony UMAP is saved under its own name.

---

## Step 7: Save (cell 41)
```python
adata_combined.write_h5ad(f'{prefix_output}/05_integrated_PBMMC_1a_1b_dimRed.h5ad')
```
- Saves everything so the next notebooks (clustering, DE) can load it.
- The `f'...'` is an **f-string**: Python fills `{prefix_output}` in with the folder path.

---

## Step 8: Judging whether integration worked
A good result has:
1. **Batches mixed** within each cell type.
2. **Cell types still separate**, with marker genes still distinct.

⚠️ **Over-correction:** batches are perfectly mixed, but different cell types have been merged together. That means real biology was removed.

**Quantitative metrics** (worth name-dropping):
- **iLISI** (integration Local Inverse Simpson's Index): higher means better batch mixing.
- **cLISI**: lower means cell types are kept separate.
- **kBET** (k-nearest-neighbour Batch Effect Test)
- **ASW** (Average Silhouette Width), computed for batch and for cell type.
- The **scIB** benchmarking package (Luecken et al., 2022) combines these.

---

## Pro tips the interviewer may probe

1. **Order-of-steps issue in this notebook:** PCA was computed **before** the batch-aware HVG selection and wasn't rerun. So Harmony corrected PCs built from the non-batch-aware HVGs. Best practice: **batch-aware HVGs → PCA → Harmony**.
2. **Correcting the embedding vs correcting the matrix.** Know the three types:
   - **Embedding:** Harmony, scVI (latent space)
   - **Graph:** BBKNN (Batch Balanced kNN)
   - **Expression matrix:** MNN correct, Seurat CCA/RPCA, ComBat
3. **Don't use batch-corrected expression values for DE.** They aren't true counts, cells become non-independent, and false positives go up. Use **raw counts with batch as a covariate**, ideally with **pseudobulk** (summing counts per sample and cell type) plus DESeq2/edgeR.
4. **Confounding:** if batch = condition (e.g. all patients in batch 1, all controls in batch 2), **no method can separate the technical signal from the biological one**. It's an experimental design problem, and the fix is balanced designs.
5. **Batch mixing vs biological conservation** is the core trade-off. Always mention both goals.

---

## Terms checklist (say these naturally)
batch effect · technical vs biological variation · sequencing depth / library size · UMI counts · sparsity / dropouts · library-size normalisation · log1p transform · variance stabilisation · highly variable genes / feature selection · dispersion · batch-aware HVG · PCA / principal components / elbow plot · embedding · kNN graph · UMAP (local structure only) · data integration · Harmony · soft k-means · diversity penalty θ · convergence · over-correction · confounding · iLISI / cLISI / kBET / ASW / scIB · pseudobulk · covariate

## A 30-second summary you can say out loud
> "I merged two bone marrow samples with an inner join on genes and labelled each cell by batch. I ran the standard pipeline (library-size normalisation, log1p, HVGs, PCA, kNN graph, UMAP) and saw the cells separating by batch, which confirmed a batch effect. I then selected HVGs in a batch-aware way and ran Harmony on the PCA embedding. Harmony iteratively clusters cells with a diversity penalty and applies cluster-specific linear corrections, producing a corrected embedding without changing the counts. I rebuilt the neighbour graph and UMAP on that embedding and confirmed the batches mixed while the cell types stayed distinct. For downstream differential expression I'd use raw counts with batch as a covariate, not corrected values."

Good luck tomorrow! 🪄
