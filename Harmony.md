# Harmony, step by step (the simple version)

## Why we need it: the problem in one picture

Say you have **two batches** (1a and 1b), and each contains **T cells** and **B cells**. Batch 2 was sequenced a bit differently, so **all of its cells are pushed slightly to the right** in PCA space. That shift is the batch effect.

```
BEFORE HARMONY (PCA space)

   T cells:   ●●●  ○○○        ● = batch 1
                                ○ = batch 2
   B cells:   ●●●  ○○○

→ The cells split by BATCH (left/right), which is wrong.
→ They should split only by CELL TYPE (top/bottom).
```

The goal is to slide the ○ cells left so they sit on top of the ● cells, **without** mixing T cells with B cells.

🧙 **The analogy used below:** Hogwarts students (batch 1) and Durmstrang students (batch 2) are in the Great Hall. Durmstrang students all wear heavy fur cloaks, so they look different. We want to group everyone by **what they're good at** (Seekers, Keepers) and **ignore the cloak**.

---

## Step 0: The input

Harmony gets two things:
1. **The PCA coordinates of every cell** (`X_pca`). It never looks at gene counts.
2. **The batch label of every cell** (`batch_1` / `batch_2`).

> 🧙 Harmony doesn't read every student's full school record (the genes). It only looks at a short summary card (the PCA) plus which school they came from (the batch).

---

## Step 1: Put cells into rough groups (soft clustering)

Harmony creates a set of **clusters** in PCA space and assigns cells to them.

The important word is **soft**. A cell isn't forced into exactly one group. It gets a **probability** for each cluster:

> Cell #42: 80% cluster A, 15% cluster B, 5% cluster C

> 🧙 Students aren't forced into one Quidditch position. Harry is "mostly Seeker, slightly Chaser." This matters for cells sitting between two cell types, because they won't get a harsh, all-or-nothing correction.

**Term:** *soft k-means clustering.*

---

## Step 2: Add the rule "every group must be mixed" (the diversity penalty, θ)

With normal clustering, the batch effect would win. The ● T cells would form one cluster and the ○ T cells another, because they sit in different places.

Harmony adds a **penalty**: if a cluster is made mostly of one batch, cells are **pushed toward clusters that are better mixed**.

- It compares the batch mix it **expects** (e.g. 50% batch 1 and 50% batch 2, if the batches are equal size) with the mix it **observes** in each cluster.
- A cluster that is 95% batch 2 gets penalised, so assignments shift until each cluster contains both batches.

> 🧙 Dumbledore tells the Sorting Hat: "Group students by skill, but I don't want a table that's all Durmstrang." The Hat is then encouraged to put Durmstrang Seekers at the same table as Hogwarts Seekers.
>
> 🌙 Rhysand builds his Inner Circle by nature and loyalty, not by birthplace. θ is how strongly he insists on that.

**θ (theta)** is the strength of this rule:
- **θ = 0:** no mixing rule, so ordinary clustering.
- **Higher θ:** more forced mixing. **Too high** risks **over-correction**, where different cell types get merged just to make clusters mixed.
- The default is θ = 2.

Now each cluster holds the **same cell type from both batches**, e.g. "the T cell cluster" with ● and ○ together.

---

## Step 3: Measure the batch shift inside each cluster

Look inside **one cluster**, say the T cell cluster:
- Where do the **batch 1** T cells sit, on average?
- Where do the **batch 2** T cells sit, on average?
- The **difference** between them is the **batch effect for this cluster**.

```
Inside the T-cell cluster:
  batch 1 centre: position 2
  batch 2 centre: position 5
  → batch 2 is shifted +3 → that's the batch effect here
```

Because they're the **same cell type**, this difference can only be technical. That's the logic that makes Harmony work.

This is done **separately for every cluster**, because the batch effect can be **different for each cell type** (T cells might shift +3, B cells only +1).

> 🧙 Among the Seekers, Durmstrang students look "3 cloak-units bulkier." Among the Keepers, only "1 cloak-unit." Each group has its own cloak size.

**Term:** a **linear model** (technically a *ridge regression*) fitted **per cluster**. Together these form a **mixture of experts**: each cluster is an "expert" for its own region of the data.

---

## Step 4: Correct each cell

Now each cell is **moved** to cancel out its batch shift.

Because each cell belongs **softly** to several clusters, its correction is a **weighted blend**:

> Cell #42 is 80% cluster A and 20% cluster B, so
> correction = 0.8 × (cluster A's shift) + 0.2 × (cluster B's shift)

Only the **batch part** is removed. The cell keeps its biological position.

```
AFTER CORRECTION

   T cells:   ●○●○●○        ← batches overlap
   B cells:   ●○●○●○        ← but T and B stay apart ✅
```

> 🧙 Every Durmstrang student's cloak is removed. A student who's 80% Seeker and 20% Chaser gets mostly the "Seeker-sized" cloak taken off.
>
> 🌙 The glamour lifts, and what's left is each person's true form.

---

## Step 5: Repeat until it stops changing

After correcting, the cells are better aligned, so **the clustering in Step 1 gets better**. Better clusters give a **more accurate estimate of the shift** in Step 3, which gives a **better correction** in Step 4.

So Harmony **loops: cluster → measure → correct → cluster again → ...** until the cells stop moving much. That's called **convergence**.

> 🧙 It's like Neville's potion-making improving over the years: each attempt fixes the last one's mistakes, until the result is good enough.

---

## Step 6: The output

Harmony gives back **new, corrected PCA coordinates**: `adata.obsm['X_pca_harmony']`.

What it does **not** change:
- ❌ Raw counts
- ❌ Normalised/log counts (`.X`)
- ❌ The original `X_pca`, which is still there

That's why you then run:
```python
sc.pp.neighbors(adata, use_rep='X_pca_harmony')   # rebuild the neighbour graph on the corrected coordinates
sc.tl.umap(adata)                                   # new UMAP, where the batches should now overlap
```

---

## The whole thing in one flow

```
PCA + batch labels
       ↓
① Soft-cluster the cells
       ↓
② Force clusters to contain every batch (θ penalty)
       ↓
③ In each cluster, measure: "how far is batch 2 from batch 1?"
       ↓
④ Move each cell back by that amount (weighted by cluster membership)
       ↓
⑤ Repeat ①–④ until the cells stop moving
       ↓
X_pca_harmony → neighbours → UMAP
```

---

## Parameters to name-drop

| Parameter | Meaning in plain words |
|---|---|
| **θ (theta)** | How strongly clusters are forced to be batch-mixed. Higher means more mixing, with a risk of over-correction. |
| **σ (sigma)** | How "soft" the clusters are. Higher means cells spread their membership over more clusters. |
| **λ (lambda)** | How cautious the correction is. Higher means smaller, more conservative shifts. |
| **nclust** | Number of clusters Harmony uses internally. These are **not** your final cell-type clusters. |
| **max_iter** | The maximum number of loops. |

---

## Say it in the interview (20 seconds)

> "Harmony takes the PCA embedding and batch labels. It softly clusters the cells, with a diversity penalty, theta, that forces each cluster to contain cells from all batches. Within each cluster, the same cell type is present in every batch, so any positional difference between batches must be technical. Harmony estimates that shift with a per-cluster linear model and subtracts it from each cell, weighted by the cell's cluster membership. It iterates until convergence and outputs a corrected embedding, `X_pca_harmony`. The expression matrix itself is never modified."

## Two follow-ups they might ask

- **"Why not just subtract one global shift for batch 2?"** Because the batch effect differs between cell types. Correcting per cluster handles that, making the overall correction effectively **non-linear** even though each local correction is linear.
- **"What's the risk?"** **Over-correction**, especially when θ is too high or when one batch contains a cell type the other batch doesn't have. Harmony may pull that unique population onto the wrong cell type. Always check marker genes after integration.
