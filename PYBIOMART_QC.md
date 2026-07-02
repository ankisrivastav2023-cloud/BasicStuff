`pybiomart` is a small but extremely useful Python package that allows you to **query the Ensembl BioMart database directly from Python**. In scRNA-seq, it is commonly used to **convert gene IDs into gene names** and **retrieve annotations** (such as chromosome location, gene biotype, or whether a gene is mitochondrial).

Let's understand why it exists and how it fits into an scRNA-seq workflow.

---

# The problem pybiomart solves

Suppose your scRNA-seq count matrix looks like this:

| Gene ID         | Cell1 | Cell2 | Cell3 |
| --------------- | ----: | ----: | ----: |
| ENSG00000198888 |    35 |    42 |    30 |
| ENSG00000198727 |    12 |     8 |    15 |
| ENSG00000141510 |     4 |     2 |     1 |

These IDs are informative to software but not very readable to humans.

You probably want to know:

| Gene ID         | Gene Symbol |
| --------------- | ----------- |
| ENSG00000198888 | MT-ND1      |
| ENSG00000198727 | MT-CO1      |
| ENSG00000141510 | TP53        |

Rather than manually looking up thousands of genes, `pybiomart` queries Ensembl and returns these annotations automatically.

---

# What is BioMart?

Think of **BioMart** as a searchable interface built on top of the Ensembl database.

Imagine Ensembl as a giant library:

```
Ensembl Database

Genes
Transcripts
Proteins
Chromosomes
Gene Names
Gene Biotypes
Gene Locations
...
```

BioMart acts like a librarian.

Instead of searching manually, you ask:

> "Give me the gene symbol and chromosome for these Ensembl IDs."

BioMart returns a table with the requested information.

`pybiomart` lets you make these requests from Python.

---

# How does pybiomart work?

The typical workflow is:

```
Your Python script
        │
        ▼
pybiomart
        │
        ▼
Ensembl BioMart server
        │
        ▼
Gene annotations
```

You send a query, and it returns a pandas DataFrame.

---

# Example

Suppose you have

```
ENSG00000198888
ENSG00000198727
ENSG00000141510
```

You can ask for

* Gene symbol
* Chromosome
* Gene biotype

The result might look like

| Ensembl Gene ID | Gene Symbol | Chromosome | Gene Biotype   |
| --------------- | ----------- | ---------- | -------------- |
| ENSG00000198888 | MT-ND1      | MT         | protein_coding |
| ENSG00000198727 | MT-CO1      | MT         | protein_coding |
| ENSG00000141510 | TP53        | 17         | protein_coding |

Now you know which genes are mitochondrial.

---

# Why mitochondrial genes matter in scRNA-seq

One of the earliest quality control (QC) steps in scRNA-seq is calculating the **percentage of reads coming from mitochondrial genes**.

Healthy cells generally have relatively low mitochondrial RNA proportions, while stressed or dying cells often have much higher proportions because cytoplasmic RNA is lost more readily than mitochondrial RNA.

Typical workflow:

```
Raw counts
      ↓
Identify mitochondrial genes
      ↓
Calculate mitochondrial counts per cell
      ↓
Compute mitochondrial percentage
      ↓
Filter poor-quality cells
```

---

# How do we identify mitochondrial genes?

This depends on how your genes are labeled.

### Case 1: Gene symbols

Many processed datasets already use gene symbols:

```
MT-ND1
MT-ND2
MT-CO1
MT-ATP6
MT-CYB
```

These are easy to identify because they begin with `"MT-"` (in human data).

For example:

```
MT-ND1
MT-ND2
MT-CO1
MT-CO2
```

You simply select genes whose names start with `"MT-"`.

---

### Case 2: Ensembl IDs

Suppose your matrix instead contains

```
ENSG00000198888
ENSG00000198727
ENSG00000198695
...
```

Now you cannot tell which genes are mitochondrial just by looking at the IDs.

This is where `pybiomart` becomes very useful.

It can retrieve the chromosome or gene symbol corresponding to each Ensembl ID.

For example:

| Ensembl ID      | Symbol | Chromosome |
| --------------- | ------ | ---------- |
| ENSG00000198888 | MT-ND1 | MT         |
| ENSG00000198727 | MT-CO1 | MT         |
| ENSG00000141510 | TP53   | 17         |

Genes on chromosome `MT` are mitochondrial genes.

---

# Typical scRNA-seq pipeline using pybiomart

Imagine you receive data from a sequencing facility.

```
Raw Matrix

Rows

ENSG00000198888
ENSG00000141510
ENSG00000198727
...
```

Step 1

Load the matrix.

↓

Step 2

Use `pybiomart` to retrieve annotations.

↓

You obtain

| Gene ID         | Symbol | Chromosome |
| --------------- | ------ | ---------- |
| ENSG00000198888 | MT-ND1 | MT         |
| ENSG00000198727 | MT-CO1 | MT         |
| ENSG00000141510 | TP53   | 17         |

↓

Step 3

Identify mitochondrial genes.

```
Chromosome == MT
```

or

```
Gene Symbol starts with MT-
```

↓

Step 4

Calculate mitochondrial percentage for each cell.

```
mitochondrial counts
-------------------------
total counts
```

↓

Step 5

Filter poor-quality cells.

---

# Example of the calculation

Suppose one cell has

| Gene   | Counts |
| ------ | -----: |
| MT-ND1 |     30 |
| MT-CO1 |     25 |
| TP53   |      5 |
| GAPDH  |     40 |

Total counts

```
30 + 25 + 5 + 40 = 100
```

Mitochondrial counts

```
30 + 25 = 55
```

Mitochondrial percentage

```
55 / 100 × 100 = 55%
```

This is unusually high and may indicate a damaged or dying cell, depending on the tissue and experiment.

---

# Why not hard-code mitochondrial genes?

You could keep a list like

```
MT-ND1
MT-ND2
MT-ND3
...
```

However:

* species differ (human vs mouse)
* gene annotations evolve
* Ensembl releases are updated
* your dataset may use Ensembl IDs instead of symbols

Using `pybiomart` makes your pipeline more robust and reproducible because it retrieves annotations from the reference genome used by Ensembl.

---

# Where does pybiomart fit in the workflow?

```text
Raw scRNA-seq matrix
        │
        ▼
Rows = Ensembl IDs
        │
        ▼
pybiomart
        │
        ▼
Gene annotations
(Symbols, chromosome, biotype)
        │
        ▼
Identify mitochondrial genes
        │
        ▼
Calculate % mitochondrial RNA
        │
        ▼
Filter low-quality cells
        │
        ▼
Normalization
        │
        ▼
Clustering
        │
        ▼
Cell-type identification
```

## Key idea

`pybiomart` **does not analyze gene expression itself**. It is an **annotation tool**. In scRNA-seq, its role is to enrich your count matrix with biological information from Ensembl—such as gene symbols, chromosome locations, and gene biotypes—so you can perform downstream analyses like identifying mitochondrial genes, computing mitochondrial RNA percentages for quality control, or grouping genes by biological properties.
