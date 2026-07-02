This is one of the most confusing parts of **single-cell RNA sequencing (scRNA-seq)** when you're starting out. The terms **Ensembl ID**, **Gene ID**, and **Ensembl** are related but not identical.

Let's build the idea from the ground up.

---

# What is a gene?

A gene is simply a region of DNA that can produce an RNA (and often a protein).

For example:

| Gene name | Function         |
| --------- | ---------------- |
| TP53      | Tumor suppressor |
| GAPDH     | Glycolysis       |
| ACTB      | Beta-actin       |
| CD3D      | T-cell marker    |

These names (TP53, GAPDH, etc.) are called **gene symbols**.

They are easy for humans to remember.

---

# Why isn't the gene name enough?

Gene names can change over time.

Sometimes

* a gene gets renamed
* two species have similar names
* different databases use different names

Imagine naming people.

Instead of saying

> John

you assign

> Person #00034567

The ID never changes.

Genes work the same way.

---

# What is Ensembl?

**Ensembl** is one of the world's largest genome databases.

It is developed by the European Molecular Biology Laboratory and the European Bioinformatics Institute.

It stores

* genome sequences
* genes
* transcripts
* proteins
* gene locations
* annotations

Think of Ensembl as a giant library of genomes.

---

# What is an Ensembl Gene ID?

Every gene in Ensembl gets a unique identifier.

For humans, they look like

```
ENSG00000141510
```

or

```
ENSG00000139618
```

where

* ENSG = Ensembl Gene
* the numbers uniquely identify the gene

Example:

| Gene Symbol | Ensembl ID      |
| ----------- | --------------- |
| TP53        | ENSG00000141510 |
| BRCA2       | ENSG00000139618 |
| ACTB        | ENSG00000075624 |

The Ensembl ID never changes even if the gene name changes.

---

# What about transcripts?

A single gene can produce multiple transcript isoforms.

Example:

```
Gene
 |
 +---- Transcript 1
 |
 +---- Transcript 2
 |
 +---- Transcript 3
```

Each transcript also has an Ensembl ID.

Example

```
Gene:
ENSG00000141510

Transcript:
ENST00000269305

Transcript:
ENST00000420246
```

Notice

* ENSG = gene
* ENST = transcript

---

# What is Gene ID?

The term **Gene ID** is unfortunately ambiguous.

Different databases use different gene IDs.

For example

### NCBI

```
7157
```

This is the Gene ID for TP53 in the National Center for Biotechnology Information database.

---

### Ensembl

```
ENSG00000141510
```

Same gene.

Different database.

---

### Gene Symbol

```
TP53
```

This is not technically an ID.

It is the human-readable name.

---

So the same gene can have

| Database    | Identifier      |
| ----------- | --------------- |
| Gene Symbol | TP53            |
| Ensembl     | ENSG00000141510 |
| NCBI        | 7157            |

They all refer to the same biological gene.

---

# Where do Ensembl IDs appear in scRNA-seq?

When you run an scRNA-seq experiment, the sequencing reads are aligned to a reference genome.

Suppose a read maps here:

```
Chromosome

-------------------------
         TP53
-------------------------
```

The software records

```
ENSG00000141510
```

rather than

```
TP53
```

because IDs are unambiguous.

---

# Why do many datasets show gene symbols instead?

Most analysis packages convert Ensembl IDs into gene symbols.

Instead of

```
ENSG00000141510
```

you see

```
TP53
```

because it's easier to interpret.

---

# What does a count matrix look like?

With Ensembl IDs:

| Gene            | Cell1 | Cell2 | Cell3 |
| --------------- | ----- | ----- | ----- |
| ENSG00000141510 | 10    | 3     | 0     |
| ENSG00000139618 | 5     | 1     | 2     |
| ENSG00000075624 | 100   | 95    | 120   |

After annotation:

| Gene  | Cell1 | Cell2 | Cell3 |
| ----- | ----- | ----- | ----- |
| TP53  | 10    | 3     | 0     |
| BRCA2 | 5     | 1     | 2     |
| ACTB  | 100   | 95    | 120   |

The numbers haven't changed.

Only the labels changed.

---

# Why do scRNA-seq pipelines prefer Ensembl IDs?

Because they are

* unique
* stable
* species-specific
* independent of changing gene names

Gene symbols may change over time or occasionally be reused, whereas Ensembl IDs are designed to remain stable.

---

# How are they used in Seurat?

Suppose you load data.

The row names might be

```
ENSG00000141510
ENSG00000139618
ENSG00000075624
```

Many users convert them into

```
TP53
BRCA2
ACTB
```

before visualization.

The expression values stay exactly the same.

---

# Summary

| Term                  | Example                    | Meaning                                                |
| --------------------- | -------------------------- | ------------------------------------------------------ |
| Gene Symbol           | TP53                       | Human-readable gene name                               |
| Ensembl               | Genome annotation database | Database that catalogs genes, transcripts, and genomes |
| Ensembl Gene ID       | ENSG00000141510            | Stable unique identifier for a gene in Ensembl         |
| Ensembl Transcript ID | ENST00000269305            | Stable identifier for a transcript (isoform)           |
| NCBI Gene ID          | 7157                       | Stable gene identifier in the NCBI database            |

## In scRNA-seq, think of it like this:

* **Ensembl** is the library (database).
* **Ensembl Gene ID (ENSG...)** is the unique barcode assigned to each gene in that library.
* **Gene symbol (TP53, ACTB, CD3D)** is the human-friendly label most researchers use in figures and papers.
* During analysis, the underlying expression values are the same whether rows are labeled with Ensembl IDs or gene symbols; only the identifiers differ.
