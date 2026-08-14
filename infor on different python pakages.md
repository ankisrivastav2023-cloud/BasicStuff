cellbender vs cell ranger and which one is used at which step in the sc rna analysis
They're not alternatives — they sit at different points in the pipeline. CellRanger comes first and produces counts; CellBender comes after and cleans those counts.

The pipeline position

FASTQ files (raw sequencer reads)
    │
    ▼
┌─────────────────────────────────────┐
│ CELL RANGER  (10x Genomics)         │  ← primary processing
│ • align reads to genome (STAR)      │
│ • correct barcodes + UMIs           │
│ • count reads per gene per barcode  │
│ • call which droplets = real cells  │
└─────────────────────────────────────┘
    │
    ├── raw_feature_bc_matrix.h5     (ALL droplets, ~500k–1M barcodes)
    └── filtered_feature_bc_matrix.h5 (its guess at real cells, ~5k)
    │
    ▼  (feed the RAW one here)
┌─────────────────────────────────────┐
│ CELLBENDER remove-background        │  ← optional denoising
│ • learns the "ambient RNA soup"     │
│ • subtracts it from each cell       │
│ • re-calls cells (empty vs real)    │
└─────────────────────────────────────┘
    │
    ▼ cellbender_output_filtered.h5
┌─────────────────────────────────────┐
│ SCANPY  (your Python course)        │
│ QC → normalise → HVG → PCA → UMAP   │
│ → clustering → marker genes         │
└─────────────────────────────────────┘
What each actually does
Cell Ranger	CellBender
Made by	10x Genomics	Broad Institute (open source)
Input	FASTQ (raw reads)	Cell Ranger's raw matrix
Output	Count matrix (genes × barcodes)	Cleaner count matrix
Job	Turn sequences into counts	Remove ambient contamination
Optional?	No — you need counts	Yes — many analyses skip it
Hardware	CPU, lots of RAM, hours	Ideally GPU, ~30–60 min
The problem CellBender solves
When cells are dissolved during dissociation, some break open and spill mRNA into the buffer. Every droplet — including empty ones — picks up a bit of this free-floating "soup". So a T cell can show a few counts of haemoglobin genes that it never expressed; they just came along for the ride.

Symptoms you'd see in Scanpy: marker genes appearing in clusters where they make no biological sense, muddy cluster boundaries, inflated doublet estimates.

CellBender models each droplet as real cell expression + ambient soup + noise and estimates the soup profile from the tens of thousands of empty droplets (which contain only soup). It then subtracts that estimate.

Two practical points that trip people up
1. Give CellBender the raw matrix, not the filtered one. It needs the empty droplets to learn the soup profile. Passing filtered_feature_bc_matrix.h5 throws away exactly the data it needs.

2. CellBender replaces Cell Ranger's cell calling. Its output is its own filtered set of cells, so you don't apply Cell Ranger's filtered barcodes on top.

For your course
Most teaching courses (including this one) start from a Cell Ranger matrix and go straight into Scanpy — CellBender is an extra step you'd add for real project data, especially single-nucleus RNA-seq, where ambient RNA is worst.

Lighter alternatives if you ever want ambient correction without the GPU: SoupX (R) and decontX (R, via celda). CellBender generally performs best but is the heaviest to run.

Your repo has a Data/CellRanger_Outputs/ folder — if you want, I can look at what's in it and show which file the course notebooks are reading.

is seaborn a python package?
Yes — seaborn is a Python package for making statistical plots.

What it is: a plotting library built on top of matplotlib. Matplotlib is powerful but verbose; seaborn wraps it so common statistical charts take one line instead of fifteen, with nicer defaults.

Installing and importing:


pip install seaborn        # or: conda install seaborn

import seaborn as sns      # sns is the near-universal alias
Typical use:


import seaborn as sns
import matplotlib.pyplot as plt

sns.violinplot(data=df, x="cluster", y="n_genes")
plt.show()
Where it shows up in scRNA-seq: Scanpy uses matplotlib internally for its sc.pl.* plotting functions, but seaborn is handy when you want a custom QC plot that Scanpy doesn't provide — e.g. a scatter of counts vs. genes coloured by mitochondrial percentage, or a violin plot comparing a metric across samples. Since adata.obs is a pandas DataFrame, it feeds straight into seaborn:


sns.scatterplot(data=adata.obs, x="total_counts", y="n_genes_by_co
