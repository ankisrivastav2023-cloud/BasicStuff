what do the violion plots in 25 are inferreing: E:\RNA-seq-courses\scRNA-seq_python_course_2024\Demonstrations\02_iris_qc_normalisation_yuyao_mod_Feb2025_demo.ipynb
Grep "violin" (in e:\RNA-seq-courses\scRNA-seq_python_course_2024\Demonstrations\02_iris_qc_normalisation_yuyao_mod_Feb2025_demo.ipynb)
82 lines of output
Grep ""execution_count": 25," (in e:\RNA-seq-courses\scRNA-seq_python_course_2024\Demonstrations\02_iris_qc_normalisation_yuyao_mod_Feb2025_demo.ipynb)
31 lines of output
Small numbering note: cell [25] is the is_outlier function definition. The violin plots are the cell just above it, [24]:


sc.pl.violin(
  ETV6_RUNX1_1,
  ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
  multi_panel=True
)
How to read a violin plot
Each violin is a distribution turned sideways. Height = the value; width = how many cells sit at that value.


      │    ╭─╮       ← few cells up here (a thin tail)
value │  ╭─╯ ╰─╮
      │ ╭╯     ╰╮    ← WIDE = most cells are here
      │ ╰╮     ╭╯
      │  ╰─╮ ╭─╯
      │    ╰─╯       ← few cells down here
      └──────────
Think of a histogram rotated 90° and mirrored, so it looks like a violin. Scanpy also scatters the individual cells on top as dots, so you literally see each cell's position.

What each panel is telling you
Panel	Measures	Low values mean	High values mean
n_genes_by_counts	how many different genes were detected in that cell	empty droplet, dying cell, or barely-captured cell	possible doublet — two cells in one droplet
total_counts	total UMIs (RNA molecules) in that cell	poor capture, low-quality cell	large cell, or again a doublet
pct_counts_mt	% of RNA that's mitochondrial	fine — healthy	dying cell (cytoplasmic RNA leaked out, mito stayed)
What you're actually looking for
The fat part = your main population of real cells. That's your baseline "normal."

A long thin tail at the bottom of the first two panels = debris and empty-ish droplets.

A long thin tail at the top of pct_counts_mt = dying cells. This is usually the most obvious signal — a dense blob near the bottom and a sparse spray of stragglers reaching up.

Two separate bulges (a pinched, hourglass-shaped violin) = two distinct populations mixed together. In total_counts that often means real cells plus a population of empty/ambient droplets.

The point of drawing them
These plots don't decide anything for you — they're the look before you cut. You're eyeballing where the healthy population ends and the junk begins, so that when you then write nmads=5, nmads=3, or pct_counts_mt > 8, you're choosing numbers that match what this specific dataset actually looks like rather than copying someone else's.

That's why the notebook shows the violins first, then defines is_outlier in the very next cell. The sequence is: look → choose thresholds → flag → re-plot with hue= to check the flags landed where you intended.

Why all three together
A single metric can mislead. A cell with few genes and few counts and high mito is clearly dying. A cell with few genes but perfectly normal counts and mito might be a genuine, simple cell type — a red blood cell or a small resting lymphocyte — that you'd be wrong to discard.

That's exactly what the scatter plot in cell [23] adds: total_counts vs n_genes_by_counts, coloured by pct_counts_mt, putting all three on one picture so you can see which cells fail on several fronts at once.
