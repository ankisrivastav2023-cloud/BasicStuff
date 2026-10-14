# Notebook 06 in super simple words

## The question
**Same type of cell, two kinds of people (leukaemia patients vs healthy donors). Which genes behave differently?**

## The 6 steps

**1. Add up the cells per person (pseudobulk)**
Each person has thousands of B cells. Add up all their B-cell counts into **one total per person**.
👉 *Why?* Cells from one person are like students from one school. 1,000 students from Hogwarts still only tell you about **one school**. The **person** is what we count, not the cell.

**2. Throw out tiny groups**
If a person has fewer than **20** cells of a type, skip that group. Too few cells give an unreliable result.

**3. Throw out quiet genes**
Remove genes that are barely detected (31,801 → 14,543 genes). You can't tell whether a whisper got louder.

**4. Run edgeR (the statistics)**
edgeR compares healthy vs leukaemia for every gene. It's built for **raw counts**. With only 3–4 people per group, it **borrows information from all genes** to make better guesses.

**5. Read the results**
- **log fold change** = how big the difference is (+ means higher in healthy, − means higher in leukaemia, in this notebook)
- **adjusted p-value (FDR)** = how much to trust it. **Below 0.05 = significant**
- **Volcano plot**: genes far left/right and high up are the most important

**6. Find the bigger story (gene set analysis)**
Check whether the changed genes belong to known pathways (e.g. "immune response"). Like noticing many suspects at a meeting are Death Eaters, which is more than chance would explain.

Then **repeat for every cell type** and count how many genes changed in each.

## The 3 sentences to say in the interview
1. *"I use **pseudobulk**: the patient is the replicate, not the cell, otherwise you get false positives."*
2. *"I use **raw counts with edgeR** and **FDR-adjusted p-values**, because we test thousands of genes."*
3. *"I always check **which group is the baseline**, so I know what 'up' really means."*

## One-line memory hook
**Sum per person → drop small groups → drop quiet genes → edgeR → FDR < 0.05 → pathways.**

I can also add this short version to the top of [06_DE_INTERVIEW_STUDY_GUIDE.md](Demonstrations/06_DE_INTERVIEW_STUDY_GUIDE.md), so everything is in one place.
