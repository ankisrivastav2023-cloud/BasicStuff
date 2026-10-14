Your guide for notebook 06 is ready: [06_DE_INTERVIEW_STUDY_GUIDE.md](Demonstrations/06_DE_INTERVIEW_STUDY_GUIDE.md). Open it and press **Ctrl+Shift+V** to read it formatted.

**How it's organised:**
- **🔥 At the top:** the 10 things to remember. If you only have 15 minutes tonight, read that box and the 90-second script at the bottom.
- **A one-picture map** of the whole notebook.
- **Steps 1–13:** each step has the code, a simple explanation, an analogy and expert questions with answers.
- **At the end:** problems in the notebook you can raise yourself, 12 rapid-fire questions, and a 90-second spoken explanation to practise.

**The one idea that matters most: pseudobulk.** Cells from the same patient are not independent, so the replicate is the *patient*, not the cell. You add up the counts per patient and cell type, then use bulk tools like edgeR. In the Harry Potter version, 1,000 students from Hogwarts don't count as 1,000 schools. If you can explain this clearly, you've covered the main question they're likely to ask.

**⚠️ One thing to check tonight: which direction the fold change points.**
- **What the code means:** it uses ETV6-RUNX1 as the baseline and compares PBMMC against it. So a positive fold change means the gene is **higher in healthy samples**.
- **What the notebook text says:** "HTR1F is upregulated in leukemic cells". That contradicts the code.
- **Why I couldn't settle it:** this computer's Python doesn't have the single-cell packages, so I couldn't run the data myself.
- **How to settle it:** Step 9 of the guide has a 2-line check you can run in the course environment.
- **In the interview,** say you always check which group is the reference. That shows real experience either way.

**Problems in the notebook worth knowing:**
- **The gene set analysis is broken.** Some printed output was accidentally pasted into a code cell, which causes a `SyntaxError`. That then causes the "repeated edges" error. Every cell after it never ran, so they have no outputs.
- **The "B (c1)" comparison mixes two things.** In the leukaemia samples it's thousands of cancer cells. In the healthy samples it's only 40–187 normal B cells. So the result is partly "cancer vs normal", not a clean comparison of the same cell type.

Good luck tomorrow! If you want a final check tonight, I can quiz you with the expert questions one at a time.
