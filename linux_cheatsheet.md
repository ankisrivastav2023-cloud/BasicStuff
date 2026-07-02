# Linux Terminal Cheat Sheet for Biologists

A beginner-friendly reference, tuned for someone who will be handling
sequencing data (FASTQ, CSV, count matrices) and running Python/R pipelines.

> **Mental model:** The terminal is just a text way of doing what you'd do
> with a mouse in Finder/File Explorer — move around folders, open/copy/delete
> files — plus the power to run analysis programs. Every command is:
> `program  [options]  [what to act on]`
> e.g. `ls -l data/` → program `ls`, option `-l`, target `data/`.

---

## 1. Where am I? — Navigating folders

| Command | What it does | Example |
|---|---|---|
| `pwd` | **P**rint **w**orking **d**irectory (where you are now) | `pwd` |
| `ls` | **L**i**s**t files in current folder | `ls` |
| `ls -l` | List with details (size, date, permissions) | `ls -l` |
| `ls -lh` | Same, sizes human-readable (KB/MB/GB) | `ls -lh` |
| `ls -a` | Show hidden files too (names starting with `.`) | `ls -a` |
| `cd folder` | **C**hange **d**irectory (go into a folder) | `cd data` |
| `cd ..` | Go **up** one folder | `cd ..` |
| `cd ~` | Go to your **home** folder | `cd ~` |
| `cd -` | Go back to the **previous** folder | `cd -` |
| `tree` | Show folder structure as a tree (may need install) | `tree data/` |

> **Tip:** Press **Tab** to auto-complete file/folder names. Type `cd Dow`
> then Tab → it fills in `Downloads/`. This prevents typos and saves time.

---

## 2. Looking at files and folders

| Command | What it does | Example |
|---|---|---|
| `cat file` | Print the whole file to screen | `cat samples.csv` |
| `less file` | Scroll through a big file (press `q` to quit) | `less bigfile.fastq` |
| `head file` | First 10 lines | `head data.csv` |
| `head -n 20 file` | First 20 lines | `head -n 20 data.csv` |
| `tail file` | Last 10 lines | `tail log.txt` |
| `tail -f file` | Watch a file update live (good for logs) | `tail -f run.log` |
| `wc -l file` | Count lines (e.g. how many reads/rows) | `wc -l samples.csv` |
| `column -t -s, file` | Pretty-print a CSV in aligned columns | `column -t -s, meta.csv` |

> **Why `head`/`less` matter:** A FASTQ or count matrix can be gigabytes.
> Never `cat` a huge file — it floods your screen. Use `head` to peek or
> `less` to scroll.

---

## 3. Creating, copying, moving, deleting

| Command | What it does | Example |
|---|---|---|
| `mkdir name` | Make a new folder | `mkdir results` |
| `mkdir -p a/b/c` | Make nested folders in one go | `mkdir -p project/data/raw` |
| `touch file` | Create an empty file (or update its timestamp) | `touch notes.txt` |
| `cp src dst` | **C**o**p**y a file | `cp data.csv backup.csv` |
| `cp -r src dst` | Copy a whole folder (`-r` = recursive) | `cp -r data/ data_backup/` |
| `mv src dst` | **M**o**v**e OR rename | `mv old.txt new.txt` |
| `rm file` | **R**e**m**ove (delete) a file — **no undo!** | `rm temp.txt` |
| `rm -r folder` | Delete a folder and everything in it | `rm -r scratch/` |

> ⚠️ **`rm` is permanent — there is no Trash/Recycle Bin.** Double-check the
> path before pressing Enter. Avoid `rm -rf` until you're confident; that
> flag force-deletes with zero warnings.

---

## 4. Finding things

| Command | What it does | Example |
|---|---|---|
| `find . -name "*.fastq"` | Find files by name pattern (here: all FASTQ) | `find . -name "*.csv"` |
| `grep "text" file` | Search **inside** a file for lines containing text | `grep "TP53" genes.txt` |
| `grep -i "text" file` | Case-insensitive search | `grep -i "chr1" data.txt` |
| `grep -c "text" file` | Count matching lines | `grep -c ">" seqs.fasta` |
| `grep -r "text" folder/` | Search inside every file in a folder | `grep -r "error" logs/` |

> **Bio example:** `grep -c ">" reads.fasta` counts sequences in a FASTA file
> (each record starts with `>`).

---

## 5. Pipes and redirection — combining commands

The `|` (pipe) sends one command's output into the next. `>` and `>>` save
output to a file. This is where the terminal becomes powerful.

| Symbol | What it does | Example |
|---|---|---|
| <code>&#124;</code> | Pipe output into next command | <code>cat data.csv &#124; head</code> |
| `>` | Save output to a file (**overwrites**) | `ls > filelist.txt` |
| `>>` | **Append** output to a file | `echo "done" >> log.txt` |
| `sort` | Sort lines | <code>cat genes.txt &#124; sort</code> |
| `uniq` | Collapse adjacent duplicate lines | <code>sort x &#124; uniq</code> |
| `uniq -c` | Count occurrences of each line | <code>sort x &#124; uniq -c</code> |

> **Worked example — count unique cell types in a metadata column:**
> ```bash
> cut -d, -f3 metadata.csv | sort | uniq -c
> ```
> reads column 3 of a CSV (`cut`), sorts it, then counts each unique value.

---

## 6. Running your analysis (Python / R / conda)

| Command | What it does | Example |
|---|---|---|
| `python script.py` | Run a Python script | `python qc.py` |
| `python` | Start interactive Python (`exit()` to quit) | `python` |
| `Rscript script.R` | Run an R script | `Rscript deseq.R` |
| `jupyter lab` | Launch JupyterLab in your browser | `jupyter lab` |
| `conda activate name` | Switch into a conda environment | `conda activate scanpy` |
| `conda deactivate` | Leave the environment | `conda deactivate` |
| `conda env list` | List your environments | `conda env list` |
| `which python` | Show which Python is currently active | `which python` |

> **Why environments matter:** scanpy and Seurat pipelines depend on specific
> package versions. A conda environment keeps each project's tools isolated so
> updating one project doesn't break another. You'll live in these daily.

---

## 7. Getting data & managing files

| Command | What it does | Example |
|---|---|---|
| `wget URL` | Download a file from a link | `wget https://.../data.tar.gz` |
| `curl -O URL` | Download a file (alternative) | `curl -O https://.../file.csv` |
| `gzip file` | Compress a file → `file.gz` | `gzip big.fastq` |
| `gunzip file.gz` | Decompress a `.gz` file | `gunzip reads.fastq.gz` |
| `zcat file.gz` | View a `.gz` file without unzipping | <code>zcat reads.fastq.gz &#124; head</code> |
| `tar -xzf file.tar.gz` | Extract a `.tar.gz` archive | `tar -xzf dataset.tar.gz` |
| `tar -czf out.tar.gz folder/` | Compress a folder into an archive | `tar -czf backup.tar.gz results/` |
| `du -sh folder/` | Show total size of a folder | `du -sh data/` |
| `df -h` | Show free disk space | `df -h` |

> **Sequencing data is almost always gzipped** (`.fastq.gz`). Tools read `.gz`
> directly, so you rarely need to unzip — use `zcat ... | head` to peek.

---

## 8. Getting help & staying safe

| Command | What it does | Example |
|---|---|---|
| `man command` | Full manual for a command (`q` to quit) | `man ls` |
| `command --help` | Quick usage summary | `grep --help` |
| `history` | Show your recent commands | `history` |
| `clear` | Clear the screen (or Ctrl+L) | `clear` |
| `Ctrl + C` | **Stop** a running command | *(key combo)* |
| `Ctrl + L` | Clear screen | *(key combo)* |
| `Ctrl + A` / `Ctrl + E` | Jump to start / end of the line | *(key combo)* |
| `↑` / `↓` arrows | Scroll through previous commands | *(keys)* |

---

## Survival rules for your first weeks

1. **Tab-complete everything.** Don't type full paths — type a few letters and hit Tab.
2. **`pwd` when lost.** If a command "can't find the file," you're probably in the wrong folder.
3. **Peek before you process.** `head`, `less`, `wc -l` before running anything heavy.
4. **`rm` has no undo.** Read the path twice. Keep backups of raw data.
5. **Ctrl+C is your escape hatch** if something hangs or runs away.
6. **Relative vs absolute paths:** `data/x.csv` is relative to where you are;
   `/home/you/data/x.csv` (starts with `/`) works from anywhere.

---

## One-page quick reference

```
pwd                 # where am I
ls -lh              # list files (readable sizes)
cd folder / cd ..   # move in / up
mkdir -p a/b        # make folders
cp -r src dst       # copy (folder)
mv src dst          # move / rename
rm file             # delete (NO undo)
head / less file    # peek at a file
wc -l file          # count lines
grep "x" file       # search inside a file
find . -name "*.gz" # find files by name
cat a | sort | uniq -c   # count unique values
conda activate env  # enter analysis environment
python script.py    # run analysis
Ctrl+C              # stop a command
man command         # get help
```

*Keep this open in a tab. Muscle memory builds fast — within a couple of weeks
these become automatic.*
