# CRISPR-gRNA

## It scans DNA sequences for SpCas9 NGG PAM motifs, extracts 20nt guide RNAs on both strands, calculates GC%, and identify poly-T.
<img width="1920" height="1017" alt="Screenshot (17)" src="https://github.com/user-attachments/assets/18d28af1-c339-475b-abbc-fbee87e38c95" />

## What it does----------------

It Scans both the sense and antisense DNA strands for NGG PAM sites.

Extracts the 20-nucleotide guide sequence adjacent to each PAM.

Calculates GC content and warns if it falls outside the optimal 40-60% range.

 Identifies sequences with poly-T tracts (TTTT), which can act as unwanted termination signals.

## How to use it-------------

Run the script from your terminal or command prompt:

python "CRISPR gRNA.py"

If you run it just like that it will prompt you to paste your DNA sequence.

You can also use command-line arguments:

-s : Pass a DNA sequence string directly.

-i : Provide a text or FASTA file containing your sequence.

## Output------------------

The script prints a clean table in the terminal showing:

Position and strand

20nt guide sequence and PAM

GC percentage

Evaluation (Optimal, Suboptimal GC, or Poly-T risk)

## Where it is used-------

This helps researchers design guide RNAs for CRISPR gene editing experiments. Instead of checking DNA sequences manually, this script automates the process of finding valid Cas9 binding sites and filtering out bad guides (like those with poor GC content or early termination signals) before testing them in a wet lab.

