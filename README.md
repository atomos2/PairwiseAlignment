# 🧬 Pairwise Sequence Alignment using Biopython

This project performs **pairwise sequence alignment** using Biopython's `PairwiseAligner`. The user can align two sequences in either **global** or **local** mode and receive statistics such as alignment score, identity percentage, gap percentage, and aligned length.

---

## 📌 Features

- Supports both **Global** and **Local** alignment
- Computes:
  - 🔹 Alignment score
  - 🔹 Identity percentage
  - 🔹 Gap percentage
  - 🔹 Aligned length
- Handles user input paths for FASTA files


---

## 🗂️ Project Structure
pairwise-sequence-alignment/
├── data/
│ ├── example_query.fasta
│ └── example_target.fasta
├── pairwise_alignment.py
├── pairwise_alignment.ipnyb
├── README.md
├── requirements.txt

## ▶️ How to Use

1. Clone the repository or download the script.
2. Make sure you have [Python 3](https://www.python.org/) installed.
3. Install the required package:
   ```bash
   pip install biopython
4. Run the script:
   ```bash
   python pairwise_alignment.py


5. Sample output: 
Enter the FASTA path of your query: data/example_query.fasta
Enter the FASTA path of your target: data/example_target.fasta
Enter the mode of pairwise alignment, 'global' or 'local': local

The top alignment for the Local alignment is:
AC-GTGGAT
|| ||||||
ACGTCGGAT

The score for the top the alignment is: 14.5
The percentage of identity is 87.5%
The percentage of gap is 12.5%
The aligned length is 8

## 👨‍💻 Author
Iftikhar Alam
LinkedIn: https://www.linkedin.com/in/iftikhar-alam-07b05b287
Email: alamiftikhar0006@gmail.com

## 📚 Reference:
1. https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/download/lectures/PCB_Lect02_Pairwise_allign.pdf
2. https://biopython.org






