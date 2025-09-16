# 🧬Bioinformatics
This project implements algorithms for sequence alignment, alignment scores and alignment paths calculation.

> ℹ️ This project is not open source and does not grant any usage rights.
> For usage terms and legal information, see [Code Ownership & Usage Terms](#-code-ownership--usage-terms).

## 🚀This project includes:
 - 🔠Sequence processing to construct random sequences and splitting data into two datasets.
 - ⚙️Implementing Universal Sequence Alignment on the first dataset progressively, using the Needlman Wunsch algorithm for aligning pairs of sequences.
 - ⬇️Saves file or 🔄️Updates existing file with the final aligned sequences from the first dataset.
 - 🏗️Training and construction of the hmm profile, transition table and emission table using the first dataset.
 - ⬇️Saves file or 🔄️Updates existing file with the hmm profile.
 - 🧮Alignment scores and alignment paths calculation to the second dataset using the Viterbi algorithm to the trained hmm profile, transition table and emission table.
 - ⬇️Saves file or 🔄️Updates existing file with the final alignment scores and alignment paths from the second dataset.

## 🧠Technologies used:
 - Python
 - Numpy


## 🎯Purpose
This was created to explore how alignment can be done on character sequences (such as proteins) using a combination of algorithms, how Hidden Markov Model profiles, transaction and emission tables are trained and constructed through the sequences, and finally how the profiles and tables are used to calculate automatically the alignment and scores of other sequences that were not involved in training and construction. It is developed solely for academic and research purposes.


## 🧰Prerequisites

Before running the application, make sure your environment is properly configured.

- Python Version 3.9+ is recommended
Required Libraries
- Numpy (version 2.0.0 or newer)

🧪 How to Run
1. **Clone the repository (or download and decompress the ZIP file).**
```bash
git clone https://github.com/theofanistzoumakas /Bioinformatics.git
cd Bioinformatics
```

2. Confirm that you have installed the required libraries
3. Run bioinformatics.py
4. View the sequences before and after the final/completed alignment, the hmm profile (transaction table) and the alignment scores and paths from the second dataset.
5. View the txt files that includes the final aligned sequences from the first dataset, the transaction table from the first dataset and the final alignment scores and alignment paths from the second dataset.

## 🔒 Code Ownership & Usage Terms
This project was created and maintained by:

Theofanis Tzoumakas (@theofanistzoumakas)
Konstantinos Pavlis (@kpavlis)
Michael-Panagiotis Kapetanios (@KapetaniosMP)

🚫 Unauthorized use is strictly prohibited.
No part of this codebase may be copied, reproduced, modified, distributed, or used in any form without explicit written permission from the owners.
