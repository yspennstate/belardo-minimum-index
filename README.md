# A proof of Belardo's conjecture on graphs with minimum index

`main.tex` and `main.pdf` are the paper.

`check_min_index.py` checks the quadratic-form identity and the length inequality of the paper on random inputs,
compares the index of random graphs B_L (L balanced) built from cubic multigraphs with the bound, tabulates the bound
against the bound rho(M_{n,n+s}) of Ben-Shimon, Eisner and Hoory, and, in `exhaustive` mode, reads connected graphs in
graph6 format and reports the smallest index and the structure of every graph that attains it. `verify` mode checks each
minimizer listed in `exhaustive_output.txt` against Theorem 3.

`run_exhaustive.sh` produces `exhaustive_output.txt`, the data of Table 1, and `run_half.sh` produces `half_output.txt`,
the runs added for Table 2.

`split_family.py` builds every simple graph B_L with L split (Definition 5), computes kappa(n,s) from the 2x2 matrix
Z(t), and compares both with the minimizers in the two output files (Theorem 6); its output is
`split_family_output.txt`. `check_inequalities.py` checks the inequalities used in the proof of Theorem 6 at the root
tau, in high precision, for j = 0, 1, 2 and many values of a; its output is `check_inequalities_output.txt`.

Requirements: Python 3.8 or later with numpy, networkx (for `split_family.py`) and mpmath (for
`check_inequalities.py`), and `geng` from nauty (version 2.8.9 was used) for the exhaustive runs.

    python check_min_index.py
    bash run_exhaustive.sh > exhaustive_output.txt
    bash run_half.sh > half_output.txt
    python check_min_index.py verify exhaustive_output.txt
    python split_family.py exhaustive_output.txt half_output.txt
    python check_inequalities.py
