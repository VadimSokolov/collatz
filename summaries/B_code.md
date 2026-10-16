# Appendix B. Code

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), appendix B. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23421852. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

Where the code of the book is. Every number computed in the book was produced by a short Python script, and numbers quoted from the literature are cited where they appear. The scripts use mpmath, NumPy, SciPy, SymPy and Matplotlib, with Numba for the census of total stopping times and for the gallery, CVXPY with the Clarabel solver for performance estimation, and scikit-learn only to load one data set; they run in seconds to minutes on a laptop and are in the folder book/scripts of the public repository, released under the MIT License. The appendix lists two scripts in full. The first computes the partial quotients of log_2 3, builds the candidate rungs of the semiconvergent ladder, keeps the records of the one sided distance and returns the first rung not excluded by a verified range, which by the one sided ladder theorem of Chapter 11 and the mediant lemma is the exact cycle bound. The second derives symbolically the slow mode noise amplification factor of Switching MM. An index lists all other scripts by chapter, including the script that writes the table of total stopping times below 2^26 that most scripts of Chapter 16 read.

## From the literature

The script reproduces the cycle bound of Eliahou (1993) at verified range 2^40; its bound at 2^68 is unchanged at the range 2^71 of Bařina (2025).

## New results

This appendix lists code and proves nothing.

<!-- results:begin -->
No new results of this appendix are catalogued.
<!-- results:end -->

## Key terms

code, Python scripts, reproducibility, exact cycle bound script, ladder rungs, mpmath, SymPy, noise amplification derivation, script index, MIT License, GitHub repository, Barina
