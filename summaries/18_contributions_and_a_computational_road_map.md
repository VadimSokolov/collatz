# Chapter 18. Contributions and a Computational Road Map

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 18. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23421852. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

What does the book establish about the Collatz conjecture (the 3x+1 problem), which proof strategies does it rule out, and what computational road map remains? The chapter does not prove the conjecture. It collects the book's contributions in a table, a ledger of obstructions listing each strategy shown to fail and what survives it, and four corrections. Its organizing idea is the integrality principle: an argument that uses only a bounded number of binary digits of the starting value holds unchanged on the 2 adic integers, where the map has a fixed point and further cycles, so it cannot prove the conjecture. The road map has a shared kernel and two lanes: a cycle lane that ends at a gate needing a divisibility theorem about 2^k minus 3^m and c_k(w), and a divergence lane that searches a string rewriting formulation for a termination certificate, which would settle both halves but may not exist. One further route goes through density: given the unrefereed positive density preprints of Mazur and Lau presented in Chapter 13, a bounded form of Tao's theorem (problem B1) would also settle the whole conjecture.

## From the literature

The theorem of Tao (2022) that almost every orbit attains almost bounded values in logarithmic density; the counting exponent 0.84 of Krasikov and Lagarias (2003); verification up to 2^71 (Bařina 2025); the cycle length bound from the continued fraction of log_2 3 (Eliahou 1993, sharpened by Hercher 2023); the exclusion of cycles with m ascending runs (Simons and de Weger 2005 for small m, Hercher 2023 up to m = 91); least denominator fractions in an interval (Hercher 2023) and the best one sided approximations (Hančl and Turek 2019); the equivalence of the conjecture with termination of a string rewriting system, and the matrix interpretation search (Yolcu, Aaronson and Heule 2023); the undecidability theorem of Conway (1972); the integrality principle as folklore of the 2 adic literature (Lagarias 1985); the unrefereed positive density preprints of Mazur (2026) and Lau (2026).

## New results

The chapter collects the book's contributions and states no new theorem; its table points to the chapters that state and credit each result.

<!-- results:begin -->
No new results of this chapter are catalogued.
<!-- results:end -->

## Key terms

Collatz conjecture, 3x+1 problem, Syracuse map, computational road map, computational roadmap, ledger of obstructions, integrality principle, 2 adic integers, divergence half, cycle half, nontrivial cycle, Cesàro identity, Cesaro identity, Lyapunov function, residue certificate, optional stopping, Diophantine approximation, continued fraction of log_2 3, mediant lemma, string rewriting, matrix interpretation, termination certificate, SAT solver, positive density, bounded Tao, Lean, Switching MM, Barina, Hancl
