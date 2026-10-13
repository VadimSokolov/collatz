# Chapter 12. The Path to Collatz

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 12. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23331100. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

What would a proof of the Collatz conjecture (the 3x+1 problem) have to contain, and where does each part of the book sit relative to it? The chapter splits the conjecture for the Syracuse map into a divergence half, whether an orbit can be unbounded, equivalent to the Cesàro conjecture on the running mean of the halving exponents relative to log_2 3, and a cycle half, whether a nontrivial cycle exists, governed by the cycle identity and one sided Diophantine approximation of log_2 3. For the divergence half it examines equidistribution modulo powers of 2, Lyapunov certificates of descent and uniform integrability of the iterate martingale, and shows where each fails: statistics over random starting values do not transfer to one orbit, and monotone certificates and optional stopping are blocked. For the cycle half the verified range (2^68, now 2^71) and the semiconvergent ladder of log_2 3 give at least 72,057,431,991 odd elements in a nontrivial cycle, no irrationality measure can close the gap, and the window lemma shows that what is missing is arithmetic information about the parity word, such as whether 2^k minus 3^m can divide c_k. It lists the constants a mechanism must reproduce and the role of Lean verification.

## From the literature

The split into a divergence half and a cycle half (Lagarias 1985, 2010); the density theorem and uniform distribution of parity vectors of Terras (1976) and Everett (1977); the parity heuristic of Crandall (1978) and the stochastic models of Lagarias and Weiss (1992); the theorem of Tao (2022) on almost bounded values of almost all orbits; the reduction of cycles to Diophantine approximation of log_2 3 by Crandall (1978), Steiner (1978) and Eliahou (1993); best one sided approximations of Hančl and Turek (2019); rational cycles in Halbeisen and Hungerbühler (1997); verification to 2^71 by Bařina (2025); the cycle length record of Hercher (2023).

## New results

The chapter states no theorem of its own; its Notes credit the proof that Diophantine information cannot exclude all cycles, and the window corollary, to the authors (Chapter 10).

<!-- results:begin -->
No new results of this chapter are catalogued.
<!-- results:end -->

## Key terms

Collatz conjecture, 3x+1 problem, Syracuse map, accelerated map, divergence half, cycle half, nontrivial cycle, cycle identity, Cesàro conjecture, Cesaro conjecture, unbounded orbit, Lyapunov certificate, uniform integrability, optional stopping, mediant lemma, semiconvergent ladder, one sided Diophantine approximation, irrationality measure of log_2 3, window lemma, verification range, Lean formalization, Barina, Hancl, Hungerbuhler
