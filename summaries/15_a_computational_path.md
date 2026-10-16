# Chapter 15. A Computational Path

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 15. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23421852. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

What can computation contribute to a proof of the Collatz conjecture, and does the majorization minimization (MM) reading say where to point it? The chapter separates verification, certificates and evidence. For the cycle half it lays out a five step programme: verification beyond 2^71.88, a Lean proof of the bound of 72,057,431,991 odd elements, a modular sieve on the numbers c_k(w) of the parity word (in runs up to m = 40 and primes below 3000 it excludes 16 shapes, all with m at most 8, of which only eight are not already excluded by size), window sharpening, and a divisibility theorem. For the divergence half it studies residue certificates and the residue sieve, whose uncovered classes grow like 2^(0.95 k), and proves that the exceptional set of 2 adic integers with infinite coefficient stopping time is closed, null and contains one point of each known negative cycle, and that every uncovered class contains infinitely many positive integers, so no finite residue computation closes the divergence half. It discusses the automated termination prover of Yolcu, Aaronson and Heule, reads Lyapunov potentials, the horizon, Switching MM and potentials with memory as constraints on any certificate, and shows why Polyak's heavy ball, gradient dominance and averaging do not transfer to the conjecture.

## From the literature

The density theorem and parity vectors of Terras (1976) and Everett (1977), with its exponential form in Lagarias (1985); the 2 adic formulation of the stopping time problem by Terras (1976) and Lagarias (1985) and the noniterative form of Bernstein (1994); the rewriting and matrix interpretation prover of Yolcu, Aaronson and Heule (2023); the verification to 2^71 of Bařina (2025) and the cycle accounting of Hercher (2023); the heavy ball of Polyak (1964), gradient dominance and the step size of Polyak (1987), and averaging of Ruppert (1988), Polyak (1990) and Polyak and Juditsky (1992).

## New results

The proposition on residue certificates is the authors' in the form stated, although each ingredient is classical; the residue sieve counts to depth 400 and the divisibility sieve experiment are also the authors'.

<!-- results:begin -->
- [Residue certificates cannot prove descent for all integers (the exceptional 2 adic set)](../RESULTS.md#residue_certificates_cannot_prove_descent_for_all_integers) (Proposition 15.1, page 267; proved): No residue certificate of finite level can certify descent for every integer, because every residue class it leaves uncovered contains infinitely many positive integers, while the 2 adic integers with infinite coefficient stopping time form a closed null set that contains a point of each known negative cycle. New in this book.
<!-- results:end -->

## Key terms

computational path, residue sieve, residue certificate, coefficient stopping time, exceptional 2 adic set, Haar measure, binary entropy, divisibility sieve, window sharpening, Lean proof of the cycle bound, termination prover, string rewriting, matrix interpretation, joint spectral radius, switched system, potential with memory, heavy ball, Polyak and Łojasiewicz inequality, Lojasiewicz, Foster drift condition, Barina
