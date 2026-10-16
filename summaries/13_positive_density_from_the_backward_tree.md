# Chapter 13. Positive Density from the Backward Tree

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 13. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23421852. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

This chapter asks what recent positive density theorems say about the Collatz conjecture (the 3x+1 problem). It presents four unrefereed preprints of 2026 by Mazur, Shaik and Lau, those of Mazur and Lau with Lean formalizations: the integers whose orbit under the Collatz map reaches any target t not divisible by 3 have positive lower density, and a positive proportion of integers n reach 1 within 10.46 log n steps. The proof runs the Syracuse map backwards: it counts predecessors in the inverse tree, weighs backward words by a power of 3 over a power of 2, and uses fine scale mixing of the Syracuse random variable modulo powers of 3. The chapter explains the argument, with Lau's propagation step and the seeds it starts from, and the corollary that upper natural density one for the integers reaching 1 would force every orbit to reach 1. It shows that the forward model of the book is 2 adic and the backward model 3 adic, so that the forward walk of log x and the backward walk of the word weights are one random walk run in two directions. It adds a proposition that, granting the preprints and the verification up to 2^71, a bounded form of Tao's theorem would imply the conjecture. Its computations give the exact backward tree from 1 below 2^24, the step counts below 3 million, the first crossings of random words, and the exact offset law to level 18, where the largest Fourier coefficient stays between 3.6/d^2 and 3.9/d^2 from level 7 on and sits at a power of 2. The book compared the statements with the preprints but has not checked the proofs. It poses problems B1 (bounded Tao), B2 (the mixing rate) and B3 (the sharp rate).

## From the literature

- the bound of Krasikov and Lagarias (2003): about x^0.84 integers up to x provably reach 1
- the almost bounded values theorem (Tao 2022), with fine scale mixing and Fourier decay of the Syracuse random variable
- the positive lower density theorems of Mazur (2026), with his almost bounded theorem in natural density (2026)
- the step constant 3/log(4/3), about 10.428, of Shaik (2026)
- the explicit constants and Fourier decay of Lau (2026), whose mathematics and Lean formalization were produced by the AxiomProver system
- the propagation of positive density from a seed, Mazur (2026), made explicit by Lau (2026)
- the density one corollary of Mazur (2026, upper density) and Lau (2026, natural density)
- the verification up to 2^71 of Bařina (2025)

## New results

Only statements that are the authors' own are catalogued, not the theorems of Mazur, Shaik and Lau.

<!-- results:begin -->
- [A bounded form of Tao's theorem would imply the Collatz conjecture, given the positive density theorems](../RESULTS.md#bounded_tao_would_imply_collatz_given_positive_density) (Proposition 13.3, page 241; conditional on the positive density theorem of Mazur (2026) or of Lau (2026), recent unrefereed preprints with Lean formalizations whose proofs the book does not check: for every integer $t\ge1$ not divisible by $3$, the set of positive integers whose orbit under the Collatz map reaches $t$ has positive lower natural density; the computer verification of Bařina (2025) that every positive integer up to $2^{71}$ has an orbit that reaches $1$, which the book cites and does not repeat): If the positive density theorems of Mazur and Lau are right, the Collatz conjecture would follow from a bounded form of Tao's theorem, namely that the least value on the orbit of n is at most 2 to the power 71 for all n in a set of logarithmic density one. New in this book.
<!-- results:end -->

## Key terms

Collatz conjecture, 3x+1 problem, Syracuse problem, positive density, positive lower density, logarithmic density, upper natural density, backward tree, inverse tree, backward word, offset, Syracuse random variable, fine scale mixing, total stopping time, first crossing, propagation, seed, residue density, 2 adic, 3 adic, bounded Tao, Mazur, Shaik, Lau, AxiomProver, Lean formalization, unrefereed preprint
