# Chapter 11. Cycles and Continued Fractions

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 11. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23421852. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

This chapter studies the cycle half of the Collatz conjecture (the 3x+1 problem): can the orbit of a positive integer close up into a nontrivial cycle? A cycle with r odd elements and s halvings satisfies the cycle identity, which puts s/r above log_2 3 and bounds the least element through the one sided distance from r log_2 3 up to the next integer. Continued fractions and a mediant lemma give a self contained proof of the classical fact that the first value of r not excluded by a verification range is a best one sided approximation. With Bařina's verification up to 2^68 or 2^71 this gives at least 72,057,431,991 odd elements for every nontrivial cycle, a bound the authors do not claim as new. The chapter also proves Diophantine insufficiency (infinitely many r survive every verification range), a window lemma and what arithmetic information would suffice, and covers an exact integer certificate, negative and rational cycles, Erdős problem 125, Lean 4 and a refined open problem.

## From the literature

- the cycle identity and the bound on the least element, classical since the 1970s (Böhm and Sontacchi 1978, Steiner 1978)
- the criterion of Crandall (1978) and the bound of at least 17,087,915 elements under the accelerated map of Eliahou (1993)
- the exclusion of m cycles (m ascending runs) for m = 1 (Steiner 1978), up to 68 (Simons and de Weger 2005) and up to 91 (Hercher 2023)
- the length bounds of Hercher (2018), Puchert (2018) and Hercher (2023), the last at least 137,528,045,312 odd elements, conditional on verification to 3·2^69 and so unconditional since Bařina (2025)
- the verification up to 2^68 and 2^71 of Bařina (2021, 2025)
- the optimal criterion of Halbeisen and Hungerbühler (1997) and the rational cycles of Lagarias (1990)
- the classical description of the best one sided approximations by intermediate fractions, stated and proved in Hančl and Turek (2019)
- Erdős problem 125 of Burr, Erdős, Graham and Li (1996), settled by AlphaProof Nexus (Tsoukalas and collaborators 2026)

## New results

Only statements that are the authors' own are catalogued.

<!-- results:begin -->
- [Diophantine insufficiency, infinitely many cycle lengths survive any verification range](../RESULTS.md#diophantine_insufficiency_cycle_lengths) (Proposition 11.20, page 205; proved): However far the verification of the Collatz conjecture is pushed, infinitely many cycle lengths survive the cycle inequality, so information about the approximation of log2 3 alone cannot exclude all nontrivial cycles, a point Halbeisen and Hungerbühler (1997) made earlier in another form. New in this book.
- [Window lemma, the least element of a cycle lies in an explicit window](../RESULTS.md#window_lemma_least_cycle_element) (Lemma 11.21, page 206; proved): Every element of a hypothetical Collatz cycle with m odd steps is at least an explicit lower bound and the least element is below the standard upper bound, so it lies in a window whose endpoints differ by a factor of about m over 3. New in this book.
- [What arithmetic would suffice to exclude all Collatz cycles](../RESULTS.md#what_arithmetic_would_suffice) (Corollary 11.22, page 207; conditional on the hypothesis that, for the function $f$ of the statement, every nontrivial cycle of the accelerated map with $m$ odd steps and $k$ steps in all has least element $M<m/(3\Lambda f(m))$, where $\Lambda=k\log2-m\log3$; no such bound is proved in the book): If an argument using the arithmetic of the parity word improved the standard upper bound on the least element of a cycle by a factor slightly above m over 3 for all large m, a finite verification would exclude every nontrivial cycle. New in this book.
<!-- results:end -->

## Key terms

Collatz conjecture, 3x+1 problem, nontrivial cycle, cycle length bound, cycle identity, continued fraction, convergents, semiconvergents, best one sided approximation, mediant lemma, Diophantine approximation, irrationality measure, linear forms in logarithms, Baker's method, verification range, rational cycles, negative cycles, window lemma, Erdos problem 125, AlphaProof Nexus, Lean 4, Barina, Hancl
