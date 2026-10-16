# Chapter 10. Three Point Inequalities and Similar Triangles

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 10. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23421852. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

This chapter asks what the three point inequality for Bregman divergences, the core of convergence proofs for majorization minimization (MM) and accelerated gradient methods, says about Switching MM (the authors' method) and the Collatz map (the 3x+1 problem). It uses relative smoothness, the method of similar triangles and performance estimation. For convex L smooth functions it proves that Switching MM with momentum beta, whose long step 1+beta is at most the golden ratio, has function value error at most L R^2/((2+beta) N) at even N (R the initial distance to a minimizer), and that the 1/J^2 rate in J similar triangles steps survives any switching with plain gradient steps. It conjectures the exact worst case L R^2/((4+2 beta) N+2) at even N for beta up to 0.68614. For the Syracuse map it proves an exact Bregman decomposition of log y_N minus log y_0 into a martingale term (in the stochastic model) and an Itakura Saito dissipation, bounded below by 0.029092 at every odd y at least 3, and explains why the optimality of the step is missing.

## From the literature

- the three point identity for Bregman divergences (Chen and Teboulle 1993) and the classical three point inequality
- relative smoothness (Bauschke, Bolte and Teboulle 2017; Lu, Freund and Nesterov 2018) and the classical 1/k rate of the MM step
- Nesterov's accelerated method (Nesterov 1983, 2004), the accelerated proximal method (Beck and Teboulle 2009), and the method of similar triangles with its classical energy proof (Gasnikov and Nesterov 2018)
- the exact worst case L R^2/(4N+2) of gradient descent (Drori and Teboulle 2014) and the interpolation conditions behind performance estimation (Taylor, Hendrickx and Glineur 2017)
- long step schedules that beat gradient descent (Grimmer 2024; Altschuler and Parrilo 2025)
- the Itakura Saito divergence as a Bregman divergence (Banerjee, Merugu, Dhillon and Ghosh 2005)

## New results

Only statements that are the authors' own are catalogued.

<!-- results:begin -->
- [Convex rate for Switching MM when the long step is at most the golden ratio](../RESULTS.md#switching_mm_convex_rate_by_three_point_inequality) (Theorem 10.5, page 173; proved): On smooth convex functions the bound for Switching MM after N steps is L R^2/((2+beta) N), against L R^2/(2N) for gradient descent by the same argument, provided the long step 1+beta is at most the golden ratio, as it is for the book's beta = log_2 3 minus 1. New in this book.
- [Exact convex worst case of Switching MM (conjecture)](../RESULTS.md#switching_mm_exact_worst_case_conjecture) (Conjecture 10.6, page 174; open): Conjectured exact worst case L R^2/((4+2 beta) N+2) of Switching MM on smooth convex functions for even N and beta up to about 0.686, supported by performance estimation computations and proved in a later chapter for N = 4, and for N = 2 when beta is at most 0.6446. New in this book.
- [Acceleration by similar triangles survives every switching rule](../RESULTS.md#acceleration_survives_every_switching_rule) (Proposition 10.7, page 176; proved): Mixing similar triangles steps with plain gradient steps in any order, even one that depends on the iterates, never increases the energy, so the 1/J^2 rate in the number J of similar triangles steps is kept. New in this book.
- [Bregman decomposition of the logarithmic change along a Syracuse orbit](../RESULTS.md#bregman_decomposition_of_a_syracuse_orbit) (Proposition 10.8, page 177; proved): Along any Syracuse orbit the change of the logarithm is exactly a linear term minus a sum of Itakura Saito divergences, and in the stochastic model the linear term is a martingale and the dissipation has mean log(4/3) per step, the size of the drift. New in this book.
- [Uniform lower bound on the Bregman dissipation of the Syracuse map](../RESULTS.md#uniform_dissipation_of_the_syracuse_map) (Proposition 10.9, page 179; proved): At every odd y at least 3 the Itakura Saito divergence of the Syracuse multiplier S(y)/y from 1 is at least 0.029092, with equality only at y = 9, and its limit inferior as y grows is 0.037682. New in this book.
<!-- results:end -->

## Key terms

three point inequality, three point identity, Bregman divergence, Itakura Saito divergence, majorization minimization, MM algorithm, relative smoothness, Nesterov acceleration, accelerated gradient method, method of similar triangles, performance estimation, semidefinite programming, long step gradient descent, golden ratio, Switching MM, Syracuse map, Collatz conjecture, 3x+1 problem
