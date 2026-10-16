# Chapter 17. Progress on the Open Problems

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 17. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23421852. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

This chapter reports progress on problems O4, S2, P3, P4, P5, P7 and P11 of chapter 19, by proof where possible, otherwise by computation. For the Collatz map it studies the cumulants of total stopping times of the accelerated map (O4, P4) and single Syracuse orbits from near 2^100000 and 2^300000 (P3). For Switching MM, the momentum method that alternates a long and a unit gradient step, it treats the exact worst case rate on smooth convex functions (P11) and the noise floor under stochastic gradients (S2, P5). It also extends the exact density of the Erdős 125 sumset A+B to 4^17 (P7). Proved: explicit worst functions bounding the convex rate below for every even N and momentum coefficient below 1; the matching upper bound for N = 2 and 4, by a performance estimation certificate, on ranges containing the Collatz momentum; and the 5/3 noise floor bound, relative to gradient descent, for arbitrary noise covariance on quadratics. Conjectured: the two branch closed form of the convex rate. Numerical evidence only: constant deficits in the cumulants of stopping times (boundedness is open), the 5/3 bound on a logistic regression problem, and no departure from independent sampling along single orbits.

## From the literature

- the worst case LR^2/(4N+2) of gradient descent on smooth convex functions, by performance estimation (Drori and Teboulle 2014), whose proof structure the chapter follows
- the interpolation conditions for smooth convex functions (Taylor, Hendrickx and Glineur 2017)
- the Huber function, the classical worst case function for gradient descent
- cumulants of first passage times from Wald's identities (Feller 1971, Gut 1988)
- the theorem of Terras (1976), which makes the Bernoulli model exact in the head of an orbit
- Erdős problem 125 (Burr, Erdős, Graham and Li 1996), solved by AlphaProof Nexus (Tsoukalas and others 2026): the lower density of A+B is zero

## New results

All five labeled statements are the authors' own: one conjecture, one theorem and three propositions, one of them a computed fit.

<!-- results:begin -->
- [Additive deficits of the second and third cumulants of total stopping times, fitted to the exact census](../RESULTS.md#stopping_time_cumulant_deficits_census_fit) (Proposition 17.1, page 290; proved): A fit to the exact census of total stopping times shows that their second and third cumulants grow at the rate of the Bernoulli model up to constant additive deficits, which accounts for the observed skewness; it is a computation on the blocks from 2^20 to 2^30, not a proof for all starting values. New in this book.
- [Conjectured exact worst case rate of Switching MM on smooth convex functions, both branches](../RESULTS.md#switching_mm_exact_convex_rate_conjecture) (Conjecture 17.2, page 292; open): A conjectured closed form for the exact worst case of Switching MM on smooth convex functions, for every even number of steps and every momentum coefficient below one, with a main branch and an overshoot branch that cross at an explicit threshold. New in this book.
- [Explicit worst case functions for Switching MM on smooth convex functions](../RESULTS.md#switching_mm_explicit_worst_functions) (Proposition 17.3, page 293; proved): Two explicit one dimensional functions, a Huber function and an overshoot function, show that the worst case of Switching MM on smooth convex functions is at least the larger of the two branches of the conjectured exact rate, for every even number of steps and every momentum coefficient below one. New in this book.
- [Exact worst case rate of Switching MM after two and four steps, by an explicit certificate](../RESULTS.md#switching_mm_exact_rate_two_and_four_steps) (Theorem 17.4, page 294; proved): For two and four steps and momentum coefficients up to about 0.64 and 0.70, an explicit certificate built from interpolation inequalities proves the conjectured exact rate of Switching MM on smooth convex functions, in ranges that contain the Collatz momentum. New in this book.
- [The five thirds noise floor bound for Switching MM under correlated gradient noise on quadratics](../RESULTS.md#switching_mm_noise_bound_correlated_noise) (Proposition 17.5, page 295; proved): On a quadratic with gradient noise of arbitrary covariance, independent from call to call, the stationary excess loss of Switching MM is at most five thirds of that of gradient descent with the same step, for every step up to one over the largest curvature and every momentum coefficient below one. New in this book.
<!-- results:end -->

## Key terms

Collatz conjecture, 3x+1 problem, total stopping time, accelerated map, Syracuse orbit, third cumulant, skewness, census, Switching MM, worst case rate, smooth convex functions, performance estimation, Huber function, noise floor, stochastic gradient, correlated noise, Erdős 125, Erdos 125, sumset, lower density, trajectory mixing
