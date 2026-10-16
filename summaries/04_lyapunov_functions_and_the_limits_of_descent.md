# Chapter 4. Lyapunov Functions and the Limits of Descent

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 4. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23421852. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

This chapter asks whether the heuristic that log x decreases on average along a 3x+1 (Collatz) orbit can be turned into a Lyapunov function, a potential certifying that every orbit reaches 1. It defines Lyapunov certificates for discrete maps, gives a discrete criterion and Foster's drift criterion for Markov chains, and shows that a certificate exists if and only if the Collatz conjecture holds. The positive results are the stride descent theorem for the Syracuse map, an exact two stride analysis, and bounds on the expected passage time of the stride model with its overshoot, from Wald's identity. The negative results are three obstructions. No strictly increasing function of x is a strict or strong Lyapunov function for one stride, because every x congruent to 3 mod 4 ascends. No monotone certificate over any fixed number of strides exists, because runs of ascending strides can be arbitrarily long. In majorization minimization (MM), the tangent line majorizer of log x has a surrogate that drifts upward while log x drifts downward in the model, and power majorizers with exponent below the Wald root 1 repair the drift in the model but not the pointwise failure.

## From the literature

Lyapunov (1892, English translation 1992) for Lyapunov functions and LaSalle (1986) for discrete systems; Foster (1953) and Meyn and Tweedie (1993) on drift criteria; the heuristic for log x of Terras (1976), Crandall (1978) and Lagarias and Weiss (1992); the density theorem of Terras (1976) and Everett (1977); the equivalence of the conjecture with finiteness of every stopping time, folklore implicit in Terras (1976), see Lagarias (1985, 2010); the classical identity for runs of ascents (Terras 1976; Hercher 2023); the MM principle of Lange, Hunter and Yang (2000) and Hunter and Lange (2004).

## New results

The new results are the stride descent theorem, the two stride analysis, the passage time bounds, the power majorizers and three obstructions to certificates.

<!-- results:begin -->
- [Passage time of the stride model with overshoot](../RESULTS.md#stride_model_passage_time) (Proposition 4.5, page 70; proved): In the random stride model of the Collatz map, the expected number of strides for log x to fall to zero or below is the starting log x plus the mean overshoot, which lies between log 2 and 2 log 2, divided by log(4/3), the size of the drift. New in this book.
- [Stride descent theorem for the Syracuse map](../RESULTS.md#stride_descent) (Theorem 4.7, page 71; proved): One Syracuse stride lowers log x exactly when the halving exponent is at least 2 and raises it exactly when x is 3 mod 4, with sharp bounds on the size of the change. New in this book.
- [No monotone strong Lyapunov function for a single Syracuse stride](../RESULTS.md#no_monotone_single_stride_lyapunov) (Theorem 4.8, page 73; proved): No strictly increasing function of x is a strict or strong Lyapunov function for the Syracuse map, because every integer congruent to 3 mod 4, half of all odd integers, ascends in one stride. New in this book.
- [Exact behaviour of two Syracuse strides](../RESULTS.md#exact_two_stride_behaviour) (Proposition 4.12, page 74; proved): Over two Syracuse strides an odd integer decreases exactly when the two halving exponents sum to at least 4, which holds for half of the odd integers, and the largest two stride ascent is log(17/7), attained only at 7. New in this book.
- [No fixed horizon monotone certificate for Collatz descent](../RESULTS.md#no_fixed_horizon_monotone_certificate) (Corollary 4.15, page 78; proved): For every fixed number J of strides, no strictly increasing function of x certifies descent over J Syracuse strides, because the odd integers congruent to minus 1 modulo 2 to the J+1 ascend for J strides in a row. New in this book.
- [The linear majorizer of log cannot certify Collatz descent](../RESULTS.md#linear_majorizer_cannot_certify_descent) (Proposition 4.19, page 79; proved): In the stochastic model the tangent line majorizer of log x has a positive mean increment per step of the standard map while log x has a negative one, so the linear surrogate cannot certify descent. New in this book.
- [Power majorizers of log and the Wald root of the stride model](../RESULTS.md#power_majorizers_of_log) (Proposition 4.20, page 81; proved): In the stride model the power majorizers of log with exponent strictly between 0 and 1 have strictly negative mean increment, interpolating between the true drift log(3/4) and the zero drift of the linear majorizer at the Wald root 1. New in this book.
<!-- results:end -->

## Key terms

Collatz conjecture, 3x+1 problem, Syracuse map, Lyapunov function, Lyapunov certificate, potential function, Foster's criterion, drift condition, passage time, overshoot, Wald identity, stride model, stopping time, monotone certificate, chained ascents, log squared potential, majorization minimization, MM algorithm, tangent line majorizer, power majorizer, Wald root, martingale duality.
