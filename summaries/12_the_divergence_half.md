# Chapter 12. The Divergence Half

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 12. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23421852. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

This chapter treats the divergence half of the Collatz conjecture (the 3x+1 problem): can the Syracuse orbit of a positive integer grow without bound? It proves a Cesàro identity linking the logarithmic growth of an orbit to the running mean of the halving exponents (the 2 adic valuations of 3y+1), and deduces that along a divergent orbit this mean has limit superior at most log_2 3, so the lower limit of its fraction of odd steps under the accelerated map is at least 1/log_2 3, about 0.631. The Cesàro liminf conjecture is shown equivalent to the divergence half, and the obvious empirical test is shown to be circular. The model chain on residues modulo 3^a, with independent Geometric(1/2) halving exponents, has a one step law that is the exact limit for random starting values and is exactly stationary after a strides. The chapter also describes Tao's theorem, the Krasikov and Lagarias counting bound, equidistribution along single orbits, and three failed routes (equidistribution, Lyapunov certificates, uniform integrability), and it notes that if the positive density theorems of 2026 described in Chapter 13 hold, convergence to 1 on a set of logarithmic density one would already imply the conjecture.

## From the literature

- the density theorem of Terras (1976, 1979) and Everett (1977): almost every integer falls below its starting value
- the exponent log 3/log 4 of Korec (1994) for almost every orbit
- the theorem of Tao (2022): almost all orbits attain almost bounded values in logarithmic density, with Syracuse random variables, a fine scale mixing estimate modulo powers of 3, and the identities behind the model chain (loss of memory modulo 3^a, the one step recursion and the stationary law)
- the counting bound of Krasikov and Lagarias (2003): at least x^0.84 integers up to x reach 1
- the branching random walk model of preimages of Lagarias and Weiss (1992)
- the verification up to 2^71 of Bařina (2025)
- the elementary fact, noted many times, that a divergent orbit needs at least a fraction 1/log_2 3 of odd steps

## New results

Only statements that are the authors' own are catalogued.

<!-- results:begin -->
- [Cesàro identity for the halving exponents of a Syracuse orbit](../RESULTS.md#cesaro_identity) (Proposition 12.3, page 223; proved): Over N strides of the Syracuse map, the logarithmic change of the orbit equals a sum of logarithms minus N times the running mean of the halving exponents times log 2, which ties growth of the orbit to that running mean. New in this book.
- [Cesàro liminf conjecture on the running mean of the halving exponents](../RESULTS.md#cesaro_liminf_conjecture) (Conjecture 12.6, page 224; open): Along every Syracuse orbit the lower limit of the running average of the halving exponents is strictly above log2 3, a statement the book shows is equivalent to the divergence half of the Collatz conjecture. New in this book.
<!-- results:end -->

## Key terms

Collatz conjecture, 3x+1 problem, divergence half, divergent orbit, Syracuse map, halving exponent, 2 adic valuation, Cesàro identity, Cesaro identity, Cesàro liminf conjecture, running mean, Terras density theorem, Tao theorem, almost all orbits, logarithmic density, Syracuse random variable, model chain, stationary law, equidistribution, total variation, Krasikov and Lagarias, branching random walk, Lyapunov certificate, uniform integrability
