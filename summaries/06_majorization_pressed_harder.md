# Chapter 6. Majorization Pressed Harder

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 6. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23421852. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

This chapter presses majorization minimization (MM) as far as it goes on the Collatz map. The multiplicative majorizer of a Syracuse stride, summed around a cycle, gives the cycle inequality, and its minorizer is the window lemma, so the cycle half becomes a question of how tight the majorizer can be made along a closed orbit. For the divergence half the chapter proves impossibility results for pointwise potentials. A potential log x plus a function of the residue modulo 2 to the d cannot decrease at every step for all large x, for any d, by a negative cycle criterion from shortest path theory: the residue graph has cycles of positive weight, the loop at minus 1 and the images of the cycles of the map on negative integers, and its maximum cycle mean is log of 3 over 2, attained by the loop. Any potential log x plus a correction that decreases at every stride needs a surplus growing linearly in the number of trailing ones of x. Block majorization over a run of ascents gives an exact block map and a renewal process under the model, but run length potentials fail at a fraction of odd integers that no weight brings below 0.143 in the model. In expectation, majorization is the supermartingale method.

## From the literature

The negative cycle criterion is the Bellman and Ford condition of shortest path theory (Ford 1956, Bellman 1958), and the maximum cycle mean is computed by the algorithm of Karp (1978). The parity vector theorem of Terras (1976) and Everett (1977). The cycle identity and the cycles of the map on negative integers are classical, see Lagarias (1985). The classical identity for runs of ascents, see Hercher (2023), who cites Terras (1976). Steiner (1978) showed that the trivial cycle is the only cycle consisting of a single block, a circuit, see also Lagarias (1985).

## New results

The new results are the residue potential theorem, the surplus along chains and the run length potential; the cycle identity and the block map are classical.

<!-- results:begin -->
- [No residue potential for the Collatz map at any modulus](../RESULTS.md#no_residue_potential_for_any_modulus) (Theorem 6.3, page 101; proved): No potential made of log x plus a function of residues modulo a power of 2 decreases at every step of the accelerated Collatz map for all large x, because the residue graph has cycles of positive weight, such as the loop at minus 1. New in this book.
- [Surplus of a pointwise potential along runs of trailing ones](../RESULTS.md#surplus_along_chains) (Proposition 6.4, page 102; proved): A potential log x plus a correction that decreases at every Syracuse stride needs an unbounded correction, with a surplus growing linearly in the number of trailing ones of x. New in this book.
- [Run length potential and its failure fraction](../RESULTS.md#run_length_potential) (Proposition 6.6, page 104; proved): With weight c above log(3/2), the run length potential log x plus c times the 2 adic valuation of x+1 decreases at all but finitely many ascending strides, but under the stride model it fails to decrease at a positive fraction of odd integers, 0.2087 at weight c equal to 1. New in this book.
<!-- results:end -->

## Key terms

majorization minimization, MM algorithm, multiplicative majorizer, cycle inequality, cycle identity, window lemma, residue potential, residue graph, negative cycle, shortest path, maximum cycle mean, Karp's algorithm, negative integer cycles, surplus, trailing ones, block map, circuit, Steiner's theorem, run length potential, renewal process, stopping time.
