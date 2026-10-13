# Chapter 14. Samplers and the Horizon

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 14. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23331100. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

Can simulation, in the form of nested sampling, vertical likelihood Monte Carlo or split sampling, say anything about the Collatz conjecture in polynomial time? The chapter treats the question as rare event simulation: put an energy on the integers, declare the high energy ones interesting, and estimate the prior volume of level sets through the layer cake identity and a weighted slice sampler. The prior is Haar measure on the 2 adic integers and the two energies are the coefficient stopping time and the total stopping time. It proves that level sets of the coefficient stopping time are unions of residue classes of exactly computable volume (Chernoff rate D = 0.0347, hazard staircase, partition function finite exactly for beta up to D), that the total stopping time is almost surely infinite under Haar measure so a finite box is needed, and that a perfect split sampler for the coefficient stopping time, with no inner Markov chain, costs O(K) per draw after an O(K^2) table and reaches residues of prior mass about one in 860 million at depth 400. Beyond the horizon set by bit length, orbits of nearby integers merge, deep level sets are a few trees of merged orbits, and local samplers rarely leave the tree where they start, as a census of total stopping times below 2^26 shows. Simulation cannot decide the conjecture.

## From the literature

Nested sampling of Skilling (2006); vertical likelihood Monte Carlo and the weighted slice sampler of Polson and Scott (2015); split sampling, with its Rao Blackwellized volume curve and product estimator, of Birge, Chang and Polson (2013); the slice sampler of Neal (2003); the recursive method of uniform generation from counts of Nijenhuis and Wilf (1978), the backward step of the exact sampler; the Terras bijection and density theorem (Terras 1976, Everett 1977); the merging of consecutive integers under the map, called classical in the book, as in Garner (1985) on consecutive integers with the same total stopping time.

## New results

The perfect sampler is the authors' as a sampler (its backward step is the method of Nijenhuis and Wilf), and so are the two energies proposition, the coalescence heuristic and the tree decomposition.

<!-- results:begin -->
- [Two energies under Haar measure, the coefficient stopping time and the total stopping time](../RESULTS.md#two_energies_coefficient_and_total_stopping_time) (Proposition 14.1, page 241; proved): The level sets of the coefficient stopping time are unions of residue classes with exactly computable volume, while the total stopping time is almost surely infinite under Haar measure, which separates what a sampler can reach from what the conjecture is about. New in this book.
- [Perfect split sampling for the coefficient stopping time](../RESULTS.md#perfect_split_sampling_coefficient_stopping_time) (Theorem 14.2, page 242; proved): An exact sampler with no inner Markov chain draws residues whose coefficient stopping time exceeds a given depth, in time linear in the depth after a quadratic table, even when their prior mass is exponentially small. New in this book.
- [Coalescence of orbits shifts the total stopping time](../RESULTS.md#orbit_coalescence_shifts_total_stopping_time) (Proposition 14.4, page 245; proved): If the orbits of two positive integers meet at a point other than 1 or 2, their total stopping times differ by exactly the difference of the numbers of steps each takes to reach it. New in this book.
<!-- results:end -->

## Key terms

nested sampling, vertical likelihood, split sampling, weighted slice sampler, layer cake identity, rare event simulation, level sets, prior volume, perfect sampling, exact sampling, Gibbs sampler, coefficient stopping time, total stopping time, Haar measure, Chernoff rate, partition function, coalescence, trees of slow integers
