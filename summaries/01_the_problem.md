# Chapter 1. The $3x+1$ Problem

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 1. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23331100. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

This chapter states the 3x+1 problem, also called the Collatz conjecture, the Syracuse problem, Kakutani's problem, Ulam's problem and Hasse's algorithm: halve an even integer, replace an odd integer n by 3n+1, and ask whether every positive integer reaches 1. It separates the two ways the conjecture can fail, a divergent orbit and a nontrivial cycle (the divergence half and the cycle half), introduces the standard, accelerated and Syracuse maps and stopping times, and works examples such as the orbit of 27, which takes 111 steps and peaks at 9232. It surveys computer verification up to 2^71 and the lower bounds on the length of a nontrivial cycle that it gives through the continued fraction of log base 2 of 3, the undecidability of generalized Collatz maps (Conway, FRACTRAN), and density theorems for stopping times, from Terras and Everett to Tao's almost bounded values. It explains the parity heuristic, a random walk for log x with drift log(sqrt(3)/2), and why it is not a proof: the parities along one orbit are deterministic. It describes the Lean statements of the conjecture, of its descent form and of the cycle problem written by the authors. It ends with the plan of the book, which does not prove the conjecture but reads it through martingales, optimization and Diophantine approximation.

## From the literature

The survey of Lagarias (1985) and the volume he edited (2010); the density theorem for stopping times of Terras (1976, 1979) and Everett (1977); the exponent log 3 over log 4 of Korec (1994); the x to the power 0.84 bound of Krasikov and Lagarias (2003); the almost bounded values theorem of Tao (2022); the undecidability results of Conway (1972, 1987) and Kurtz and Simon (2007); computer verification by Oliveira e Silva (2010) and Bařina (2021, 2025); cycle length bounds from the continued fraction of log base 2 of 3 by Crandall (1978), Eliahou (1993) and Hercher (2023); the parity heuristic of Crandall (1978) and Lagarias and Weiss (1992).

## New results

This chapter is an introduction and survey, and its labeled statements are the conjecture itself and a theorem of Tao (2022).

<!-- results:begin -->
No new results of this chapter are catalogued.
<!-- results:end -->

## Key terms

3x+1 problem, 3n+1 problem, Collatz problem, Collatz conjecture, Syracuse problem, Kakutani's problem, Ulam's problem, Hasse's algorithm, standard map, accelerated map, orbit, trajectory, total stopping time, delay record, divergent orbit, nontrivial cycle, verified range, logarithmic density, FRACTRAN, Lean, proof assistant, machine checked proof, Erdős, Erdos, Bařina, Barina
