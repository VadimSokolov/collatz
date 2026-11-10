# On the Collatz Conjecture <small>Martingales, Momentum and Continued Fractions</small>

By Nick Polson and [Vadim Sokolov](https://orcid.org/0000-0002-6618-2965).

[![PDF](https://vsokolov.org/books/icon-pdf.svg)](https://github.com/VadimSokolov/collatz/releases/download/v1.4.0/collatz.pdf "PDF") [![Source on GitHub](https://vsokolov.org/books/icon-github.svg)](https://github.com/VadimSokolov/collatz/tree/v1.4.0 "Source on GitHub") [![Paperback, available at Amazon](https://vsokolov.org/books/available-at-amazon.png)](https://www.amazon.com/dp/B0HN8TS5N8 "Paperback on Amazon")

## Description

The Collatz conjecture says that, from any positive integer, repeatedly halving an even number and sending an odd number n to 3n+1 always reaches 1. This book explains why the problem has resisted proof.

The conjecture splits into a divergence half, a law of large numbers along a single deterministic orbit, and a cycle half, a question about rational approximation of log2 3 that Diophantine information alone cannot finish. The book proves that no Lyapunov function increasing in the number certifies descent over any fixed number of steps, and that no finite residue sieve certifies descent for every integer. Infinitely many cycle lengths survive the cycle inequality, a point Halbeisen and Hungerbühler made earlier in another form; a sharp minimum inequality for rational cycles, nearly attained by Christoffel cycles, measures what a property of integers must add to exclude them. A seventeen state Turing machine that beeps forever exactly when the conjecture holds shows what one more quantifier costs. A tail criterion reduces the conjecture, given Lau's density theorem, to a polynomial bound.

Switching MM, a method that applies the momentum log2 3 minus 1 on alternate steps, has an explicit linear rate at step 1/L on smooth strongly convex functions, and on quadratics with additive gradient noise, at steps up to 1/L, its noise variance stays below 5/3 times that of gradient descent for every momentum coefficient below one. By the gain weight theorem, on quadratics with small steps no stable linear periodic first order method has a lower slow mode noise floor than gradient descent at matched speed. Beyond quadratics, at small steps, it runs at a fixed multiple of gradient descent's temperature.

Every result is proved in full or marked as conditional, and the Python scripts behind every computed number, figure and table are in the repository.

## How to cite

Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.24113412.

```bibtex
@book{polsonsokolov2026collatz,
  author = {Polson, Nick and Sokolov, Vadim},
  title = {{On the Collatz Conjecture: Martingales, Momentum and Continued Fractions}},
  year = {2026},
  publisher = {Independently published},
  doi = {10.5281/zenodo.24113412},
  url = {https://vsokolov.org/books/collatz/}
}
```

## Main results

- [No fixed horizon monotone certificate for Collatz descent](https://vsokolov.org/books/collatz/results.html#no_fixed_horizon_monotone_certificate) (Corollary 4.15, page 76): For every fixed number J of strides, no strictly increasing function of x certifies descent over J Syracuse strides, because the odd integers congruent to minus 1 modulo 2 to the J+1 ascend for J strides in a row.
- [Noise bound for Switching MM uniform in the normalized step](https://vsokolov.org/books/collatz/results.html#smm_variance_bound_uniform_in_normalized_step) (Proposition 8.7, page 128): On a quadratic with additive gradient noise, on every eigendirection with normalized step up to one, the stationary variance of Switching MM lies between that of gradient descent and 5/3 times it, uniformly over momentum coefficients below one, whereas every step momentum amplifies slow mode noise by a factor that grows without bound as the coefficient tends to one.
- [Gain weight theorem for the noise floor of periodic first order methods](https://vsokolov.org/books/collatz/results.html#gain_weight_theorem) (Theorem 8.12, page 133): For any linear periodic first order method on a quadratic that is stable at zero normalized step and has positive total gain, the slow mode noise floor at matched speed relative to gradient descent is a Cauchy and Schwarz ratio of the step gains, at least one, and equal to one only when all gains are equal.
- [Exact linear rate of Switching MM on smooth strongly convex functions](https://vsokolov.org/books/collatz/results.html#smm_linear_rate_on_smooth_strongly_convex_functions) (Theorem 9.9, page 169): For strongly convex and L smooth functions, two steps of Switching MM at step one over L contract the distance to the minimizer by the product of the classical one step factors of a long and a short gradient step, which explicit functions show is the exact worst case over the class, and the function gap by its square, using the one step function value bound of Taylor, Hendrickx and Glineur (2018).
- [Diophantine insufficiency, infinitely many cycle lengths survive any verification range](https://vsokolov.org/books/collatz/results.html#diophantine_insufficiency_cycle_lengths) (Proposition 11.20, page 210): However far the verification of the Collatz conjecture is pushed, infinitely many cycle lengths survive the cycle inequality, so information about the approximation of log2 3 alone cannot exclude all nontrivial cycles, a point Halbeisen and Hungerbühler (1997) made earlier in another form.
- [Residue certificates cannot prove descent for all integers (the exceptional 2 adic set)](https://vsokolov.org/books/collatz/results.html#residue_certificates_cannot_prove_descent_for_all_integers) (Proposition 15.1, page 273): No residue certificate of finite level can certify descent for every integer, because every residue class it leaves uncovered contains infinitely many positive integers, while the 2 adic integers with infinite coefficient stopping time form a closed null set that contains a point of each known negative cycle.
- [Sharp minimum inequality for cycles of the Syracuse map](https://vsokolov.org/books/collatz/results.html#sharp_minimum_inequality) (Theorem 18.3, page 309): The smallest element of a cycle is below a bound that depends only on the number of odd elements and the closing defect, with constant 1/(2 log 2) in place of 1, uniformly in the word.
- [Christoffel cycles attain the sharp minimum inequality](https://vsokolov.org/books/collatz/results.html#christoffel_cycles_are_extremal) (Theorem 18.5, page 310): The positive rational cycle of the lower Christoffel word of slope s/r has 3 M_w Lambda between 2 to the minus 2 delta times beta_r = r/(2 log 2) + 1/4 + O(1/r) and the sharp bound B_r(delta), and along the rungs the ratio to the sharp bound tends to 1, so no inequality in r, s and the minimum alone can improve the constant 1/(2 log 2).

All 132 results of the book (85 proved, 2 computed, 5 conditional, 40 open) are on the page [Results by chapter](https://vsokolov.org/books/collatz/results.html), each with its statement, status and number.

***

<small>License: Text, figures and PDF: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: [MIT](https://opensource.org/license/mit).</small>

<small>This page is generated from the book's repository ([main](https://github.com/VadimSokolov/collatz/tree/main)). The numbers and pages are those of the PDF, [version 1.4.0](https://github.com/VadimSokolov/collatz/releases/tag/v1.4.0).</small>

<!-- Generated by bookrepo.py from book.yaml and results.yaml. Edit those files, not this one. -->
