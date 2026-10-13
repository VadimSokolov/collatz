# Chapter 7. Switching MM

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 7. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23331100. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

This chapter asks whether the momentum coefficient of a first order method can be fixed in advance, without knowing the condition number, and what such a default achieves. It derives the coefficient log(3/2)/log 2, which equals log2 3 minus 1 and is about 0.585, in two ways: from the Collatz map, as the ratio of the ascent of an odd step to the descent of a halving (the log displacement identity of the Syracuse map), and from elementary information theory, as a ratio of Kullback and Leibler divergences of nested uniform laws. The resulting method, Switching MM (SMM), applies Nesterov extrapolation only on odd iterations and a plain gradient step on even ones, with step size one over L. The main results are exact. For every differentiable function, two SMM steps from an even index are a gradient step of length (1+beta) alpha followed by one of length alpha. On strongly convex quadratics SMM is therefore a period two cyclic step size method with an explicit stability range and linear rate, and in the slow modes it acts like gradient descent with a step enlarged by (2+beta)/2, so it is not an accelerated method. A defect bound covers nonquadratic objectives with Lipschitz Hessian. The chapter states that the arithmetic bears nothing on the Collatz conjecture.

## From the literature

Heavy ball momentum (Polyak 1964, 1987) and accelerated gradient (Nesterov 1983, 2004); the limit cycle of tuned heavy ball (Lessard, Recht and Packard 2016); the Adam default momentum 0.9 (Kingma and Ba 2015); differential equation views of acceleration (Su, Boyd and Candès 2016; Shi, Du, Jordan and Su 2022); cyclic and Chebyshev step sizes (Richardson 1911; Young 1953; Golub and Varga 1961) and the silver stepsize schedule (Altschuler and Parrilo 2025); MM algorithms (Lange, Hunter and Yang 2000; Hunter and Lange 2004); switched linear systems (Liberzon 2003).

## New results

<!-- results:begin -->
- [The extrapolated point of Switching MM is a long gradient step](../RESULTS.md#extrapolated_point_is_a_long_gradient_step) (Lemma 7.5, page 117; proved): For every differentiable objective, two consecutive steps of Switching MM from an even index are exactly a long gradient step, longer by the factor one plus the momentum coefficient, followed by a plain gradient step. New in this book.
- [Two step equivalence of Switching MM with cyclic gradient steps on quadratics](../RESULTS.md#smm_two_step_equivalence_on_quadratics) (Theorem 7.6, page 118; proved): On a strongly convex quadratic, two steps of Switching MM from an even index equal a plain gradient step followed by a longer one, so the method is a cyclic step size method and every two steps its error is multiplied by an explicit degree two polynomial of the Hessian. New in this book.
- [Stability range and linear convergence of Switching MM on quadratics](../RESULTS.md#smm_stability_range_on_quadratics) (Corollary 7.8, page 119; proved): Switching MM converges linearly on every strongly convex quadratic for any momentum coefficient below one and step at most one over L, and it is stable exactly for normalized steps below an explicit limit. New in this book.
- [Slow mode effective step of Switching MM on quadratics](../RESULTS.md#smm_slow_mode_effective_step) (Corollary 7.9, page 120; proved): In the slow directions Switching MM behaves like gradient descent with a step enlarged by the factor two plus beta over two, about 1.29 at the Collatz momentum, so it is not an accelerated method. New in this book.
- [Third order defect of the cyclic step representation of Switching MM](../RESULTS.md#smm_defect_of_cyclic_representation) (Proposition 7.10, page 122; proved): For objectives with a Lipschitz Hessian, the failure of the odd iterates of Switching MM to follow a cyclic gradient path is of third order in the step size, and it vanishes for quadratics. New in this book.
<!-- results:end -->

## Key terms

Switching MM, SMM, momentum coefficient, log2 3 minus 1, Nesterov accelerated gradient, NAG, heavy ball, cyclic step sizes, Chebyshev step sizes, silver stepsize schedule, effective step multiplier, stability range, Syracuse map, Kullback and Leibler divergence, KL divergence, MM algorithm, switched linear system, low tuning optimizer
