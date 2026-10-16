# Chapter 3. Stochastic Models and Martingale Duality

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 3. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23421852. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

This chapter asks what the probabilistic heuristic for the 3x+1 problem predicts and what martingale theory can prove about it. It defines three stochastic models of the Collatz iteration: the Bernoulli fair coin model for the accelerated map, the Geometric stride model for the Syracuse map, and a renewal (Markov chain) model for the standard map, where an odd step is always followed by a halving. In the Bernoulli model the iterates form a multiplicative martingale, because one half times one half plus one half times three halves equals one, while their logarithms form a supermartingale with drift log(sqrt(3)/2), about minus 0.1438, which is the Jensen gap (martingale duality). Wald's identity gives the slope 1/|drift|, about 6.952, of the expected passage time, the model's prediction for total stopping times, classical first passage results give its dispersion, Ville's maximal inequality bounds an excursion by a factor K with probability at most 1/K, and the model is absorbed almost surely. In the renewal model for the standard map the parity autocorrelation is exactly (minus one half) to the power h. The stopped martingale is not uniformly integrable, which shows why it cannot be used to settle absorption, even in the model. The predictions are compared with orbits, a census of all n below 2^30 places most of the ten percent variance deficit of total stopping times at the end of the orbit, and the chapter asks what a certificate for individual orbits would need.

## From the literature

The fair coin picture of Terras (1976) and Everett (1977); the parity heuristic of Crandall (1978) and its drift log(sqrt(3)/2), recorded by Lagarias (1985); the random walk models of Lagarias and Weiss (1992); the structure theorem of Kontorovich and Sinai (2002) and the survey of Kontorovich and Lagarias (2010); the first passage asymptotics of Heyde (1967) and Siegmund (1969) and the normal limit for this model of Borovkov and Pfeifer (2001); Wald's identity (Wald 1944); the maximal inequality for nonnegative supermartingales of Ville (1939); the Lundberg inequality of risk theory, as treated by Asmussen and Albrecher (2010).

## New results

The Bernoulli drift is classical; the authors' statements are elementary closed forms for the renewal model and the failure of uniform integrability.

<!-- results:begin -->
- [Martingale structure of the renewal model for the standard Collatz map](../RESULTS.md#martingale_structure_of_the_standard_map) (Proposition 3.7, page 48; proved): In the renewal model of the standard Collatz map, where every tripling is followed by a halving, the iterates are not a martingale, but they become one once each coming tripling is booked in advance. New in this book.
- [Exact parity autocorrelation function of the renewal model for the standard Collatz map](../RESULTS.md#parity_autocorrelation_function) (Proposition 3.12, page 56; proved): In the renewal model of the standard Collatz map one step in three is odd and the parity autocorrelation at lag h is minus one half raised to the power h. New in this book.
- [Long run variance of the log path in the renewal model for the standard Collatz map](../RESULTS.md#long_run_variance_of_the_log_path) (Corollary 3.14, page 58; proved): The negative parity correlations of the renewal model cut the long run variance of the log path to one third of its value for independent parities with the same frequency, one in three, of odd steps. New in this book.
- [The stopped martingale of the Bernoulli model is not uniformly integrable](../RESULTS.md#failure_of_uniform_integrability) (Proposition 3.17, page 61; proved): The martingale of the Bernoulli model, stopped when the iterate first drops to 1 or below, is not uniformly integrable, because the unabsorbed paths keep expected value at least the starting value minus 1 while their probability tends to zero. New in this book.
<!-- results:end -->

## Key terms

stochastic model, fair coin heuristic, parity heuristic, Bernoulli model, stride model, renewal model, martingale duality, supermartingale, Jensen gap, Wald identity, Lundberg inequality, Ville inequality, first passage time, parity autocorrelation, long run variance, variance deficit, absorption, uniform integrability
