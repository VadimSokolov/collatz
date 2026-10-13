# Chapter 5. Majorization and Momentum

<!-- cite:begin -->
From *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions* by Nick Polson and Vadim Sokolov (2026), chapter 5. Cite as Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1, https://doi.org/10.5281/zenodo.23331100. Book page: https://vsokolov.org/books/collatz/. Source: https://github.com/VadimSokolov/collatz.
<!-- cite:end -->

This chapter turns from descent to momentum methods in optimization and asks which analogies between momentum and the Collatz orbit mean something. It reviews Polyak's heavy ball method and Nesterov's accelerated gradient method on strongly convex quadratics, including the optimal constant momentum coefficient for condition number kappa, and the continuous time view in which momentum is a damped oscillator certified by an energy with a kinetic term. The Collatz orbit is read as momentum injection by the odd step followed by correction by the halvings, with the orbit of 27 as an example. The ratio of the largest stride ascent to the mean stride descent is about 1.776, a ratio of increments and not a condition number, and mixing constants defined at different levels produces spurious agreements. The number beta equal to log base 2 of 3 minus 1, about 0.585, appears as an ascent to descent ratio, the excess over one of the Cesàro threshold for halvings per odd step, a rewriting of drift constants, and a ratio of Kullback and Leibler divergences. Parrondo's paradox illustrates state dependent switching, and beta becomes the momentum coefficient of Switching MM.

## From the literature

The heavy ball method of Polyak (1964, 1987). The accelerated gradient method of Nesterov (1983, 2004). The optimal constant momentum on quadratics, classical, see Nesterov (2004). The failure of heavy ball to converge globally, from Lessard, Recht and Packard (2016). The differential equation of Su, Boyd and Candès (2016) and the high resolution equations of Shi, Du, Jordan and Su (2022). Lyapunov analysis of acceleration by Wilson, Recht and Jordan (2021) and adaptive restart by O'Donoghue and Candès (2015). Parrondo's games as published by Harmer and Abbott (1999), with Parrondo, Harmer and Abbott (2000) and the review of Harmer and Abbott (2002). An identity of Tao (2022) on the residue chain.

## New results

The one labeled proposition of the chapter, the optimal constant momentum on quadratics, is classical, and the constants and identities are elementary.

<!-- results:begin -->
No new results of this chapter are catalogued.
<!-- results:end -->

## Key terms

momentum, heavy ball, Nesterov acceleration, accelerated gradient, condition number, optimal momentum coefficient, damped oscillator, momentum injection, stride constants, Cesàro threshold, Cesaro, Kullback and Leibler divergence, Parrondo's paradox, Switching MM.
