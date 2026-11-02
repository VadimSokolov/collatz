Chapter 20, "What the Book Missed", for On the Collatz Conjecture (Polson and Sokolov)
Third draft, 27 October 2026.

chapter20_what_we_missed.tex   standalone chapter (book trim 7 x 10 in, Palatino); to merge,
                               drop the preamble and input the body after Chapter 17
chapter20_what_we_missed.pdf   compiled chapter (17 pages)
fig_ch20_cycle.pdf             Figure 20.1 (Christoffel cycle, constants along the ladder)
fig_ch20_phi.pdf               Figure 20.2 (Phi/B_r for uniform random admissible words)

code/
  ladder.py          continued fraction of log2 3 and the ladder of rungs
  br_sum.c           exact 128 bit binning of {t theta}; rigorous brackets for B_r up to 1.4e11
                     build: gcc -O3 -DHI=0x95c01a39fbd6879fULL -DLO=0xa00b120a068badd1ULL br_sum.c -o br_sum -lm
                     run:   ./br_sum 72057431991 137528045312   (about 4 minutes)
  sharp.py           Christoffel rational cycles, Table 20.1
  check_random.py    exact rational test of the minimum inequality on 3000 random cycles
  check_max.py       exact rational test of the maximum inequality
  thresholds.py      Table 20.2 (reads data/br_rungs.txt for the large rungs)
  admissible.py      exact DP count and uniform sampling of admissible words; "python3 admissible.py fig"
                     regenerates the samples in data/ (about 5 minutes)
  fig_phi.py         Figure 20.2 from the saved samples
  bandcount.py       exact counts of admissible words by time in the band, Table 20.3;
                     "python3 bandcount.py 2301" takes about 20 minutes; per rung results in data/
  mod8.py            Table 20.4
  persist.py         persistence constants, Proposition 20.12
  fig_ch20.py        Figure 20.1
data/                output of br_sum for the large rungs; saved Phi/B_r samples; band counts

Open computations for Vadim:
  Problem 20.2  branch and bound over words from the minimum: verification range at which
                q23 = 137,528,045,312 is excluded (Theorem 20.3 gives 2^74.96)
  Problem 20.4  count admissible words with Phi >= c B_r (bandcount.py with Phi discretised
                in place of n_t; or a biased split sampler for larger r), and the rate rho(c)
