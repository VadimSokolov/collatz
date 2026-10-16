# Results by chapter <small>On the Collatz Conjecture</small>

The results of [On the Collatz Conjecture: Martingales, Momentum and Continued Fractions](https://vsokolov.org/books/collatz/) that are new in the book, that come from the authors' own papers, or that correct one of those papers, as the book states them: 66 entries (48 proved, 3 conditional, 15 open). Results of other authors that the book presents are credited in the [chapter summaries](https://github.com/VadimSokolov/collatz/tree/main/summaries).

Status: proved; conditional (proved assuming the named statements, which the book does not prove); open (stated as an open question); counterexample (a proved counterexample); withdrawn (found wrong after release; the entry says what was wrong and which version corrects it).

## Chapter 3. Stochastic Models and Martingale Duality

<a id="martingale_structure_of_the_standard_map"></a>

### Martingale structure of the renewal model for the standard Collatz map

Proposition 3.7, page 48 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $(B_n)_{n\ge0}$ be the Markov chain on $\{0,1\}$ with $\mathbb{P}(B_{n+1}=0\mid B_n=1)=1$ and $\mathbb{P}(B_{n+1}=0\mid B_n=0)=\mathbb{P}(B_{n+1}=1\mid B_n=0)=\tfrac12$, let $x_0>0$, and put $x_{n+1}=3x_n$ if $B_n=1$ and $x_{n+1}=x_n/2$ if $B_n=0$, with $\mathcal{F}_n=\sigma(B_0,\ldots,B_n)$ (the renewal model for the standard Collatz map, which sends an even $m$ to $m/2$ and an odd $m$ to $3m+1$, where $B_n$ models the parity of the $n$th iterate). Write $\delta^{*}=\log(\sqrt3/2)$. Then (i) the stationary law of the chain satisfies $\mathbb{P}(B_n=1)=\tfrac13$, and the fraction of odd steps (indices $i$ with $B_i=1$) among the first $n$ steps converges to $\tfrac13$ almost surely from any initial state; (ii) $(x_n)$ is not a martingale, since $\mathbb{E}[x_{n+1}\mid\mathcal{F}_n]=3x_n$ on $\{B_n=1\}$, but $M_n=3^{B_n}x_n$ is a martingale; (iii) $(\log x_n)$ is not a supermartingale, since its increment is $+\log3$ on $\{B_n=1\}$, but $\log M_n$ is a supermartingale with conditional drift $0$ on $\{B_n=1\}$ and $\delta^{*}$ on $\{B_n=0\}$; (iv) almost surely $\frac1n\log x_n\to\tfrac13\log3 - \tfrac23\log2=\tfrac13\log\tfrac34=\tfrac23\delta^{*}\approx -0.0959$; (v) if $\tau_1<\tau_2<\cdots$ are the successive times at which $B_n=1$, then $X_j=x_{\tau_j}$ satisfies $X_{j+1}=3X_j2^{-K_{j+1}}$, where $K_{j+1}=\tau_{j+1} - \tau_j - 1$ are independent with the Geometric$(\tfrac12)$ law $\mathbb{P}(K=k)=2^{-k}$ for $k\ge1$, and $(X_j)$ is a martingale because $\mathbb{E}[3\cdot2^{-K}]=3\sum_{k\ge1}4^{-k}=1$.

Source: [`book/chapters/ch03_stochastic.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch03_stochastic.tex), label `prop:ch3:standard`.

<a id="parity_autocorrelation_function"></a>

### Exact parity autocorrelation function of the renewal model for the standard Collatz map

Proposition 3.12, page 56 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $(B_n)_{n\ge0}$ be the Markov chain on $\{0,1\}$ with $\mathbb{P}(B_{n+1}=0\mid B_n=1)=1$ and $\mathbb{P}(B_{n+1}=0\mid B_n=0)=\mathbb{P}(B_{n+1}=1\mid B_n=0)=\tfrac12$ (the renewal model for the standard Collatz map, where $B_n$ models the parity of the $n$th iterate), started from its stationary law. Then $(B_n)$ is stationary with $\mathbb{P}(B_n=1)=\frac{1}{1+\mathbb{E}[K]}=\frac13$ and $\operatorname{Var}(B_n)=\frac29$, where $K$, the number of zeros between successive ones of the chain, has the Geometric$(\tfrac12)$ law $\mathbb{P}(K=k)=2^{-k}$ for $k\ge1$ (so $\mathbb{E}[K]=2$), and its autocorrelation function is $\rho(h)=\operatorname{Corr}(B_n,B_{n+h})=\bigl(-\tfrac12\bigr)^{h}$ for $h\ge0$.

Source: [`book/chapters/ch03_stochastic.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch03_stochastic.tex), label `thm:autocorr`.

<a id="long_run_variance_of_the_log_path"></a>

### Long run variance of the log path in the renewal model for the standard Collatz map

Corollary 3.14, page 58 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $(B_n)_{n\ge0}$ be the Markov chain on $\{0,1\}$ with $\mathbb{P}(B_{n+1}=0\mid B_n=1)=1$ and $\mathbb{P}(B_{n+1}=0\mid B_n=0)=\mathbb{P}(B_{n+1}=1\mid B_n=0)=\tfrac12$, started from any initial state, let $x_0>0$, and put $x_{n+1}=3x_n$ if $B_n=1$ and $x_{n+1}=x_n/2$ if $B_n=0$ (the renewal model for the standard Collatz map). Let $\rho(h)=(-\tfrac12)^h$ be the autocorrelation function of $(B_n)$. The autocorrelation factor is $\mathcal{A}:=1+2\sum_{h\ge1}\rho(h)=1+2\cdot\frac{-1/2}{1+1/2}=\frac13$, and the variance of the log path satisfies $\frac{\operatorname{Var}(\log x_n)}{n}\longrightarrow\mathcal{A}\cdot\frac29(\log6)^2=\frac{2}{27}(\log6)^2\approx0.2378$ as $n\to\infty$.

Source: [`book/chapters/ch03_stochastic.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch03_stochastic.tex), label `cor:longrun`.

<a id="failure_of_uniform_integrability"></a>

### The stopped martingale of the Bernoulli model is not uniformly integrable

Proposition 3.17, page 61 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $x_0>1$ and let $\xi_1,\xi_2,\ldots$ be independent random variables with $\mathbb{P}(\xi_i=\tfrac12)=\mathbb{P}(\xi_i=\tfrac32)=\tfrac12$, and put $x_n=x_0\xi_1\cdots\xi_n$ (the Bernoulli model for the accelerated Collatz map) and $\tau=\inf\{n\ge0:x_n\le1\}$. Then $\tau<\infty$ almost surely, the stopped process $(x_{n\wedge\tau})_{n\ge0}$ is a nonnegative martingale converging almost surely to $x_\tau\le1$, and it is not uniformly integrable. More precisely, for every $n$, $\mathbb{E}\bigl[x_n\mathbf{1}_{\{\tau>n\}}\bigr]\ge x_0 - 1$, while $\mathbb{P}(\tau>n)\to0$ as $n\to\infty$.

Source: [`book/chapters/ch03_stochastic.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch03_stochastic.tex), label `prop:uimartingale`.

## Chapter 4. Lyapunov Functions and the Limits of Descent

<a id="stride_model_passage_time"></a>

### Passage time of the stride model with overshoot

Proposition 4.5, page 70 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Consider the stride model of the Collatz map. Let $K_1,K_2,\dots$ be independent random variables with $\mathbb{P}(K_j=k)=2^{-k}$ for $k\ge1$ (the Geometric$(\tfrac12)$ law), fix $x_0>1$, and put $x_j=x_{j-1}\cdot3\cdot2^{-K_j}$, so that $\log x_j=\log x_0+\Delta_1+\dots+\Delta_j$ with $\Delta_j=\log3-K_j\log2$. Let $\tau$ be the first $j$ with $\log x_j\le0$, and let $O=-\log x_\tau\ge0$ be the overshoot below $0$ at time $\tau$. Then $\log2\le\mathbb{E}[O]<2\log2$ and \[\frac{\log x_0+\log2}{|\log(3/4)|}\le\mathbb{E}_{x_0}[\tau]<\frac{\log x_0+2\log2}{|\log(3/4)|}.\]

Source: [`book/chapters/ch04_lyapunov.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch04_lyapunov.tex), label `prop:ch4:passage`.

<a id="stride_descent"></a>

### Stride descent theorem for the Syracuse map

Theorem 4.7, page 71 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $S$ be the Syracuse map on the odd positive integers, $S(x)=(3x+1)/2^{k}$ with $k=v_2(3x+1)$, where $v_2(m)$ is the exponent of the largest power of $2$ dividing $m$. Let $x\ge3$ be odd and $k=v_2(3x+1)$. (i) $\log S(x)-\log x=\log(3+1/x)-k\log2$. (ii) If $k\ge2$, then $S(x)<x$, $\log S(x)-\log x\le\log\bigl(\tfrac34+\tfrac1{4x}\bigr)<0$, and $\log S(x)-\log x\le\log(7/9)$, with equality exactly at $x=9$. (iii) If $k=1$, which happens exactly when $x\equiv3\pmod4$, then $S(x)=(3x+1)/2>x$ and $\log\tfrac32<\log S(x)-\log x=\log\bigl(\tfrac32+\tfrac1{2x}\bigr)\le\log\tfrac53$, with equality on the right exactly at $x=3$.

Source: [`book/chapters/ch04_lyapunov.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch04_lyapunov.tex), label `thm:stride`.

<a id="no_monotone_single_stride_lyapunov"></a>

### No monotone strong Lyapunov function for a single Syracuse stride

Theorem 4.8, page 73 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $S$ be the Syracuse map on the odd positive integers, $S(x)=(3x+1)/2^{v_2(3x+1)}$, where $v_2(m)$ is the exponent of the largest power of $2$ dividing $m$. Let $W$ be a strictly increasing function from the odd positive integers to $(0,\infty)$, and let $\rho\le1$. Then $W(S(x))>\rho\,W(x)$ for every $x\equiv3\pmod4$. In particular no strictly increasing $W$ is a strict Lyapunov function for $(S,\{1\})$ (one with $W(S(x))<W(x)$ for every odd $x\ne1$) or a strong Lyapunov function for $(S,\{1\})$ (one with $W>0$ and $W(S(x))\le\rho W(x)$ for every odd $x\ne1$, for some constant $\rho<1$), and the set of odd integers on which $W(S(x))\le\rho W(x)$ fails contains the residue class $3\bmod4$, which has natural density $\tfrac12$ among the odd integers; for $\rho=1$ the set on which it fails is exactly this class.

Source: [`book/chapters/ch04_lyapunov.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch04_lyapunov.tex), label `thm:nolyapunov`.

<a id="exact_two_stride_behaviour"></a>

### Exact behaviour of two Syracuse strides

Proposition 4.12, page 74 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $S$ be the Syracuse map on the odd positive integers, $S(x)=(3x+1)/2^{v_2(3x+1)}$, where $v_2(m)$ is the exponent of the largest power of $2$ dividing $m$. For odd $x$ write $k_1=v_2(3x+1)$, $k_2=v_2(3S(x)+1)$ and $S^{(2)}(x)=S(S(x))$. Let $x\ge3$ be odd. (i) $S^{(2)}(x)<x$ if and only if $k_1+k_2\ge4$. (ii) The odd integers with a two stride descent, that is with $S^{(2)}(x)<x$, have natural density exactly $\tfrac12$ among the odd integers. (iii) $\max_{x\ge3\text{ odd}}\bigl(\log S^{(2)}(x)-\log x\bigr)=\log\tfrac{17}{7}\approx0.8873$, attained only at $x=7$, where $S(7)=11$ and $S^{(2)}(7)=17$.

Source: [`book/chapters/ch04_lyapunov.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch04_lyapunov.tex), label `prop:ch4:twostride`.

<a id="no_fixed_horizon_monotone_certificate"></a>

### No fixed horizon monotone certificate for Collatz descent

Corollary 4.15, page 78 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $S$ be the Syracuse map on the odd positive integers, $S(x)=(3x+1)/2^{v_2(3x+1)}$, where $v_2(m)$ is the exponent of the largest power of $2$ dividing $m$, and let $S^{(J)}$ be its $J$th iterate. Fix $J\ge1$, $\rho\le1$ and a strictly increasing function $W$ from the odd positive integers to $(0,\infty)$. Then $W(S^{(J)}(x))>\rho\,W(x)$ for every odd $x\equiv-1\pmod{2^{J+1}}$, that is, for every $x=2^{J+1}q-1$ with $q\ge1$. This set is infinite and has natural density $2^{-J}$ among the odd integers. Moreover, along it, \[\frac{S^{(J)}(x)}{x}\longrightarrow\Bigl(\frac32\Bigr)^{J}\qquad(q\to\infty).\]

Source: [`book/chapters/ch04_lyapunov.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch04_lyapunov.tex), label `cor:ch4:multistride`.

<a id="linear_majorizer_cannot_certify_descent"></a>

### The linear majorizer of log cannot certify Collatz descent

Proposition 4.19, page 79 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $T_0$ be the standard Collatz map, $T_0(x)=x/2$ for even $x$ and $T_0(x)=3x+1$ for odd $x$, and let $g(y\mid x)=\log x+(y-x)/x$ be the tangent line majorizer of $\log$ at $x$, so that the surrogate increment $g(T_0(x)\mid x)-g(x\mid x)$ equals $-\tfrac12$ for even $x$ and $2+1/x$ for odd $x$. Under the stochastic model for the standard map, in which the parities of successive steps form a stationary renewal chain $(B_n)$ on $\{0,1\}$ (with $B_n=1$ for an odd step; an odd step is always followed by an even step, and an even step is followed by an odd or an even step with probability $\tfrac12$ each), so that a step is odd with stationary probability $\tfrac13$, the mean surrogate increment per step of $T_0$ is $\tfrac13\cdot2+\tfrac23\cdot\bigl(-\tfrac12\bigr)=+\tfrac13>0$, up to the $O(1/x)$ term, whereas the mean true increment per step is $\tfrac13\log3-\tfrac23\log2=\tfrac13\log(3/4)\approx-0.0959<0$. Per Syracuse stride, that is one odd step followed by $k$ even steps, in the stride model with $k$ distributed as Geometric$(\tfrac12)$ on $\{1,2,\dots\}$, the per step surrogate increments summed over the stride, with the majorizer recentred at every step, have mean $+1$ up to the $O(1/x)$ term, and the mean true increment is $\log(3/4)\approx-0.2877$.

Source: [`book/chapters/ch04_lyapunov.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch04_lyapunov.tex), label `prop:mmfails`.

<a id="power_majorizers_of_log"></a>

### Power majorizers of log and the Wald root of the stride model

Proposition 4.20, page 81 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Work in the stride model of the Collatz map, in which the log increment of one stride is $\Delta=\log3-K\log2$ with $K$ distributed as Geometric$(\tfrac12)$ on $\{1,2,\dots\}$, and write $M(p)=\mathbb{E}[e^{p\Delta}]=3^{p}/(2^{1+p}-1)$. For $0<p\le1$ let $g_p(y\mid x)=\log x+\frac{(y/x)^{p}-1}{p}$; each $g_p$ dominates $\log y$ and is tangent to it at $y=x$, and $p=1$ is the tangent line (linear) majorizer of $\log$. Then the mean stride increment of the surrogate $g_p$, ignoring the $1/x$ correction, is $(M(p)-1)/p$, and: (i) $M(0)=M(1)=1$, and $M(p)<1$ for $0<p<1$; (ii) $p\mapsto(M(p)-1)/p$ is increasing on $(0,1]$, tends to $\mathbb{E}[\Delta]=\log(3/4)$ as $p\downarrow0$, and equals $0$ at $p=1$. In particular every power majorizer with $0<p<1$ has strictly negative mean surrogate increment, and the linear majorizer $p=1$ sits exactly on the boundary.

Source: [`book/chapters/ch04_lyapunov.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch04_lyapunov.tex), label `prop:ch4:power`.

## Chapter 6. Majorization Pressed Harder

<a id="no_residue_potential_for_any_modulus"></a>

### No residue potential for the Collatz map at any modulus

Theorem 6.3, page 101 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $T$ be the accelerated Collatz map, $T(x)=x/2$ for even $x$ and $T(x)=(3x+1)/2$ for odd $x$. Fix $d\ge1$ and let $G_d$ be the directed graph on the residues modulo $2^d$ in which the residue $a$ has an edge to each residue $a'$ with $a'\equiv T(a)\pmod{2^{d-1}}$ (the two lifts of $T(a)$), with weight $w(a)=\log(3/2)$ if $a$ is odd and $w(a)=-\log2$ if $a$ is even. (i) There are a function $\varphi$ on $\mathbb{Z}/2^d\mathbb{Z}$ and a bound $X$ such that $V(T(x))<V(x)$ for all integers $x\ge X$, where $V(x)=\log x+\varphi(x\bmod2^d)$, if and only if every cycle of $G_d$ has negative total weight. (ii) For every $d\ge1$ the graph $G_d$ contains cycles of positive weight: the loop at $-1$, of weight $\log(3/2)$, and the images of the cycles of $T$ on the negative integers, whose weights are $r\log3-s\log2>0$ for a cycle with $r$ odd steps among $s$ steps. The maximum cycle mean of $G_d$ is exactly $\log(3/2)$ for every $d$. (iii) Consequently no potential $\log x+\varphi(x\bmod2^d)$ with bounded $\varphi$ decreases at every step for all large $x$, for any $d$, and the same holds for potentials required to decrease over any fixed number $J$ of steps.

Source: [`book/chapters/ch05c_mmcollatz.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch05c_mmcollatz.tex), label `thm:ch5c:residue`.

<a id="surplus_along_chains"></a>

### Surplus of a pointwise potential along runs of trailing ones

Proposition 6.4, page 102 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $S$ be the Syracuse map on the odd positive integers, $S(x)=(3x+1)/2^{v_2(3x+1)}$, where $v_2(m)$ is the exponent of the largest power of $2$ dividing $m$, and let $V(x)=\log x+\varphi(x)$ for a real function $\varphi$ on the odd positive integers. Suppose $V(S(x))<V(x)$ for every odd $x\ge3$. Then for every $j\ge2$ and every odd $q\ge1$, \[\varphi(2^jq-1)-\varphi(2\cdot3^{j-1}q-1)\;>\;(j-1)\log\frac32 .\] In particular $\varphi$ is unbounded above or below, and $\sup\varphi-\inf\varphi\ge(j-1)\log(3/2)$ for every $j$.

Source: [`book/chapters/ch05c_mmcollatz.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch05c_mmcollatz.tex), label `prop:ch5c:surplus`.

<a id="run_length_potential"></a>

### Run length potential and its failure fraction

Proposition 6.6, page 104 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $S$ be the Syracuse map on the odd positive integers, $S(x)=(3x+1)/2^{k}$ with $k=v_2(3x+1)$, where $v_2(m)$ is the exponent of the largest power of $2$ dividing $m$, let $\varepsilon(x)=\log\bigl(1+\tfrac1{3x}\bigr)$, and for $c>\log(3/2)$ let $V_c(x)=\log x+c\,v_2(x+1)$. Then $V_c$ decreases at every ascending stride (a stride with $k=1$) from an odd $x>1/\bigl(2(e^c-3/2)\bigr)$, hence at every ascending stride when $c>\log(5/3)$ and at all but finitely many when $c\le\log(5/3)$. The descending strides ($k\ge2$) at which $V_c$ fails to decrease are exactly those $x\mapsto y=S(x)$, with $k$ halvings, that land on a $y$ with $v_2(y+1)-1\ge(k\log2-\log3-\varepsilon(x))/c$. Under the stride model, in which $k$ and $v_2(y+1)$ are independent Geometric$(\tfrac12)$ variables on $\{1,2,\dots\}$, the fraction of odd $x$ at which $V_c$ fails to decrease is \[p(c)=\sum_{k\ge2}2^{-k}\,2^{-\lceil(k\log2-\log3)/c\rceil},\] which is $0.2087$ at $c=1$ and decreases to $0.1431$ as $c\downarrow\log(3/2)$, and $V_c$ fails at a positive fraction of odd integers for every $c$.

Source: [`book/chapters/ch05c_mmcollatz.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch05c_mmcollatz.tex), label `prop:ch5c:runpot`.

## Chapter 7. Switching MM

<a id="extrapolated_point_is_a_long_gradient_step"></a>

### The extrapolated point of Switching MM is a long gradient step

Lemma 7.5, page 117 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f:\mathbb{R}^d\to\mathbb{R}$ be differentiable, let $\alpha>0$, and write $G_h(x)=x-h\nabla f(x)$ for a gradient step of length $h$. Switching MM with step $\alpha$ and momentum coefficient $\beta$ starts from $x_{-1}=x_0$, takes the plain gradient step $x_{t+1}=G_\alpha(x_t)$ at even $t$, and at odd $t$ takes the extrapolated step $x_{t+1}=G_\alpha(y_t)$ with $y_t=x_t+\beta(x_t-x_{t-1})$. Let $t$ be even. Then $x_{t+1}=G_\alpha(x_t)$, and the extrapolated point at the following odd iteration is $y_{t+1}=x_{t+1}+\beta(x_{t+1}-x_t)=x_t-(1+\beta)\alpha\nabla f(x_t)=G_{(1+\beta)\alpha}(x_t)$. Consequently $x_{t+2}=G_\alpha\bigl(G_{(1+\beta)\alpha}(x_t)\bigr)$, and in the form that exhibits the momentum, $x_{t+2}=x_{t+1}-\beta\alpha\nabla f(x_t)-\alpha\nabla f(y_{t+1})$.

Source: [`book/chapters/ch06_smm.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch06_smm.tex), label `lem:ch6:ystep`.

<a id="smm_two_step_equivalence_on_quadratics"></a>

### Two step equivalence of Switching MM with cyclic gradient steps on quadratics

Theorem 7.6, page 118 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f(x)=\tfrac12x^\top Ax-b^\top x$ with $A$ symmetric positive definite and minimizer $x^\ast=A^{-1}b$, let $\alpha>0$ and $\beta\ge0$, write $G_h(x)=x-h\nabla f(x)$, and let $(x_t)$ be generated by Switching MM with step $\alpha$ and momentum coefficient $\beta$ (it starts from $x_{-1}=x_0$, takes the plain gradient step $x_{t+1}=x_t-\alpha\nabla f(x_t)$ at even $t$, and at odd $t$ takes the extrapolated step $x_{t+1}=y_t-\alpha\nabla f(y_t)$ with $y_t=x_t+\beta(x_t-x_{t-1})$). For every even $t\ge0$, $x_{t+2}=G_{(1+\beta)\alpha}(x_{t+1})=G_{(1+\beta)\alpha}\bigl(G_\alpha(x_t)\bigr)$, that is, two steps of Switching MM from an even index are a gradient step of length $\alpha$ followed by a gradient step of length $(1+\beta)\alpha$, and every iterate of Switching MM is an iterate of this cyclic schedule. In terms of the error $e_t=x_t-x^\ast$, $e_{t+2}=(I-\alpha A)\bigl(I-(1+\beta)\alpha A\bigr)e_t$ for even $t$.

Source: [`book/chapters/ch06_smm.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch06_smm.tex), label `thm:smmquad`.

<a id="smm_stability_range_on_quadratics"></a>

### Stability range and linear convergence of Switching MM on quadratics

Corollary 7.8, page 119 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f(x)=\tfrac12x^\top Ax-b^\top x$ with $A$ symmetric positive definite whose largest eigenvalue is $L$, let $\alpha>0$ and $0\le\beta<1$, and let Switching MM with step $\alpha$ and momentum coefficient $\beta$ (it starts from $x_{-1}=x_0$, takes the plain gradient step $x_{t+1}=x_t-\alpha\nabla f(x_t)$ at even $t$, and at odd $t$ takes the extrapolated step $x_{t+1}=y_t-\alpha\nabla f(y_t)$ with $y_t=x_t+\beta(x_t-x_{t-1})$) have the two step multiplier $\rho(a)=(1-a)\bigl(1-(1+\beta)a\bigr)=1-(2+\beta)a+(1+\beta)a^2$ on an eigendirection with eigenvalue $\lambda$ and normalized step $a=\alpha\lambda$. Then for real $a$, $|\rho(a)|<1$ if and only if $0<a<(2+\beta)/(1+\beta)$. In particular Switching MM converges linearly on every strongly convex quadratic whenever $\alpha\le1/L$, and in fact whenever $\alpha L<(2+\beta)/(1+\beta)$. For all $a>0$, $\rho(a)\ge1-(2+\beta)^2/\bigl(4(1+\beta)\bigr)$, a minimum attained at $a=(2+\beta)/\bigl(2(1+\beta)\bigr)$.

Source: [`book/chapters/ch06_smm.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch06_smm.tex), label `cor:ch6:stability`.

<a id="smm_slow_mode_effective_step"></a>

### Slow mode effective step of Switching MM on quadratics

Corollary 7.9, page 120 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f(x)=\tfrac12x^\top Ax-b^\top x$ with $A$ symmetric positive definite, let $\alpha>0$ and $\beta\ge0$, and let Switching MM with step $\alpha$ and momentum coefficient $\beta$ (it starts from $x_{-1}=x_0$, takes the plain gradient step $x_{t+1}=x_t-\alpha\nabla f(x_t)$ at even $t$, and at odd $t$ takes the extrapolated step $x_{t+1}=y_t-\alpha\nabla f(y_t)$ with $y_t=x_t+\beta(x_t-x_{t-1})$) have the two step multiplier $\rho(a)=(1-a)\bigl(1-(1+\beta)a\bigr)$ on an eigendirection with eigenvalue $\lambda$ and normalized step $a=\alpha\lambda$ (that is, the error component along that eigenvector satisfies $e_{t+2}=\rho(a)e_t$ for even $t$). As $a\to0$, $\rho(a)=1-(2+\beta)a+O(a^2)$, and the per step contraction factor $|\rho(a)|^{1/2}$ equals $1-\tfrac{2+\beta}{2}a+O(a^2)$. Thus in the slow directions Switching MM behaves like gradient descent with the effective step $m_{\mathrm{SMM}}\,\alpha$, where $m_{\mathrm{SMM}}=\tfrac{2+\beta}{2}$, which is $1.2925$ at $\beta=\log_23-1$.

Source: [`book/chapters/ch06_smm.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch06_smm.tex), label `cor:ch6:slow`.

<a id="smm_defect_of_cyclic_representation"></a>

### Third order defect of the cyclic step representation of Switching MM

Proposition 7.10, page 122 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f:\mathbb{R}^d\to\mathbb{R}$ be twice differentiable with $\nabla^2f$ Lipschitz with constant $M$ (in the Euclidean and induced operator norms), let $\alpha>0$ and $\beta\ge0$, write $G_h(x)=x-h\nabla f(x)$, and let $(x_t)$ be generated by Switching MM with step $\alpha$ and momentum coefficient $\beta$ (it starts from $x_{-1}=x_0$, takes the plain gradient step $x_{t+1}=x_t-\alpha\nabla f(x_t)$ at even $t$, and at odd $t$ takes the extrapolated step $x_{t+1}=y_t-\alpha\nabla f(y_t)$ with $y_t=x_t+\beta(x_t-x_{t-1})$). For even $t$, with $h=x_{t+1}-x_t=-\alpha\nabla f(x_t)$, $x_{t+2}-G_{(1+\beta)\alpha}(x_{t+1})=\alpha\bigl[(1+\beta)\nabla f(x_{t+1})-\beta\nabla f(x_{t+1}-h)-\nabla f(x_{t+1}+\beta h)\bigr]$, and $\bigl\|x_{t+2}-G_{(1+\beta)\alpha}(x_{t+1})\bigr\|\le\tfrac{M}{2}\,\beta(1+\beta)\,\alpha^3\,\|\nabla f(x_t)\|^2$. The bracket vanishes identically when $\nabla f$ is affine.

Source: [`book/chapters/ch06_smm.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch06_smm.tex), label `prop:ch6:defect`.

## Chapter 8. Noise, Variance and Experiments

<a id="smm_stationary_variance_additive_noise"></a>

### Exact stationary variance of Switching MM under additive gradient noise

Theorem 8.4, page 128 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f(x)=\tfrac12x^\top Ax-b^\top x$ with $A$ symmetric positive definite and minimizer $x^\ast$, and suppose the algorithm sees the noisy gradient $g(x)=\nabla f(x)+\xi$, where a fresh $\xi$ is drawn at every call, independent of the past, with mean zero and covariance $\sigma^2I$. Run Switching MM with step $\alpha>0$ and momentum coefficient $0\le\beta<1$ with $g$ in place of $\nabla f$ (it starts from $x_{-1}=x_0$, takes the plain step $x_{t+1}=x_t-\alpha g(x_t)$ at even $t$, and at odd $t$ takes the extrapolated step $x_{t+1}=y_t-\alpha g(y_t)$ with $y_t=x_t+\beta(x_t-x_{t-1})$). Fix an eigenvalue $\lambda$ of $A$ with normalized step $a=\alpha\lambda\in(0,1]$, let $e_t$ be the component of $x_t-x^\ast$ along its eigenvector and $\xi_t$ the component of the noise drawn at iteration $t$ (independent, mean zero, variance $\sigma^2$), and put $\rho=(1-a)\bigl(1-(1+\beta)a\bigr)$. (i) For every even $t$, $e_{t+2}=\rho\,e_t+w_t$ with $w_t=-(1-a)(1+\beta)\alpha\,\xi_t-\alpha\,\xi_{t+1}$, where the $w_t$ for even $t$ are independent, independent of $e_t$, with mean zero and $\operatorname{Var}(w_t)=\alpha^2\sigma^2\bigl[(1-a)^2(1+\beta)^2+1\bigr]$. (ii) As $t\to\infty$ through even values, $\operatorname{Var}(e_t)$ converges to $V_{\mathrm{even}}=\alpha^2\sigma^2\bigl[(1-a)^2(1+\beta)^2+1\bigr]/(1-\rho^2)$, and through odd values to $V_{\mathrm{odd}}=(1-a)^2V_{\mathrm{even}}+\alpha^2\sigma^2$. (iii) The phase averaged stationary variance $V_{\mathrm{SMM}}=\tfrac12(V_{\mathrm{even}}+V_{\mathrm{odd}})$ satisfies $V_{\mathrm{SMM}}(a)/(\alpha^2\sigma^2)=Q(a,\beta)/(1-\rho^2)$, where $Q(a,\beta)=(2+2\beta+\beta^2)-(1+\beta)(2+3\beta)\,a+(1+\beta)(1+3\beta)\,a^2-\beta(1+\beta)\,a^3$ and $1-\rho^2=a\bigl(2+\beta-(1+\beta)a\bigr)\bigl(2-(2+\beta)a+(1+\beta)a^2\bigr)$.

Source: [`book/chapters/ch07_noise.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch07_noise.tex), label `thm:smmvar`.

<a id="smm_slow_mode_noise_amplification"></a>

### Slow mode noise amplification of Switching MM

Corollary 8.5, page 129 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f(x)=\tfrac12x^\top Ax-b^\top x$ with $A$ symmetric positive definite, with the noisy gradient oracle $g(x)=\nabla f(x)+\xi$ (a fresh $\xi$ at every call, independent of the past, mean zero, covariance $\sigma^2I$), and run Switching MM with step $\alpha>0$ and momentum coefficient $0\le\beta<1$ with $g$ in place of $\nabla f$ (plain step at even $t$, extrapolated step $x_{t+1}=y_t-\alpha g(y_t)$ with $y_t=x_t+\beta(x_t-x_{t-1})$ at odd $t$). On the eigendirection of an eigenvalue $\lambda$ with normalized step $a=\alpha\lambda\in(0,1]$, let $V_{\mathrm{even}}$ and $V_{\mathrm{odd}}$ be the limits of the variance of the error component along even and odd times, $V_{\mathrm{SMM}}=\tfrac12(V_{\mathrm{even}}+V_{\mathrm{odd}})$, and $V_{\mathrm{GD}}(a)=\alpha^2\sigma^2/\bigl(a(2-a)\bigr)$ the stationary variance of noisy gradient descent $e_{t+1}=(1-a)e_t-\alpha\xi_t$ with the same step. As $a\to0$, both $aV_{\mathrm{even}}$ and $aV_{\mathrm{odd}}$ converge to $\alpha^2\sigma^2\,\dfrac{(1+\beta)^2+1}{2(2+\beta)}=\alpha^2\sigma^2\,\dfrac{1+\beta+\beta^2/2}{2+\beta}$. Consequently $\lim_{a\to0}V_{\mathrm{SMM}}(a)/V_{\mathrm{GD}}(a)=V(\beta):=\dfrac{2+2\beta+\beta^2}{2+\beta}=\beta+\dfrac{2}{2+\beta}$.

Source: [`book/chapters/ch07_noise.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch07_noise.tex), label `cor:ch7:slow`.

<a id="smm_variance_bound_uniform_in_normalized_step"></a>

### Noise bound for Switching MM uniform in the normalized step

Proposition 8.7, page 130 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f(x)=\tfrac12x^\top Ax-b^\top x$ with $A$ symmetric positive definite, with the noisy gradient oracle $g(x)=\nabla f(x)+\xi$ (a fresh $\xi$ at every call, independent of the past, mean zero, covariance $\sigma^2I$), and run Switching MM with step $\alpha>0$ and momentum coefficient $0\le\beta<1$ with $g$ in place of $\nabla f$. On the eigendirection of an eigenvalue $\lambda$ with normalized step $a=\alpha\lambda\in(0,1]$, let $V_{\mathrm{SMM}}(a)$ be the phase averaged stationary variance of the error component (the average of its limits along even and odd times), $V_{\mathrm{GD}}(a)=\alpha^2\sigma^2/\bigl(a(2-a)\bigr)$ the stationary variance of noisy gradient descent with the same step, and $V(\beta)=\beta+2/(2+\beta)$. Then for $0\le\beta<1$ and every $a\in(0,1]$, $V_{\mathrm{GD}}(a)\le V_{\mathrm{SMM}}(a)\le V(\beta)\,V_{\mathrm{GD}}(a)<\tfrac53\,V_{\mathrm{GD}}(a)$, with $V_{\mathrm{SMM}}(a)=V_{\mathrm{GD}}(a)$ only when $\beta=0$ or $a=1$, and $V_{\mathrm{SMM}}(a)=V(\beta)\,V_{\mathrm{GD}}(a)$ only when $\beta=0$.

Source: [`book/chapters/ch07_noise.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch07_noise.tex), label `prop:ch7:uniform`.

<a id="noise_floor_at_matched_speed"></a>

### Noise floor at matched speed for momentum methods and Switching MM

Proposition 8.9, page 132 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Consider noisy gradient methods on a strongly convex quadratic with additive gradient noise of variance $\sigma^2$ and step $\alpha$, in the slow mode limit $a=\alpha\lambda\to0$ along an eigendirection. A method has the effective step multiplier $m$ if its slow mode contracts like gradient descent with step $m\alpha$ (and gradient descent with step $m\alpha$ has slow mode variance $m$ times that at step $\alpha$), and the slow mode variance amplification $V$ if its stationary variance is $V$ times that of gradient descent at step $\alpha$. Its noise floor at matched speed relative to gradient descent is $V/m$. For momentum coefficient $0\le\beta<1$, this equals $1$ for every step heavy ball and for every step Nesterov momentum (where $V=m=1/(1-\beta)$), and for Switching MM (where $m=(2+\beta)/2$ and $V=(2+2\beta+\beta^2)/(2+\beta)$) it equals $\dfrac{2(2+2\beta+\beta^2)}{(2+\beta)^2}=1+\dfrac{\beta^2}{(2+\beta)^2}\ge1$, which is $1.051$ at $\beta=\log_23-1$ and increases to $10/9$ as $\beta\to1$. More generally, for any period two cyclic gradient step schedule with steps $h_1\alpha$ and $h_2\alpha$ the ratio is $\tfrac12(h_1^2+h_2^2)\big/\bigl(\tfrac12(h_1+h_2)\bigr)^2\ge1$, with equality only for equal steps.

Source: [`book/chapters/ch07_noise.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch07_noise.tex), label `prop:ch7:matched`.

<a id="gain_weight_theorem"></a>

### Gain weight theorem for the noise floor of periodic first order methods

Theorem 8.12, page 135 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Fix a period $p\ge1$ and a memory length $K\ge1$. On one eigendirection of a quadratic with normalized step $a=\alpha\lambda$, with the state $z_t=(e_t,\dots,e_{t-K+1})$ of past errors, a linear periodic first order method evaluates at step $t$, with residue $j=t\bmod p$, the gradient at the extrapolated point $y_t=\gamma_j\cdot z_t$ and updates $e_{t+1}=\delta_j\cdot z_t-h_j\bigl(a\,y_t+\alpha\xi_t\bigr)$, where $\gamma_j,\delta_j\in\mathbb{R}^K$ and the step multipliers $h_j>0$ are fixed, and the gradient noise $\xi_t$ has mean zero, is uncorrelated across steps and has variance $\sigma^2$. Consistency requires $\delta_j\cdot\iota=\gamma_j\cdot\iota=1$ for $\iota=(1,\dots,1)$. Let $S_j(a)$ be the step matrix mapping $z_t$ to $z_{t+1}$ without noise (its first row is $(\delta_j-a\,h_j\gamma_j)^\top$ and its other rows shift the entries of $z_t$ down by one place), and write the product of the $p$ step matrices of one period as $B(a)=B_0+aB_1+O(a^2)$. Assume the method is stable at $a=0$: the eigenvalue $1$ of $B_0$, with right eigenvector $\iota$, is simple and all other eigenvalues of $B_0$ lie strictly inside the unit disc. Let $w$ be the left eigenvector of $B_0$ for the eigenvalue $1$ with $w\cdot\iota=1$. The gain of the $i$th step of the period is $c_i=h_i\,w\cdot\Pi_iu_1$, where $u_1$ is the first unit vector and $\Pi_i$ is the product, at $a=0$, of the step matrices of steps $i+1,\dots,p$ of the period, and $M=\sum_{i=1}^pc_i$, which is assumed to be positive. Then, as $a\to0$: (i) the per step contraction of the slow mode is $1-(M/p)\,a+O(a^2)$; (ii) the phase averaged stationary variance of $e_t$ is $\alpha^2\sigma^2\,\dfrac{\sum_ic_i^2}{2Ma}\,(1+o(1))$; (iii) consequently the noise floor relative to gradient descent run at the same slow mode speed, with step $(M/p)\alpha$, is $\mathcal F=\dfrac{p\sum_{i=1}^pc_i^2}{\bigl(\sum_{i=1}^pc_i\bigr)^2}\ge1$, with equality if and only if all gains are equal.

Source: [`book/chapters/ch07_noise.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch07_noise.tex), label `thm:ch7:gainweight`.

<a id="rate_and_noise_for_period_two_steps"></a>

### Rate and noise Pareto frontier for period two step schedules

Proposition 8.14, page 137 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f$ be a quadratic whose Hessian has its spectrum in $[\mu,L]$, with $\kappa=L/\mu$, and consider gradient descent with cyclic steps $h_1,h_2$ (alternating absolute step lengths) of ratio $r=h_2/h_1\ge1$, under additive gradient noise. Put $G(r)=\dfrac{(1+r)^2}{4r}\ge1$. The worst case two step contraction over $[\mu,L]$ is $\max_{\lambda\in[\mu,L]}|(1-h_1\lambda)(1-h_2\lambda)|$, and the slow mode noise floor at matched speed is the slow mode stationary variance divided by that of gradient descent with the constant step of the same slow mode speed. (i) The slow mode noise floor at matched speed is $N=\dfrac{2(h_1^2+h_2^2)}{(h_1+h_2)^2}=2-\dfrac1{G(r)}$. (ii) For fixed $r$ with $G(r)\le G_{\mathrm{ch}}(\kappa)=\dfrac{2(\kappa+1)^2}{\kappa^2+6\kappa+1}$, the smallest possible worst case two step contraction over $[\mu,L]$ is $R^\ast(r)=1-\dfrac{4\kappa\,G(r)}{(\kappa+1)^2}$, attained when the steps are scaled so that $h_1+h_2=4G(r)/(\mu+L)$. For $G(r)>G_{\mathrm{ch}}(\kappa)$ the smallest worst case contraction exceeds that of the Chebyshev ratio. (iii) Consequently the Pareto frontier between rate and noise floor is $R=1-\dfrac{4\kappa}{(2-N)(\kappa+1)^2}$ for $1\le N\le N_{\mathrm{ch}}(\kappa)$, where $N_{\mathrm{ch}}=\dfrac{3\kappa^2+2\kappa+3}{2(\kappa+1)^2}$ is the noise floor of the two step Chebyshev schedule (the period two schedule minimizing the worst case two step contraction over $[\mu,L]$), whose ratio tends to $(1+\sqrt2)^2=5.83$ and whose floor tends to $3/2$ as $\kappa\to\infty$.

Source: [`book/chapters/ch07_noise.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch07_noise.tex), label `prop:ch7:pareto`.

<a id="smm_multiplicative_noise"></a>

### Switching MM under multiplicative gradient noise

Proposition 8.15, page 141 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Consider on one eigendirection of a quadratic with eigenvalue $\lambda$ and normalized step $a=\alpha\lambda$ the oracle $g(y)=(\lambda+\eta_t)\,e_y+\xi_t$ for the gradient at a point with error $e_y$ along that eigenvector, where $\eta_t,\xi_t$ are independent of the past, with mean zero, uncorrelated with each other, $\operatorname{Var}(\alpha\eta_t)=\omega^2$ and $\operatorname{Var}(\xi_t)=\sigma^2$. For Switching MM with step $\alpha$ and momentum coefficient $0\le\beta<1$ the second moments of the error obey a closed linear recursion: with $\nu_1=(1-a)^2+\omega^2$ and $\nu_2=\bigl(1-(1+\beta)a\bigr)^2+(1+\beta)^2\omega^2$, the stationary second moment at even times is $\alpha^2\sigma^2\bigl(\nu_1(1+\beta)^2+1\bigr)/(1-\nu_1\nu_2)$, finite if and only if $\nu_1\nu_2<1$. As $a\to0$: (i) under the least squares scaling $\omega^2=\vartheta a^2$ with $\vartheta$ fixed (the relative curvature noise, $\vartheta=2$ for a Gaussian regressor), the slow mode ratio of the stationary second moment of Switching MM to that of gradient descent with the same step and noise is unchanged, equal to $(2+2\beta+\beta^2)/(2+\beta)$ for every $\vartheta$; (ii) under the diffusive scaling $\omega^2=\epsilon a$, for any linear periodic first order method of period $p$ on this eigendirection (with state $z_t=(e_t,\dots,e_{t-K+1})$ and residue $j=t\bmod p$ it sets $e_{t+1}=\delta_j\cdot z_t-h_j\alpha\bigl((\lambda+\eta_t)y_t+\xi_t\bigr)$ with $y_t=\gamma_j\cdot z_t$, where $\gamma_j,\delta_j\in\mathbb{R}^K$ with $\gamma_j\cdot\iota=\delta_j\cdot\iota=1$ for $\iota=(1,\dots,1)$ and $h_j>0$ are fixed) that is stable at $a=0$ (the product $B_0$ of its step matrices over one period at $a=0$ has the simple eigenvalue $1$ and all its other eigenvalues strictly inside the unit disc) and has $M>0$, where $c_i=h_i\,w\cdot\Pi_iu_1$ is the gain of step $i$ ($w$ the left eigenvector of $B_0$ for the eigenvalue $1$ with $w\cdot\iota=1$, $u_1$ the first unit vector, $\Pi_i$ the product at $a=0$ of the step matrices of steps $i+1,\dots,p$) and $M=\sum_ic_i$, the stationary second moment of the error satisfies $a\,\mathbb{E}[e^2]\to\alpha^2\sigma^2\sum_ic_i^2\big/\bigl(2M-\epsilon\sum_ic_i^2\bigr)$, so that the method is mean square stable for small $a$ exactly when $\epsilon<\epsilon^\ast=2M/\sum_ic_i^2$.

Source: [`book/chapters/ch07_noise.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch07_noise.tex), label `prop:ch7:mult`.

## Chapter 9. Nonlinear Difference Equations and Momentum MM

<a id="certified_growth_exponents_for_the_inverse_collatz_tree"></a>

### Certified growth exponents for the inverse Collatz tree from a simplified difference inequality method

Theorem 9.3, page 158 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $T$ be the accelerated Collatz map, $T(n)=n/2$ for even $n$ and $T(n)=(3n+1)/2$ for odd $n$. For an integer $n\not\equiv0\pmod 3$ and $X\ge1$ let $\pi_n(X)=\#\{m\le Xn:\ \text{the orbit of }m\text{ under }T\text{ passes through }n\}$. For an integer $j\ge1$ and $n\not\equiv0\pmod3$ with $2^jn\equiv2\pmod3$ put $c_j(n)=(2^{j+1}n-1)/3$, the odd preimage of $2^jn$, and call $j$ fertile for $n$ if $2^jn\equiv2$ or $8\pmod9$. Fix $k\ge1$ and let $\mathcal R_k$ be the set of residues modulo $3^k$ prime to $3$. For $g>0$ and a positive vector $(w_r)_{r\in\mathcal R_k}$ define $(\Phi_gw)_r=\min_{\ell}\sum_{j\ \mathrm{fertile\ for}\ \ell}\bigl(3/2^{j+1}\bigr)^g\,w_{[c_j(\ell)]}$, where the minimum is over the three lifts $\ell$ of $r$ to residues modulo $3^{k+1}$, the sum is over all $j\ge1$, and fertility of $j$ and the class $[c_j(\ell)]\in\mathcal R_k$ (the residue of $c_j(n)$ modulo $3^k$ for any integer $n\equiv\ell$) depend only on $\ell$ and $j$. If $\Phi_gw\ge(1+\eta)w$ componentwise for some $\eta>0$, then there is $C>0$ such that $\pi_n(X)\ge C\,w_{[n]}\,X^g$ for all $X\ge1$ and all $n\not\equiv0\pmod3$ not lying on a cycle of $T$, where $[n]$ is the residue of $n$ modulo $3^k$.

Source: [`book/chapters/ch05b_difference.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch05b_difference.tex), label `thm:ch5b:exponent`.

<a id="simplified_counting_bound_for_integers_reaching_one"></a>

### Counting bound with exponent 0.477 for integers reaching 1

Corollary 9.4, page 159 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. For all sufficiently large $x$, at least $x^{0.477}$ integers $m\le x$ have orbits under the accelerated Collatz map $T$ ($T(n)=n/2$ for even $n$ and $T(n)=(3n+1)/2$ for odd $n$) that reach $1$.

Source: [`book/chapters/ch05b_difference.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch05b_difference.tex), label `cor:ch5b:count`.

<a id="smm_monotone_on_alternate_iterates"></a>

### Switching MM decreases the objective on alternate iterates for every smooth function

Proposition 9.5, page 161 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f:\mathbb{R}^d\to\mathbb{R}$ be $L$ smooth (its gradient is Lipschitz with constant $L$) and bounded below by $f_{\inf}$, and run Switching MM with step $1/L$ and momentum coefficient $\beta\in[0,1)$ (it starts from $x_{-1}=x_0$, takes the plain gradient step $x_{t+1}=x_t-\nabla f(x_t)/L$ at even $t$, and at odd $t$ takes the extrapolated step $x_{t+1}=y_t-\nabla f(y_t)/L$ with $y_t=x_t+\beta(x_t-x_{t-1})$). For even $t$, $f(x_{t+2})\le f(x_t)-\dfrac{1-\beta^2}{2L}\|\nabla f(x_t)\|^2-\dfrac1{2L}\|\nabla f(y_{t+1})\|^2$, where $y_{t+1}=x_{t+1}+\beta(x_{t+1}-x_t)$. Consequently $\min_{t<2K,\ t\ \mathrm{even}}\|\nabla f(x_t)\|^2\le\dfrac{2L\bigl(f(x_0)-f_{\inf}\bigr)}{(1-\beta^2)K}$.

Source: [`book/chapters/ch05b_difference.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch05b_difference.tex), label `prop:ch5b:monotone`.

<a id="smm_exact_joint_spectral_radius"></a>

### Exact joint spectral radius of the descent and injection matrices of Switching MM

Proposition 9.6, page 162 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. For parameters $a$ and $\beta$ (in Switching MM on an eigendirection of a quadratic, $a=\alpha\lambda>0$ and $0\le\beta<1$), let $D=\begin{pmatrix}1-a&0\\1&0\end{pmatrix}$ and $N=\begin{pmatrix}(1-a)(1+\beta)&-(1-a)\beta\\1&0\end{pmatrix}$ be the matrices that map the state $(e_t,e_{t-1})$ to $(e_{t+1},e_t)$ for a descent step and for an injection step, and let $\hat\rho(\{D,N\})=\lim_{k\to\infty}\max_{M_i\in\{D,N\}}\|M_k\cdots M_1\|^{1/k}$ be their joint spectral radius. Write $D=vu^\top$ with $u=(1,0)^\top$ and $v=(1-a,1)^\top$, and put $\tau_j=u^\top N^jv$ for $j\ge0$. Then $\hat\rho(\{D,N\})=\max\Bigl(\rho(N),\ \sup_{j\ge0}|\tau_j|^{1/(j+1)}\Bigr)$, where $\rho(N)$ is the spectral radius of $N$.

Source: [`book/chapters/ch05b_difference.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch05b_difference.tex), label `prop:ch5b:jsr`.

<a id="smm_transient_bound_conjecture"></a>

### Transient bound conjecture for the injection matrix of Switching MM

Conjecture 9.7, page 163 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $N=\begin{pmatrix}(1-a)(1+\beta)&-(1-a)\beta\\1&0\end{pmatrix}$, $u=(1,0)^\top$ and $v=(1-a,1)^\top$, and let $\tau_j=u^\top N^jv$ for $j\ge0$ (the first coordinate after $j$ injection steps of Switching MM, with normalized step $a$ and momentum coefficient $\beta$, started from the state left by one descent step). Conjecture: for all $a\in(0,1]$, $\beta\in[0,1)$ and $j\ge0$, $|\tau_j|\le1-a$.

Source: [`book/chapters/ch05b_difference.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch05b_difference.tex), label `conj:ch5b:transient`.

<a id="every_switching_rule_is_stable_on_quadratics"></a>

### Every switching rule between descent and injection is stable on quadratics, conditional on the transient bound

Corollary 9.8, page 163 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: conditional on the transient bound $|\tau_j|<1$ for every $j\ge0$ at the fixed normalized step $a$, where $\tau_j=u^\top N^jv$ is defined in the statement; it is the content of the transient bound conjecture ($|\tau_j|\le1-a$ for all $a\in(0,1]$, $\beta\in[0,1)$ and $j\ge0$), which has been checked numerically on a finite grid of parameters and is not proved.

Statement. Let $0\le\beta<1$ and fix $a\in(0,1)$. On an eigendirection of a quadratic with eigenvalue $\lambda$ and normalized step $a=\alpha\lambda$, let $D(a)=\begin{pmatrix}1-a&0\\1&0\end{pmatrix}$ and $N(a)=\begin{pmatrix}(1-a)(1+\beta)&-(1-a)\beta\\1&0\end{pmatrix}$ be the matrices that map the state $(e_t,e_{t-1})$ to $(e_{t+1},e_t)$ for a descent step and for an injection step of Switching MM, let $\hat\rho(\{D(a),N(a)\})$ be their joint spectral radius, and let $\tau_j=u^\top N(a)^jv$ with $u=(1,0)^\top$ and $v=(1-a,1)^\top$. Suppose $|\tau_j|<1$ for all $j\ge0$. Then on a quadratic with spectrum in $(0,L]$ and step $\alpha=1/L$, every switching sequence between descent and injection steps, whether periodic, random, or chosen adaptively from the state, converges to the minimizer on the eigendirection with $a=\alpha\lambda$, with exponential rate at most $\hat\rho(\{D(a),N(a)\})<1$.

Source: [`book/chapters/ch05b_difference.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch05b_difference.tex), label `cor:ch5b:anyswitch`.

<a id="smm_linear_rate_on_smooth_strongly_convex_functions"></a>

### Exact linear rate of Switching MM on smooth strongly convex functions

Theorem 9.9, page 164 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f$ be $\mu$ strongly convex and $L$ smooth with minimizer $x^\ast$, put $\kappa=L/\mu$, and run Switching MM with step $\alpha=1/L$ and momentum coefficient $\beta\in[0,1)$ (it starts from $x_{-1}=x_0$, takes the plain gradient step $x_{t+1}=x_t-\alpha\nabla f(x_t)$ at even $t$, and at odd $t$ takes the extrapolated step $x_{t+1}=y_t-\alpha\nabla f(y_t)$ with $y_t=x_t+\beta(x_t-x_{t-1})$). For every even $t$, $\|x_{t+2}-x^\ast\|\le q\,\|x_t-x^\ast\|$ with $q=q(\beta,\kappa)=\Bigl(1-\dfrac1\kappa\Bigr)\max\Bigl(1-\dfrac{1+\beta}{\kappa},\ \beta\Bigr)$, and $f(x_{t+2})-f^\ast\le q^2\,\bigl(f(x_t)-f^\ast\bigr)$.

Source: [`book/chapters/ch05b_difference.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch05b_difference.tex), label `thm:ch5b:rate`.

<a id="smm_speed_relative_to_gradient_descent"></a>

### Speed of Switching MM relative to gradient descent on smooth strongly convex functions

Corollary 9.10, page 165 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f$ be $\mu$ strongly convex and $L$ smooth with condition number $\kappa=L/\mu$, run Switching MM with step $1/L$ and momentum coefficient $\beta\in[0,1)$, and let $q=(1-1/\kappa)\max\bigl(1-(1+\beta)/\kappa,\ \beta\bigr)$ be the worst case, over this class of functions, of the factor by which two steps contract the distance to the minimizer. For $\kappa\ge(1+\beta)/(1-\beta)$ the per step rate of Switching MM is $q^{1/2}=\bigl((1-1/\kappa)(1-(1+\beta)/\kappa)\bigr)^{1/2}$, so that $-\kappa\log q^{1/2}\to(2+\beta)/2$ as $\kappa\to\infty$, against $1$ for gradient descent with the same step.

Source: [`book/chapters/ch05b_difference.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch05b_difference.tex), label `cor:ch5b:speed`.

<a id="modified_equation_for_period_two_steps"></a>

### Modified equation for period two gradient steps

Proposition 9.11, page 165 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f:\mathbb{R}^d\to\mathbb{R}$ be three times differentiable, write $G_h(x)=x-h\nabla f(x)$, and consider the composition of gradient steps of lengths $h_1$ and $h_2$, in either order. To second order in $h=h_1+h_2$ (that is, up to terms of order $h^3$) it agrees with the time $h$ flow of the modified equation $\dot x=-\nabla f(x)-\dfrac{h_1^2+h_2^2}{2(h_1+h_2)}\,\nabla^2f(x)\nabla f(x)=-\nabla\Bigl(f+\dfrac{\bar h}{4}\|\nabla f\|^2\Bigr)$, where $\bar h=\dfrac{h_1^2+h_2^2}{h_1+h_2}$.

Source: [`book/chapters/ch05b_difference.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch05b_difference.tex), label `prop:ch5b:modified`.

## Chapter 10. Three Point Inequalities and Similar Triangles

<a id="switching_mm_convex_rate_by_three_point_inequality"></a>

### Convex rate for Switching MM when the long step is at most the golden ratio

Theorem 10.5, page 173 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f:\mathbb{R}^d\to\mathbb{R}$ be convex and differentiable with an $L$ Lipschitz gradient ($f$ is $L$ smooth), let $x_*$ be a minimizer of $f$ and $f_*=f(x_*)$, let $\varphi=(1+\sqrt5)/2$ be the golden ratio, and let $0\le\beta\le\varphi-1$. Let $x_0,x_1,x_2,\ldots$ be the iterates of Switching MM with step size $1/L$ and momentum $\beta$ started at $x_0$: put $x_{-1}=x_0$, and for $t\ge0$ let $x_{t+1}=x_t-\frac{1}{L}\nabla f(x_t)$ if $t$ is even, and $x_{t+1}=y_t-\frac{1}{L}\nabla f(y_t)$ with $y_t=x_t+\beta(x_t-x_{t-1})$ if $t$ is odd. Let $R=\|x_0-x_*\|$ (Euclidean norm). Then for every $j\ge0$, $(2+\beta)\bigl(f(x_{2j+2})-f_*\bigr)\le\frac{L}{2}\bigl(\|x_{2j}-x_*\|^2-\|x_{2j+2}-x_*\|^2\bigr)$, and consequently, for every even $N$, $f(x_N)-f_*\le\frac{LR^2}{(2+\beta)N}$.

Source: [`book/chapters/ch07b_threepoint.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch07b_threepoint.tex), label `thm:ch7b:smmrate`.

<a id="switching_mm_exact_worst_case_conjecture"></a>

### Exact convex worst case of Switching MM (conjecture)

Conjecture 10.6, page 174 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $\beta_2=0.68614\ldots$ be the positive real root of $2\beta^3+5\beta^2-3=0$. Let $f:\mathbb{R}^d\to\mathbb{R}$ be convex and differentiable with an $L$ Lipschitz gradient ($f$ is $L$ smooth) and a minimizer $x_*$ (with $f_*=f(x_*)$), and let $x_0,x_1,x_2,\ldots$ be the iterates of Switching MM with step size $1/L$ and momentum $\beta$ started at $x_0$: put $x_{-1}=x_0$, and for $t\ge0$ let $x_{t+1}=x_t-\frac{1}{L}\nabla f(x_t)$ if $t$ is even, and $x_{t+1}=y_t-\frac{1}{L}\nabla f(y_t)$ with $y_t=x_t+\beta(x_t-x_{t-1})$ if $t$ is odd. Conjecture: for $0\le\beta\le\beta_2$ and every even $N$, the worst case of Switching MM over convex $L$ smooth functions with $\|x_0-x_*\|\le R$ (the supremum over all dimensions $d$, all such $f$ and all such $x_0$) is $\sup\bigl(f(x_N)-f_*\bigr)=\frac{LR^2}{(4+2\beta)N+2}$.

Source: [`book/chapters/ch07b_threepoint.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch07b_threepoint.tex), label `conj:ch7b:tight`.

<a id="acceleration_survives_every_switching_rule"></a>

### Acceleration by similar triangles survives every switching rule

Proposition 10.7, page 176 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f:\mathbb{R}^d\to\mathbb{R}$ be convex and differentiable with an $L$ Lipschitz gradient ($f$ is $L$ smooth), let $x_*$ be a minimizer of $f$ and $f_*=f(x_*)$. A similar triangles step acts on a state $(x,z,A)$ with $A\ge0$: it chooses $a>0$ with $La^2=A+a$, sets $A^+=A+a$, $y=(Ax+az)/A^+$, $z^+=z-a\nabla f(y)$ and $x^+=(Ax+az^+)/A^+$ (which equals $y-\frac{1}{L}\nabla f(y)$), and takes $(x^+,z^+,A^+)$ as the new state. A plain gradient step is $x^+=x-\frac{s}{L}\nabla f(x)$ with $0<s<2$, and leaves $z$ and $A$ unchanged. Start from $x=z=x_0$ and $A=0$ and apply any sequence of iterations, each of which is a similar triangles step or a plain gradient step; the order of the two kinds of step may be arbitrary and may depend on the iterates. Then the energy $E=A\bigl(f(x)-f_*\bigr)+\frac{1}{2}\|z-x_*\|^2$ never increases, and after $J$ similar triangles steps $f(x)-f_*\le\frac{2L\|x_0-x_*\|^2}{J^2}$. In particular the switched method that alternates the two kinds of step, beginning with a similar triangles step, satisfies $f(x_N)-f_*\le\frac{8L\|x_0-x_*\|^2}{N^2}$ after $N$ iterations.

Source: [`book/chapters/ch07b_threepoint.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch07b_threepoint.tex), label `prop:ch7b:switch`.

<a id="bregman_decomposition_of_a_syracuse_orbit"></a>

### Bregman decomposition of the logarithmic change along a Syracuse orbit

Proposition 10.8, page 177 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $S(y)=(3y+1)/2^{v_2(3y+1)}$ be the Syracuse map on odd positive integers, where $v_2$ is the $2$ adic valuation, and let $D(t)=t-1-\log t\ge0$ for $t>0$ (the Itakura and Saito divergence of $t$ from $1$). For an odd positive integer $y_0$ let $y_{i+1}=S(y_i)$, $k_i=v_2(3y_i+1)$ and $M_i=y_{i+1}/y_i=(3+1/y_i)/2^{k_i}$. Then for every odd $y_0$ and every $N$, $\log y_N-\log y_0=\sum_{i<N}(M_i-1)-\sum_{i<N}D(M_i)$. In the stochastic model in which the $M_i$ are independent copies of $M=3/2^k$ with $k$ geometric, $\mathbb{P}(k=j)=2^{-j}$ for $j\ge1$, the partial sums $\sum_{i<N}(M_i-1)$ form a martingale in $N$, $\mathbb{E}M=1$, and the expected dissipation per step equals the drift exactly, $\mathbb{E}D(M)=-\mathbb{E}\log M=\log(4/3)$.

Source: [`book/chapters/ch07b_threepoint.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch07b_threepoint.tex), label `prop:ch7b:bregman`.

<a id="uniform_dissipation_of_the_syracuse_map"></a>

### Uniform lower bound on the Bregman dissipation of the Syracuse map

Proposition 10.9, page 179 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $S(y)=(3y+1)/2^{v_2(3y+1)}$ be the Syracuse map on odd positive integers, where $v_2$ is the $2$ adic valuation, and let $D(t)=t-1-\log t$ for $t>0$. Then for every odd $y\ge3$, $D\bigl(S(y)/y\bigr)\ge c_*=D\bigl(\tfrac{7}{9}\bigr)=\log\tfrac{9}{7}-\tfrac{2}{9}=0.029092\ldots$, with equality only at $y=9$, and $\liminf_{y\to\infty}D\bigl(S(y)/y\bigr)=D\bigl(\tfrac{3}{4}\bigr)=\log\tfrac{4}{3}-\tfrac{1}{4}=0.037682\ldots$, where $y$ runs over odd integers.

Source: [`book/chapters/ch07b_threepoint.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch07b_threepoint.tex), label `prop:ch7b:dissipation`.

## Chapter 11. Cycles and Continued Fractions

<a id="diophantine_insufficiency_cycle_lengths"></a>

### Diophantine insufficiency, infinitely many cycle lengths survive any verification range

Proposition 11.20, page 205 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $\theta=\log_2 3$, let $\log$ denote the natural logarithm, and for an integer $r\geq1$ let $\delta(r)=\lceil r\theta\rceil-r\theta$ be the distance from $r\theta$ up to the next integer. Say that $r$ is excluded at verification range $X$ if $\delta(r)\geq r/(3X\log2)$. (If every element of every nontrivial cycle of the Syracuse map on the positive integers is at least $X$, then such a cycle with $r$ odd elements can exist only for an $r$ that is not excluded at $X$. The Syracuse map sends an odd $y$ to $(3y+1)/2^{v_2(3y+1)}$, where $v_2$ is the 2 adic valuation.) For every real $X\geq1$ there are infinitely many integers $r\geq1$ with $0<\delta(r)<\dfrac{r}{3X\log2}$, that is, infinitely many $r$ that are not excluded at verification range $X$.

Source: [`book/chapters/ch08_cycles.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch08_cycles.tex), label `prop:ch8:insufficient`.

<a id="window_lemma_least_cycle_element"></a>

### Window lemma, the least element of a cycle lies in an explicit window

Lemma 11.21, page 206 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $T$ be the accelerated Collatz map on the positive integers, $T(x)=x/2$ for even $x$ and $T(x)=(3x+1)/2$ for odd $x$. Let $x_1,\dots,x_r$ be the odd elements of a cycle of $T$ with $m=r$ odd steps and $k$ steps of $T$ in all (one turn of the cycle), and put $\Lambda=k\log2-m\log3>0$, where $\log$ is the natural logarithm. Then every element $x$ of the cycle, odd or even, satisfies $\frac{3^m-2^m}{2^k-3^m}\le x$, and $M=\min_i x_i$ satisfies $M<\frac{m}{3\Lambda}\le F(m)$, where $F(m)=\frac{m}{3\,\delta(m)\log2}$, $\delta(m)=\lceil m\theta\rceil-m\theta$ and $\theta=\log_2 3$. Here $\frac{3^m-2^m}{2^k-3^m}=\frac{1-(2/3)^m}{e^{\Lambda}-1}$, and $\frac{m}{3\Lambda}=F(m)$ exactly when $k=\lceil m\theta\rceil$. In particular the least element lies in the window $\Bigl[\frac{1-(2/3)^m}{e^{\Lambda}-1},\ \frac{m}{3\Lambda}\Bigr)$, whose endpoints have ratio $\frac m3\cdot\frac{e^{\Lambda}-1}{\Lambda}\bigl(1-(2/3)^m\bigr)^{-1}$, which is $\frac m3(1+o(1))$ as $\Lambda\to0$ and $m\to\infty$.

Source: [`book/chapters/ch08_cycles.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch08_cycles.tex), label `lem:ch8:window`.

<a id="what_arithmetic_would_suffice"></a>

### What arithmetic would suffice to exclude all Collatz cycles

Corollary 11.22, page 207 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: conditional on the hypothesis that, for the function $f$ of the statement, every nontrivial cycle of the accelerated map with $m$ odd steps and $k$ steps in all has least element $M<m/(3\Lambda f(m))$, where $\Lambda=k\log2-m\log3$; no such bound is proved in the book.

Statement. Let $T$ be the accelerated Collatz map on the positive integers, $T(x)=x/2$ for even $x$ and $T(x)=(3x+1)/2$ for odd $x$, and call a cycle of $T$ nontrivial if it is not $\{1,2\}$. Let $\log$ be the natural logarithm, $\theta=\log_2 3$, $\delta(m)=\lceil m\theta\rceil-m\theta$ and $F(m)=\frac{m}{3\,\delta(m)\log2}$. For a cycle with $m$ odd steps and $k$ steps in all, let $M$ be its least element and $\Lambda=k\log2-m\log3$. Suppose that for some function $f$ every nontrivial cycle with $m$ odd steps and $k$ steps in all satisfied the sharpened upper bound $M<\frac{m}{3\Lambda f(m)}$. (i) If $f(m)\ge\frac m3\cdot\frac{e^{\Lambda}-1}{\Lambda}\bigl(1-(2/3)^m\bigr)^{-1}$ for every pair $(m,k)$ with $m\theta<k<2m$, then there are no nontrivial cycles. (ii) If $f(m)\ge\frac m3(1+\eta)$ for some $\eta>0$ and all $m\ge m_0$, where $(1+\eta)\bigl(1-(2/3)^{m_0}\bigr)>1$, let $\Lambda_0>0$ solve $\frac{e^{\Lambda_0}-1}{\Lambda_0}=(1+\eta)\bigl(1-(2/3)^{m_0}\bigr)$. Then every nontrivial cycle with $m\ge m_0$ odd steps has least element $M<C=\frac{1}{(1+\eta)\Lambda_0}$, which is about $\frac1{2\eta}$ when $\eta$ is small and $m_0$ is large. Consequently, if every element of every nontrivial cycle of positive integers is at least $X$ (a verified range, which holds for instance when every positive integer below $X$ reaches $1$) and $X\ge\max\bigl(C,\max_{m<m_0}F(m)\bigr)$, then there are no nontrivial cycles.

Source: [`book/chapters/ch08_cycles.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch08_cycles.tex), label `cor:ch8:window`.

## Chapter 12. The Divergence Half

<a id="cesaro_identity"></a>

### Cesàro identity for the halving exponents of a Syracuse orbit

Proposition 12.3, page 223 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $S$ be the Syracuse map, which sends an odd positive integer $y$ to $(3y+1)/2^{v_2(3y+1)}$, where $v_2$ is the 2 adic valuation. For an odd positive integer $x_0$ let $y_j=S^{(j)}(x_0)$ be the iterates ($y_0=x_0$), let $k_j=v_2(3y_{j-1}+1)\geq1$ be the halving exponents, so that $y_j=(3y_{j-1}+1)/2^{k_j}$, and let $\bar k_N=\frac1N\sum_{j=1}^Nk_j$. Let $\log$ be the natural logarithm. Then for every odd positive integer $x_0$ and every $N\geq1$, $\log y_N-\log x_0=\sum_{j=1}^{N}\log\Bigl(3+\frac{1}{y_{j-1}}\Bigr)-N\bar k_N\log2$. Equivalently, $\bar k_N=\log_2 3+\frac{1}{N}\Bigl(\log_2x_0-\log_2y_N+\sum_{j=1}^{N}\log_2\Bigl(1+\frac{1}{3y_{j-1}}\Bigr)\Bigr)$. Consequently, if the orbit of $x_0>1$ reaches $1$ at stride $N$ (that is, $y_N=1$), then $\bar k_N>\log_2 3+\frac{\log_2x_0}{N}>\log_2 3$.

Source: [`book/chapters/ch09_divergence.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch09_divergence.tex), label `prop:cesaro`.

<a id="cesaro_liminf_conjecture"></a>

### Cesàro liminf conjecture on the running mean of the halving exponents

Conjecture 12.6, page 224 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $S$ be the Syracuse map, which sends an odd positive integer $y$ to $(3y+1)/2^{v_2(3y+1)}$, where $v_2$ is the 2 adic valuation. For an odd positive integer $x_0$ let $y_j=S^{(j)}(x_0)$ be the iterates ($y_0=x_0$) and $k_j=v_2(3y_{j-1}+1)$ the halving exponents. Conjecture: for every odd positive integer $x_0$, $\liminf_{N\to\infty}\frac1N\sum_{j=1}^Nk_j>\log_2 3$.

Source: [`book/chapters/ch09_divergence.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch09_divergence.tex), label `conj:cesaro`.

## Chapter 13. Positive Density from the Backward Tree

<a id="bounded_tao_would_imply_collatz_given_positive_density"></a>

### A bounded form of Tao's theorem would imply the Collatz conjecture, given the positive density theorems

Proposition 13.3, page 241 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: conditional on the positive density theorem of Mazur (2026) or of Lau (2026), recent unrefereed preprints with Lean formalizations whose proofs the book does not check: for every integer $t\ge1$ not divisible by $3$, the set of positive integers whose orbit under the Collatz map reaches $t$ has positive lower natural density; the computer verification of Bařina (2025) that every positive integer up to $2^{71}$ has an orbit that reaches $1$, which the book cites and does not repeat.

Statement. Let $T_0$ be the Collatz map on the positive integers, $T_0(n)=n/2$ for even $n$ and $T_0(n)=3n+1$ for odd $n$. The orbit of $n$ reaches $1$ if $T_0^{(k)}(n)=1$ for some $k\ge0$, and $\mathrm{Col}_{\min}(n)=\min_{k\ge0}T_0^{(k)}(n)$ is the least value on the orbit of $n$. For a set $A$ of positive integers, the lower and upper natural densities are the $\liminf$ and $\limsup$ as $x\to\infty$ of $x^{-1}\#\{n\le x:n\in A\}$, the lower and upper logarithmic densities are the $\liminf$ and $\limsup$ of $(\log x)^{-1}\sum_{n\le x,\,n\in A}1/n$, and $A$ has logarithmic density one if both of the latter equal $1$. Assume the positive density theorem of Mazur or of Lau: for every integer $t\ge1$ not divisible by $3$, the set of positive integers whose orbit under $T_0$ reaches $t$ has positive lower natural density. Then (a) if the set of $n$ whose orbit reaches $1$ has logarithmic density one, every orbit reaches $1$: the positive density theorem gives this when the set has upper natural density one (the corollary of Mazur and of Lau), and for every set the lower natural density is at most the lower logarithmic density, which is at most the upper logarithmic density, which is at most the upper natural density, so logarithmic density one implies upper natural density one. And (b) since every positive integer up to $2^{71}$ reaches $1$ (verified by Bařina), if $\mathrm{Col}_{\min}(n)\le2^{71}$ for all $n$ in a set of logarithmic density one, then the set of $n$ whose orbit reaches $1$ has logarithmic density one, so every orbit reaches $1$ and the Collatz conjecture holds. Tao's theorem proves the hypothesis of (b) with $2^{71}$ replaced by any function $f(n)$ that tends to infinity.

Source: [`book/chapters/ch09b_backward.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch09b_backward.tex), label `prop:ch9b:bounded`.

## Chapter 15. A Computational Path

<a id="residue_certificates_cannot_prove_descent_for_all_integers"></a>

### Residue certificates cannot prove descent for all integers (the exceptional 2 adic set)

Proposition 15.1, page 267 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $T$ be the accelerated Collatz map on the $2$ adic integers $\mathbb{Z}_2$, given by $T(x)=x/2$ for even $x$ and $T(x)=(3x+1)/2$ for odd $x$, and for $x\in\mathbb{Z}_2$ and $j\ge1$ let $m_j(x)$ be the number of odd values among $x,T(x),\ldots,T^{(j-1)}(x)$. A residue class modulo $2^k$ is uncovered at level $k$ if its members have $3^{m_j}>2^j$ for every $1\le j\le k$, that is, if their coefficient stopping time $\tau=\min\{j\ge1:3^{m_j}<2^j\}$ exceeds $k$; a residue certificate of level $k$ covers exactly the other classes (every sufficiently large member of such a class falls below itself within $k$ steps). Let $E_k\subset\mathbb{Z}_2$ be the union of the classes uncovered at level $k$, and let $E=\{x\in\mathbb{Z}_2:3^{m_j(x)}>2^j\text{ for every }j\ge1\}$ be the set of $2$ adic integers with infinite coefficient stopping time. Then: (i) $E$ is a closed subset of $\mathbb{Z}_2$ of Haar measure zero, $E=\bigcap_{k\ge1}E_k$, and $E$ contains $-1$, $-5$ and $-17$, one point of each of the three known cycles of $T$ on the negative integers. (ii) Every uncovered class at every level contains infinitely many positive integers. Hence no residue certificate of finite level covers all integers $n\ge2$. (iii) If a positive integer $n\ge2$ lies in $E$, its orbit under $T$ never falls below $n$. Hence if every integer $n\ge2$ has a finite stopping time (the descent form of the Collatz conjecture, where the stopping time of $n$ is $\min\{k\ge1:T^{(k)}(n)<n\}$), then no integer $n\ge2$ lies in $E$; conversely, if no integer $n\ge2$ lies in $E$ and the coefficient stopping time conjecture of Terras holds (the stopping time of every $n\ge2$ equals its coefficient stopping time), then every $n\ge2$ has a finite stopping time.

Source: [`book/chapters/ch10b_computation.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch10b_computation.tex), label `prop:ch10b:residue`.

## Chapter 16. Samplers and the Horizon

<a id="two_energies_coefficient_and_total_stopping_time"></a>

### Two energies under Haar measure, the coefficient stopping time and the total stopping time

Proposition 16.1, page 275 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $\mu$ be the Haar probability measure on the $2$ adic integers $\mathbb{Z}_2$ and $T$ the accelerated Collatz map, $T(x)=x/2$ for even $x$ and $T(x)=(3x+1)/2$ for odd $x$. For $x\in\mathbb{Z}_2$ let $m_k(x)$ be the number of odd values among $x,T(x),\ldots,T^{(k-1)}(x)$, let $\tau(x)=\min\{k\ge1:3^{m_k(x)}<2^k\}$ be the coefficient stopping time, and let $\sigma_\infty(x)=\min\{k\ge0:T^{(k)}(x)=1\}$ be the total stopping time counted in steps of $T$ (both equal $\infty$ if the set is empty). Let $F(m)=\mu\{\tau\le m\}$, call a binary word alive if every prefix of length $i\ge1$ with $l$ ones satisfies $3^l>2^i$, and let $N_m(j)$ be the number of alive words of length $m$ with $j$ ones. Then: (i) For every $m$ the level set $\{\tau>m\}\subset\mathbb{Z}_2$ is a union of residue classes modulo $2^m$, hence open and closed, and its Haar measure is $Z_\tau(m)=1-F(m)=2^{-m}\sum_jN_m(j)$, which is computable exactly in $O(m^2)$ arithmetic operations by the recursion $N_k(j)=[3^j>2^k]\,(N_{k-1}(j)+N_{k-1}(j-1))$ with $N_0(0)=1$ and $N_{k-1}(-1)=0$. The sets $\{\tau>m\}$ decrease in $m$ and their intersection $\{\tau=\infty\}$ has measure zero. (ii) For every $m$ the set $\{x\in\mathbb{Z}_2:\sigma_\infty(x)\le m\}$ is finite, and it contains at most one positive integer in each residue class modulo $2^m$. Consequently $\{\sigma_\infty>m\}$ has Haar measure one for every $m$, and so does the intersection $\{\sigma_\infty=\infty\}$ of these sets: the energy $\sigma_\infty$ is almost surely infinite under $\mu$, and any simulation with this energy must use a uniform prior on a finite range, which the book takes to be the integers in $[2,2^B)$, with $Z_\sigma(m)=\#\{2\le n<2^B:\sigma_\infty(n)>m\}/(2^B-2)$.

Source: [`book/chapters/ch10c_samplers.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch10c_samplers.tex), label `prop:ch10b:horizons`.

<a id="perfect_split_sampling_coefficient_stopping_time"></a>

### Perfect split sampling for the coefficient stopping time

Theorem 16.2, page 276 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $\mu$ be the Haar probability measure on $\mathbb{Z}_2$, let $T$ be the accelerated Collatz map ($T(x)=x/2$ for even $x$ and $(3x+1)/2$ for odd $x$), let $m_k(x)$ be the number of odd values among $x,T(x),\ldots,T^{(k-1)}(x)$, and let $\tau(x)=\min\{k\ge1:3^{m_k(x)}<2^k\}$ be the coefficient stopping time ($\tau=\infty$ if there is no such $k$), a function of the residue class of $x$. Put $Z_0=1$ and $Z_k=\mu\{\tau>k\}$, and let $R_k$ be the number of residue classes modulo $2^k$ with $\tau=k$ (so $R_k>0$ exactly when $Z_k<Z_{k-1}$). Fix a depth $K\ge1$, let the levels be the integers $0\le k\le K$, put $\tau_K=\min(\tau,K+1)$, and let the discrete weights be $\omega_0=1$ and $\omega_k=1/Z_k-1/Z_{k-1}$ for $1\le k\le K$, so that $W(l)=\sum_{k<l}\omega_k=1/Z_{l-1}$ for $1\le l\le K+1$. Then the augmented law $\pi_\omega(x,k)\propto[\tau_K(x)>k]\,\mu(dx)\,\omega_k$ on $\mathbb{Z}_2\times\{0,\ldots,K\}$ has the finite mass $1+\sum_{k=1}^K\lambda_k$, where $\lambda_k=Z_k\omega_k=1-Z_k/Z_{k-1}$ is the hazard, and its Gibbs sampler can be run exactly, with no Markov chain inside the level sets: (i) given $x$, the level $k$ is drawn from the weights $\omega_k$ restricted to $\{0,\ldots,\tau_K(x)-1\}$, a finite discrete law; (ii) given $k$, a point with law $\mu$ restricted to $\{\tau>k\}$ is drawn exactly, as far as $\tau_K$ can see it, which is to say its residue modulo $2^K$. For this, call a binary word alive if every prefix of length $i\ge1$ with $l$ ones satisfies $3^l>2^i$, and let $N_i(j)$ be the number of alive words of length $i$ with $j$ ones. An alive word of length $k$ is drawn uniformly by backward sampling: the terminal count $j$ is chosen with probability $N_k(j)/\sum_{j'}N_k(j')$, and the parities are generated from the last to the first, with probability $N_{i-1}(j)/N_i(j)$ of a zero and $N_{i-1}(j-1)/N_i(j)$ of a one (reducing $j$ by one after a one). The word is followed by $K-k$ independent fair parities, and the word of length $K$ is converted to a residue modulo $2^K$ by the Terras bijection between residues modulo $2^K$ and parity words of length $K$. The cost is $O(K)$ per draw after an $O(K^2)$ table. The level marginal of the chain is proportional to $1$ at $k=0$ and to $\lambda_k$ for $1\le k\le K$, and the $x$ marginal is proportional to $\mu$ weighted by $1/Z_{\tau_K(x)-1}$, so that the event $\{\tau>K\}$ has stationary probability $1/(1+\sum_{k\le K}\lambda_k)$. Every level $k\le K$ with $R_k>0$ is visited with the positive frequency $\lambda_k/(1+\sum_{k'\le K}\lambda_{k'})$. Without the truncation the augmented law has infinite mass, since $\prod_{k\le K}(1-\lambda_k)=Z_K\to0$ forces $\sum_k\lambda_k=\infty$, and every level is then visited with long run frequency zero.

Source: [`book/chapters/ch10c_samplers.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch10c_samplers.tex), label `thm:ch10b:perfect`.

<a id="orbit_coalescence_shifts_total_stopping_time"></a>

### Coalescence of orbits shifts the total stopping time

Proposition 16.4, page 279 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $T$ be the accelerated Collatz map ($T(n)=n/2$ for even $n$ and $T(n)=(3n+1)/2$ for odd $n$) and let $\sigma_\infty(n)=\min\{k\ge0:T^{(k)}(n)=1\}$ be the total stopping time of a positive integer $n$ in steps of $T$ ($\infty$ if there is no such $k$). Let $n$ and $n'$ be positive integers whose orbits meet at a point $x\ge3$, reached from $n$ after $s$ steps and from $n'$ after $s'$ steps. Then $\sigma_\infty(n')=\sigma_\infty(n)-s+s'$. In particular, if $s=s'$ (as when the two orbits occupy the pair $8k+4$, $8k+5$ with $k\ge1$ at the same step, and then coincide three steps later), then $n'$ lies in every level set $\{\sigma_\infty>m\}$ that contains $n$.

Source: [`book/chapters/ch10c_samplers.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch10c_samplers.tex), label `prop:ch10b:acceptance`.

## Chapter 17. Progress on the Open Problems

<a id="stopping_time_cumulant_deficits_census_fit"></a>

### Additive deficits of the second and third cumulants of total stopping times, fitted to the exact census

Proposition 17.1, page 290 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $T$ be the accelerated Collatz map, $T(m)=m/2$ for even $m$ and $T(m)=(3m+1)/2$ for odd $m$, and let $\sigma_\infty(n)$ be the total stopping time of $n$, the number of steps of $T$ from $n$ to $1$ (the census computes it for every $n<2^{30}$). For each $k$ with $20\le k\le29$ let $n$ be uniform in the dyadic block $[2^k,2^{k+1})$, let $\kappa_2(\sigma_\infty)$ and $\kappa_3(\sigma_\infty)$ be the exact second and third cumulants of $\sigma_\infty(n)$, and let $\mathbb{E}\log n$ be the mean of $\log n$ over the block. Fitting $\kappa_2(\sigma_\infty)=a_2\,\mathbb{E}\log n-D_2$ and $\kappa_3(\sigma_\infty)=a_3\,\mathbb{E}\log n-D_3$ over these ten blocks gives slopes $a_2=101.06$ and $a_3=4430.2$, within $0.32$ and $0.13$ percent of the model values $A_2=\sigma^2/|\mu|^3=101.39$ and $A_3=3\sigma^4/|\mu|^5=4435.7$ of the Bernoulli model (drift $\mu=\log(\sqrt3/2)$ and step variance $\sigma^2=(\tfrac12\log3)^2$), and constant deficits $D_2=133.3$ and $D_3=34\,060$, in units of squared and cubed steps. The largest residuals of the fits over the ten blocks are $1.1$ and $89$.

Source: [`book/chapters/ch10e_progress.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch10e_progress.tex), label `prop:ch16:deficit`.

<a id="switching_mm_exact_convex_rate_conjecture"></a>

### Conjectured exact worst case rate of Switching MM on smooth convex functions, both branches

Conjecture 17.2, page 292 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $0\le\beta<1$, let $N\ge2$ be even, and run Switching MM with step $1/L$ and momentum coefficient $\beta$ from $x_0$ (from each even iterate $x_{2j}$ it takes a gradient step of length $(1+\beta)/L$ to a point $y$, then a gradient step of length $1/L$ from $y$ to $x_{2j+2}$). Put $D_N=(4+2\beta)N+2$ and $D'_N=4N-2+2(N-2)\beta$. Conjecture: over all convex $L$ smooth functions $f$ with a minimizer $x_*$, $f_*=f(x_*)$, and $\|x_0-x_*\|\le R$, the worst case is $\sup\bigl(f(x_N)-f_*\bigr)=LR^2\max\bigl\{1/D_N,\ \beta^2/D'_N\bigr\}$. The two branches cross at $\beta_N=\bigl(\sqrt{9N^2-2N+1}-N-1\bigr)/(2N)$, the positive root of $N\beta^2+(N+1)\beta-(2N-1)=0$, so that the main branch $LR^2/D_N$ is exact for all even $N$ exactly when $\beta\le\beta_2=(\sqrt{33}-3)/4=0.68614$.

Source: [`book/chapters/ch10e_progress.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch10e_progress.tex), label `conj:ch16:full`.

<a id="switching_mm_explicit_worst_functions"></a>

### Explicit worst case functions for Switching MM on smooth convex functions

Proposition 17.3, page 293 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $0\le\beta<1$, let $N\ge2$ be even, and run Switching MM with step $1/L$ and momentum coefficient $\beta$ from $x_0$ (from each even iterate $x_{2j}$ it takes a gradient step of length $(1+\beta)/L$ to a point $y$, then a gradient step of length $1/L$ from $y$ to $x_{2j+2}$). Put $D_N=(4+2\beta)N+2$ and $D'_N=4N-2+2(N-2)\beta$. Then the supremum of $f(x_N)-f_*$ over all convex $L$ smooth functions $f$ with a minimizer $x_*$, $f_*=f(x_*)$, and $\|x_0-x_*\|\le R$ satisfies $\sup\bigl(f(x_N)-f_*\bigr)\ge LR^2\max\bigl\{1/D_N,\ \beta^2/D'_N\bigr\}$.

Source: [`book/chapters/ch10e_progress.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch10e_progress.tex), label `prop:ch16:lower`.

<a id="switching_mm_exact_rate_two_and_four_steps"></a>

### Exact worst case rate of Switching MM after two and four steps, by an explicit certificate

Theorem 17.4, page 294 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $N\in\{2,4\}$ and $0\le\beta<1$. Run Switching MM with step $1/L$ and momentum coefficient $\beta$ in its gradient step form, $x_{k+1}=x_k-(s_k/L)\nabla f(x_k)$ for $0\le k<N$, with $s_k=1+\beta$ for even $k$ and $s_k=1$ for odd $k$ (a gradient step of length $(1+\beta)/L$ followed by one of length $1/L$, repeated; $x_N$ is the $N$th iterate of Switching MM). Let $f$ be convex and $L$ smooth with a minimizer $x_*$, $f_*=f(x_*)$, and $\|x_0-x_*\|\le R$. Put $D_N=(4+2\beta)N+2$, $W_k=2(s_0+\dots+s_k)$, $\lambda_k=W_k/\bigl(D_N-W_k+2\beta\,[k\text{ even}]\bigr)$ for $0\le k<N$ (where $[k\text{ even}]$ is $1$ for even $k$ and $0$ for odd $k$), $\mu_0=\lambda_0$, $\mu_k=\lambda_k-\lambda_{k-1}$ for $0<k<N$, and $\mu_N=1-\lambda_{N-1}$. Multiplying the smooth convex interpolation inequality between the points $x_k$ and $x_{k+1}$ by $\lambda_k$ and the one between the minimizer and $x_k$ by $\mu_k$, and summing, gives $f(x_N)-f_*\le LR^2/D_N-q$, where $q$ is a quadratic form in $x_0-x_*$ and $\nabla f(x_0),\dots,\nabla f(x_N)$. The form $q$ is nonnegative identically when $N=2$ and $0\le\beta\le b_2=0.64461$, or $N=4$ and $0\le\beta\le b_4=0.70438$, where $b_2$ and $b_4$ are the positive roots of $2\beta^4+11\beta^3+17\beta^2+\beta-11$ and $4\beta^4+21\beta^3+31\beta^2-\beta-23$. In these ranges, which contain the Collatz value $\beta=\log_2 3-1$, the exact worst case of $f(x_N)-f_*$ over this class is therefore $LR^2/D_N$, the main branch of the conjectured exact rate $LR^2\max\{1/D_N,\beta^2/D'_N\}$ with $D'_N=4N-2+2(N-2)\beta$.

Source: [`book/chapters/ch10e_progress.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch10e_progress.tex), label `thm:ch16:cert`.

<a id="switching_mm_noise_bound_correlated_noise"></a>

### The five thirds noise floor bound for Switching MM under correlated gradient noise on quadratics

Proposition 17.5, page 295 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: proved.

Statement. Let $f(x)=\tfrac12x^\top Ax-b^\top x$ with $A$ symmetric positive definite with largest eigenvalue $L$ and minimizer $x_*$, and suppose the algorithm sees the noisy gradient $g(x)=\nabla f(x)+\xi$, where a fresh $\xi$ is drawn at every call, independent of the past, with mean zero and an arbitrary covariance matrix $Q$ (instead of $\sigma^2I$). Run gradient descent $x_{t+1}=x_t-\alpha g(x_t)$ and Switching MM with momentum coefficient $\beta$, both with $g$ in place of $\nabla f$ and the same step $\alpha$ (Switching MM starts from $x_{-1}=x_0$, takes the plain step $x_{t+1}=x_t-\alpha g(x_t)$ at even $t$, and at odd $t$ takes the extrapolated step $x_{t+1}=y_t-\alpha g(y_t)$ with $y_t=x_t+\beta(x_t-x_{t-1})$). Then for every step $0<\alpha\le1/L$ and every $0\le\beta<1$, the stationary excess loss $\lim\mathbb{E}\bigl[f(x_t)-f(x_*)\bigr]$ of Switching MM, phase averaged over even and odd $t$, is at most $\tfrac53$ times that of gradient descent with the same step.

Source: [`book/chapters/ch10e_progress.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch10e_progress.tex), label `prop:ch16:correlated`.

## Chapter 19. Open Problems and Computational Projects

<a id="trajectory_mixing_of_syracuse_residues"></a>

### Trajectory mixing of Syracuse residues modulo powers of 3

Problem O2, page 315 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $S$ be the Syracuse map on the odd positive integers, $S(y)=(3y+1)/2^{v_2(3y+1)}$, where $v_2(m)$ is the exponent of the largest power of $2$ dividing $m$, and for an odd $x_0$ let $y_j=S^{(j)}(x_0)$. Fix $a\ge1$, let $U_a$ be the residues modulo $3^a$ prime to $3$, and let $\pi_a$ be the stationary law of the Markov chain on $U_a$ that moves from $r$ to $(3r+1)2^{-K}\bmod3^a$, with $K$ an independent Geometric$(\tfrac12)$ variable, $\mathbb{P}(K=k)=2^{-k}$ for $k\ge1$. Fix a constant $C>1/|\log(3/4)|\approx3.48$ and let $x_0$ be drawn uniformly from the odd integers up to $X$. Problem: show that the empirical distribution of $y_j\bmod3^a$ over the strides before the orbit reaches $1$, stopped after $\lfloor C\log X\rfloor$ strides if the orbit has not reached $1$ by then, converges in probability to $\pi_a$ as $X\to\infty$.

Source: [`book/chapters/ch11_open.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch11_open.tex), label `prob:ch11:o2`.

<a id="failure_fraction_of_potentials_of_growing_depth"></a>

### Failure fraction of Lyapunov potentials of growing depth for the Syracuse map

Problem O3$'$, page 316 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $S$ be the Syracuse map on the odd positive integers, $S(x)=(3x+1)/2^{v_2(3x+1)}$, where $v_2(m)$ is the exponent of the largest power of $2$ dividing $m$. A potential is a function $V(x)=\log x+\varphi(x)$ on the odd positive integers, and it fails at $x$ if $V(S(x))\ge V(x)$. Group the strides of an orbit into blocks, each a maximal run of strides $y\mapsto S(y)$ with $v_2(3y+1)=1$ followed by one stride with $v_2(3y+1)\ge2$; the block that starts at $x$ has $v_2(x+1)$ strides, its run length. No potential whose correction $\varphi$ depends on boundedly many binary digits of $x$ decreases at every stride, and the run length potentials $\log x+c\,v_2(x+1)$ fail at a fraction at least $0.143$ of the odd integers. Problem: for $h\ge1$, determine the infimum, over the potentials of depth $h$, whose correction reads the run lengths of the first $h$ blocks of the orbit of $x$, of the density of the odd integers at which the potential fails; decide whether it tends to $0$ as $h\to\infty$, and if so, how fast.

Source: [`book/chapters/ch11_open.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch11_open.tex), label `prob:ch11:o3p`.

<a id="bounded_variance_deficit_of_total_stopping_times"></a>

### Boundedness of the variance deficit and the third cumulant of total stopping times

Problem O4, page 316 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $T$ be the accelerated Collatz map, $T(n)=n/2$ for even $n$ and $T(n)=(3n+1)/2$ for odd $n$, and for $n\ge2$ whose orbit reaches $1$ let $\sigma_\infty(n)$ be the number of steps of $T$ from $n$ to $1$. For $k\ge1$ and $n$ in the dyadic block $[2^k,2^{k+1})$, let $a_{\mathrm t}(n)$ be the number of odd values among $T^{(k)}(n),\ldots,T^{(\sigma_\infty(n)-1)}(n)$ and $Z_{\mathrm t}(n)=\tfrac12\log3\,\bigl(2a_{\mathrm t}(n)-(\sigma_\infty(n)-k)\bigr)$, the centred parity sum over the steps after the first $k$. Put $\sigma^2=\tfrac14(\log3)^2$ and, with $n$ uniform on the block, define the additive variance deficit $D_{\mathrm v}(k)=\sigma^2\,\mathbb{E}(\sigma_\infty-k)-\operatorname{Var}Z_{\mathrm t}$. Problem: prove that $D_{\mathrm v}(k)$ is bounded in $k$, either assuming that every orbit reaches $1$, or unconditionally for the integers of each block whose orbits reach $1$ within a fixed multiple of $k$ steps; and explain the third cumulant: the skewness of $\sigma_\infty$ on the block is $0.67$ at $k=29$, against $0.96$ for the inverse Gaussian first passage law of the model.

Source: [`book/chapters/ch11_open.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch11_open.tex), label `prob:ch11:o4`.

<a id="arithmetic_sharpening_of_the_cycle_window"></a>

### Arithmetic sharpening of the window bound on the least element of a Collatz cycle

Problem C1, page 317 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. First stated in Polson and Sokolov (2026), Five Bayesian problems raised by AlphaProof Nexus, manuscript.

Status: open.

Statement. Let $T$ be the accelerated Collatz map, $T(n)=n/2$ for even $n$ and $T(n)=(3n+1)/2$ for odd $n$. For a cycle of $T$ on the positive integers with $m$ odd steps, $k$ steps in all and least element $x_{\min}$, put $\Lambda=k\log2-m\log3>0$. Every nontrivial cycle satisfies $x_{\min}<m/(3\Lambda)$ (the window lemma of the book). Problem: find a function $f$, derived from the arithmetic of the parity vector of the cycle rather than from the size of its elements, such that every nontrivial cycle satisfies $x_{\min}<m/(3\Lambda f(m))$. A function with $f(m)\ge\frac m3(1+\eta)$ for a fixed $\eta>0$ and all large $m$ would confine the least element of every nontrivial cycle with many odd steps below a constant, about $1/(2\eta)$ for small $\eta$, and so reduce the cycle half of the Collatz conjecture to a finite verification.

Source: [`book/chapters/ch11_open.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch11_open.tex), label `prob:ch11:c1`.

<a id="verification_ranges_for_later_rungs_under_hercher_accounting"></a>

### Verification ranges that exclude the later rungs under Hercher's sharper accounting

Problem C3, page 318 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $\theta=\log_23$ and, for an integer $r\ge1$, $\delta(r)=\lceil r\theta\rceil-r\theta$; the rungs of the one sided ladder of $\theta$ are the $r$ with $\delta(r)<\delta(r')$ for every $1\le r'<r$, the best one sided approximations of $\theta$ from above. Let $T$ be the accelerated Collatz map, $T(n)=n/2$ for even $n$ and $T(n)=(3n+1)/2$ for odd $n$, and call $X$ a verification range if every positive integer below $X$ reaches $1$ under $T$. Hercher (2023, Theorem 27) shows that if $X$ is a verification range, a nontrivial cycle of $T$ with $r$ odd elements and $s$ elements in all satisfies $0<s/r-\theta<1/(4X\log2)$, so that the least cycle length it does not exclude is a rung, and his computer assisted Corollary 29 shows that the verification range $3\cdot2^{69}$ already excludes cycles whose length is the rung $72\,057\,431\,991$. Problem: extend that computation to the later rungs: for each rung from $137\,528\,045\,312$ on, find the least verification range at which Hercher's accounting excludes it.

Source: [`book/chapters/ch11_open.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch11_open.tex), label `prob:ch11:c3`.

<a id="bounded_form_of_tao_theorem"></a>

### Bounded form of Tao's theorem: orbit minimum at most 2^71 on a set of logarithmic density one

Problem B1, page 318 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $T_0$ be the Collatz map on the positive integers, $T_0(n)=n/2$ for even $n$ and $T_0(n)=3n+1$ for odd $n$, and let $\mathrm{Col}_{\min}(n)=\min_{k\ge0}T_0^{(k)}(n)$ be the least value on the orbit of $n$. A set $A$ of positive integers has logarithmic density one if $(\log x)^{-1}\sum_{n\le x,\,n\in A}1/n\to1$ as $x\to\infty$. Tao proved that for every function $f$ with $f(n)\to\infty$ one has $\mathrm{Col}_{\min}(n)<f(n)$ for all $n$ in a set of logarithmic density one, and every positive integer up to $2^{71}$ is known to reach $1$ (verification of Bařina). Problem: show that $\mathrm{Col}_{\min}(n)\le2^{71}$ for all $n$ in a set of logarithmic density one. By the verification up to $2^{71}$, if the positive density theorems of Mazur and Lau are right, this implies the Collatz conjecture, and the conjecture implies it.

Source: [`book/chapters/ch11_open.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch11_open.tex), label `prob:ch11:b1`.

<a id="true_decay_rate_of_syracuse_fourier_coefficients"></a>

### True decay rate of the Fourier coefficients of the Syracuse random variable modulo powers of 3

Problem B2, page 318 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $a_1,a_2,\ldots$ be independent random variables with $\mathbb{P}(a_i=k)=2^{-k}$ for $k\ge1$, and let $\pi_d$ be the law on $\mathbb{Z}/3^d\mathbb{Z}$ of $\sum_{i=1}^{d}3^{i-1}2^{-(a_1+\cdots+a_i)}$, with $2$ inverted modulo $3^d$ (the law of Tao's Syracuse random variable modulo $3^d$). Its Fourier coefficient at the frequency $\xi$ is $\widehat\pi_d(\xi)=\sum_{y\bmod3^d}\pi_d(y)e^{-2\pi i\xi y/3^d}$. Tao proved that $\max_{3\nmid\xi}|\widehat\pi_d(\xi)|$ eventually decays faster than any power of $d$, and the exact computation of the book shows that through $d=18$ it stays between $3.6\,d^{-2}$ and $3.9\,d^{-2}$ from $d=7$ on and sits at $\xi=\pm2^k$ with $d\le k\le1.28\,d$. Problem: determine how fast $\max_{3\nmid\xi}|\widehat\pi_d(\xi)|$ tends to zero as $d\to\infty$. Is it $\exp(-c(\log d)^2)$, $\exp(-cd^{\gamma})$ for some $\gamma<1$, or eventually geometric?

Source: [`book/chapters/ch11_open.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch11_open.tex), label `prob:ch11:b2`.

<a id="sharpness_of_shaik_threshold_for_stopping_times"></a>

### Sharpness of Shaik's threshold 3/log(4/3) for total stopping times

Problem B3, page 319 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $T_0$ be the Collatz map on the positive integers, $T_0(n)=n/2$ for even $n$ and $T_0(n)=3n+1$ for odd $n$, and let $\sigma_0(n)$ be the total stopping time of $n$, the least $k\ge0$ with $T_0^{(k)}(n)=1$ (infinite if there is none). Shaik's preprint (2026) obtains positive lower density for convergence within $c\log x$ steps of $T_0$, for the integers below $x$, for every $c>3/\log(4/3)\approx10.428$. Problem: show that for every $c<3/\log(4/3)$ the set of $n$ with $\sigma_0(n)\le c\log n$ has natural density zero, so that Shaik's threshold is sharp. A renewal calculation in the book predicts that below the threshold the share of $n$ meeting the bound tends to zero.

Source: [`book/chapters/ch11_open.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch11_open.tex), label `prob:ch11:b3`.

<a id="exact_function_gap_rate_of_switching_mm_for_small_condition_numbers"></a>

### Exact worst case function gap rate of Switching MM for small condition numbers

Problem S1, page 319 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $f:\mathbb{R}^d\to\mathbb{R}$ be $\mu$ strongly convex and $L$ smooth, with $\kappa=L/\mu$, let $0\le\beta<1$, and let $G=\nabla f$. Switching MM with deterministic parity injection and step $\alpha$, started at $x_{-1}=x_0$, takes the gradient step $x_{t+1}=x_t-\alpha G(x_t)$ at even $t$ and, at odd $t$, $y_t=x_t+\beta(x_t-x_{t-1})$ and $x_{t+1}=y_t-\alpha G(y_t)$, here with $\alpha=1/L$. Its even iterates contract the distance to the minimizer by $q=(1-1/\kappa)\max\bigl(1-(1+\beta)/\kappa,\beta\bigr)$ per period of two steps, which is the exact worst case over the class, and the function gap contracts by $q^2$ per period, which is exact for $\kappa\ge(1+\beta)/(1-\beta)$. Problem: for $\kappa<(1+\beta)/(1-\beta)$, where performance estimation gives a smaller worst case (for example $0.049$ against $q^2=0.086$ at $\beta=\log_23-1$ and $\kappa=2$), find the exact worst case contraction of the function gap per period in closed form and describe the function that attains it.

Source: [`book/chapters/ch11_open.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch11_open.tex), label `prob:ch11:s1`.

<a id="five_thirds_noise_bound_of_switching_mm_beyond_quadratics"></a>

### Five thirds noise bound of Switching MM beyond quadratics

Problem S2, page 319 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $f:\mathbb{R}^d\to\mathbb{R}$ be $\mu$ strongly convex and $L$ smooth, accessed through a stochastic gradient oracle $G(x)=\nabla f(x)+\xi$ whose noise $\xi$ is fresh at every call, has mean zero and bounded variance, and let $0\le\beta<1$. Switching MM with deterministic parity injection and step $\alpha$, started at $x_{-1}=x_0$, takes the gradient step $x_{t+1}=x_t-\alpha G(x_t)$ at even $t$ and, at odd $t$, $y_t=x_t+\beta(x_t-x_{t-1})$ and $x_{t+1}=y_t-\alpha G(y_t)$. For quadratic objectives with additive noise of covariance $\sigma^2I$, the book proves that on every eigendirection with normalized step at most $1$ the stationary variance of this method is at most $\beta+2/(2+\beta)<\tfrac53$ times that of gradient descent with the same step. Problem: show that for every such $f$ and every $\beta<1$, the stationary mean square error of this method is at most $\tfrac53$ times that of gradient descent with the same step, up to terms that vanish with the step size.

Source: [`book/chapters/ch11_open.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch11_open.tex), label `prob:ch11:s2`.

<a id="adaptive_switching_rule_for_switching_mm"></a>

### A principled adaptive switching rule for Switching MM

Problem S4, page 320 of Polson and Sokolov (2026), *On the Collatz Conjecture: Martingales, Momentum and Continued Fractions*, edition 1. New in this book.

Status: open.

Statement. Let $f$ be an $L$ smooth objective. Switching MM alternates plain gradient steps $x_{t+1}=x_t-\alpha\nabla f(x_t)$ with momentum injection steps $y_t=x_t+\beta(x_t-x_{t-1})$, $x_{t+1}=y_t-\alpha\nabla f(y_t)$; with deterministic parity injection (an injection at every odd $t$) its stationary variance under additive gradient noise is, on quadratics and at steps up to $1/L$, at most $\tfrac53$ times that of gradient descent with the same step, uniformly in $0\le\beta<1$, but in the absence of noise it converges more slowly than momentum applied at every step, and its adaptive variant injects only when the most recent reduction of a merit gap was inadequate. Problem: find a principled switching rule, using only $L$ and observed function values, that retains that uniform noise bound while closing the noiseless gap to momentum applied at every step.

Source: [`book/chapters/ch11_open.tex`](https://github.com/VadimSokolov/collatz/blob/v1.1.0/book/chapters/ch11_open.tex), label `prob:ch11:s4`.

***

<small>License: Text, figures and PDF: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: [MIT](https://opensource.org/license/mit).</small>

<small>This page is generated from the book's repository ([main](https://github.com/VadimSokolov/collatz/tree/main)). The numbers and pages are those of the PDF, [version 1.1.0](https://github.com/VadimSokolov/collatz/releases/tag/v1.1.0).</small>

<!-- Generated by bookrepo.py from book.yaml and results.yaml. Edit those files, not this one. -->
