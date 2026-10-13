# Insert \index{} tags at the first prose occurrence of key terms in each section.
import re, glob, sys
TERMS = [  # (index key, regex)
 ('Collatz conjecture', r'Collatz conjecture'), ('Syracuse map', r'Syracuse map'),
 ('accelerated map', r'accelerated map'), ('standard map', r'standard map'),
 ('parity vector', r'parity vectors?'), ('Terras, R.', r'Terras'), ('stopping time', r'stopping times?'),
 ('coefficient stopping time', r'coefficient stopping time'), ('martingale', r'martingales?'),
 ('supermartingale', r'supermartingales?'), ('Jensen inequality', r"Jensen's inequality|Jensen gap"),
 ('Wald identity', r"Wald's (?:second )?identit(?:y|ies)"), ('Ville inequality', r'Ville'),
 ('renewal theory', r'renewal'), ('autocorrelation', r'autocorrelations?'),
 ('uniform integrability', r'uniform(?:ly)? integrab(?:le|ility)'), ('Lyapunov function', r'Lyapunov functions?'),
 ('Foster criterion', r"Foster's criterion"), ('majorization', r'majoriz(?:er|ers|ation)'),
 ('MM algorithm', r'MM (?:algorithm|scheme|framework|descent)'), ('heavy ball', r'heavy ball'),
 ('Nesterov acceleration', r'Nesterov'), ('momentum', r'momentum'), ('condition number', r'condition numbers?'),
 ('Switching MM', r'Switching MM'), ('gain weight theorem', r'gain weight'), ('Chebyshev step sizes', r'Chebyshev'),
 ('cyclic step sizes', r'cyclic step'), ('noise floor', r'noise floors?'), ('multiplicative noise', r'multiplicative noise|curvature noise'),
 ('Pareto frontier', r'Pareto'), ('joint spectral radius', r'joint spectral radius'), ('switched system', r'switched systems?'),
 ('performance estimation', r'performance estimation'), ('modified equation', r'modified equations?'),
 ('implicit regularization', r'implicit (?:gradient )?regularization'), ('difference inequality', r'difference inequalit(?:y|ies)'),
 ('inverse tree', r'inverse tree'), ('Krasikov, I.', r'Krasikov'), ('Lagarias, J. C.', r'Lagarias'), ('Tao, T.', r'Tao'),
 ('Eliahou, S.', r'Eliahou'), ('Hercher, C.', r'Hercher'), ('Steiner, R. P.', r'Steiner'), ('Baker method', r"Baker's (?:method|theory|theorem)"),
 ('continued fraction', r'continued fractions?'), ('convergent', r'convergents?'), ('semiconvergent', r'semiconvergents?|intermediate fractions?'),
 ('mediant lemma', r'mediant lemma'), ('one sided approximation', r'one sided (?:best )?approximations?'),
 ('irrationality measure', r'irrationality measures?'), ('window lemma', r'window lemma'), ('cycle identity', r'cycle identity'),
 ('Erdos problem 125@Erd\\H{o}s problem 125', r'Erd\\H\{o\}s \$?125'), ('AlphaProof Nexus', r'AlphaProof Nexus'), ('Lean', r'\bLean\b'),
 ('certificate', r'certificates?'), ('verification range', r'verification ranges?|verified range'),
 ('Cesaro identity@Ces\\`aro identity', r"Ces\\`aro identity"), ('Cesaro conjecture@Ces\\`aro conjecture', r"Ces\\`aro (?:liminf )?conjecture"),
 ('equidistribution', r'equidistribut(?:ion|e)'), ('spectral gap', r'spectral gap'), ('mixing', r'\bmix(?:es|ing)\b'),
 ('Parrondo paradox', r'Parrondo'), ('Kullback and Leibler divergence', r'KL divergence|Kullback'), ('Conway, J. H.', r'Conway'),
 ('FRACTRAN', r'FRACTRAN'), ('Barina, D.', r'Bar\\v\{i\}na'), ('balanced words', r'balanced words?|Christoffel'),
 ('rewriting system', r'rewriting'), ('Bernoulli model', r'Bernoulli model|Bernoulli\$\(\\tfrac12\)\$ model'),
 ('drift', r'\bdrift\b'), ('central limit theorem', r'central limit theorem'), ('large deviations', r'large deviations?|Chernoff'),
 ('2 adic integers', r'\$2\$ adic|\b2 adic'), ('Haar measure', r'Haar measure'), ('density theorem', r'density theorem'),
 ('residue sieve', r'residue (?:sieve|certificates?)'), ('negative cycles', r'negative (?:integers|cycles)'),
 ('rational cycles', r'rational cycles?'), ('Diophantine approximation', r'Diophantine'),
 ('Ornstein', r'(?!x)x'),
]
MASK_CMDS = r'\\(?:label|cref|Cref|ref|eqref|cite|citep|citet|includegraphics|texttt|url|input|index|chapter|section|subsection|paragraph|caption|begin|end|bibliography)\*?(?:\[[^\]]*\])?\{[^{}]*\}'
def mask_line(line):
    m = [False] * len(line)
    for mm in re.finditer(MASK_CMDS, line):
        for i in range(mm.start(), mm.end()): m[i] = True
    inmath = False; i = 0
    while i < len(line):
        if line[i] == '\\' and i + 1 < len(line): 
            if inmath: m[i] = m[i+1] = True
            i += 2; continue
        if line[i] == '$': inmath = not inmath; m[i] = True
        elif inmath: m[i] = True
        if line[i] == '%': 
            for j in range(i, len(line)): m[j] = True
            break
        i += 1
    return m
SKIP_ENVS = ('equation', 'align', 'gather', 'multline', 'lstlisting', 'figure', 'table', 'tabular', 'algorithm', 'algorithmic', 'array', 'longtable', 'center')
for f in sys.argv[1:]:
    L = open(f).read().split('\n'); out = []; done = set(); depth = []
    for line in L:
        st = line.strip()
        for e in re.findall(r'\\begin\{(\w+)\*?\}', line):
            if e in SKIP_ENVS: depth.append(e)
        if re.match(r'\\(chapter|section|subsection)', st): done = set()
        skip = bool(depth) or st.startswith('\\') and not st.startswith('\\emph') and not st.startswith('\\Cref') and not st.startswith('\\cref') \
               or '\\[' in line or '\\]' in line or '\\caption' in line or 'begin{' in line and '[' in line
        if not skip:
            msk = mask_line(line); ins = []
            for key, pat in TERMS:
                if key in done: continue
                for mm in re.finditer(pat, line):
                    if any(msk[mm.start():mm.end()]): continue
                    ins.append((mm.end(), key)); done.add(key); break
            for pos, key in sorted(ins, reverse=True):
                line = line[:pos] + '\\index{' + key + '}' + line[pos:]
        out.append(line)
        for e in re.findall(r'\\end\{(\w+)\*?\}', line):
            if e in SKIP_ENVS and depth and depth[-1] == e: depth.pop()
    open(f, 'w').write('\n'.join(out))
print('done')
