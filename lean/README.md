# Lean files

The two files of the section "Machine checked statements" of the book (chapter 10), in Lean 4 with the core library only (no Mathlib):

- `CollatzProved.lean`: the accelerated map `T`, its iterates `Tn` and the number of odd steps `odd`, and the theorems B1 (the affine identity), B2 (the cycle inequality $3^m<2^k$), B3 (every starting value from 1 to 10^4 reaches 1 within 400 steps) and B4 (the trivial cycle).
- `CollatzOpen.lean`: the open statements O1 (the conjecture), O3 (its descent form) and O4 (no nontrivial cycles), each closed by `sorry`.

The book takes these statements from the authors' note on AlphaProof Nexus. The files here were written again in October 2026 from the statements the book prints, which they reproduce character for character, with proofs that follow the routes the book describes: B1 by induction through a stronger lemma that names the remainder by its recursion, B2 from B1, B3 by a Boolean evaluator with a soundness lemma, run natively, and B4 by induction.

## Checking

With [elan](https://github.com/leanprover/elan) installed, in this folder (`lean-toolchain` selects Lean 4.15.0):

```bash
lean CollatzProved.lean   # no output: no errors, no warnings
lean CollatzOpen.lean     # one warning per open statement: declaration uses 'sorry'
```

## Axioms

`#print axioms` with Lean 4.15.0 reports: B1 depends on `propext` and `Quot.sound`; B2 also on `Classical.choice`; B3 on `propext`, `Quot.sound` and `Lean.ofReduceBool`, since its finite check is evaluated natively, which amounts to trusting the compiled evaluator; B4 on no axioms.
