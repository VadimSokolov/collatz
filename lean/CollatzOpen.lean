/-
The open statements of the section "Machine checked statements" of the book
"On the Collatz Conjecture" (N. Polson and V. Sokolov, 2026): the conjecture, its descent
form and the cycle problem. Lean 4.15 core library only, no Mathlib. Each `sorry` marks an
open problem, so `lean CollatzOpen.lean` prints one warning per statement.
-/

def T (x : Nat) : Nat :=
  if x % 2 = 0 then x / 2 else (3 * x + 1) / 2
def Tn : Nat → Nat → Nat
  | 0, x => x
  | k + 1, x => T (Tn k x)

theorem O1_collatz : ∀ x, x ≥ 1 → ∃ k, Tn k x = 1 := sorry
theorem O3_descent : ∀ x, x ≥ 2 → ∃ k, k > 0 ∧ Tn k x < x := sorry
theorem O4_no_cycles : ¬ ∃ x k, x > 2 ∧ k > 0 ∧ Tn k x = x := sorry
