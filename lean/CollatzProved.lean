/-
The proved statements of the section "Machine checked statements" of the book
"On the Collatz Conjecture" (N. Polson and V. Sokolov, 2026): B1 to B4.
Lean 4.15 core library only, no Mathlib. Check with `lean CollatzProved.lean`.
-/

def T (x : Nat) : Nat :=
  if x % 2 = 0 then x / 2 else (3 * x + 1) / 2
def Tn : Nat → Nat → Nat
  | 0, x => x
  | k + 1, x => T (Tn k x)
def odd : Nat → Nat → Nat
  | 0, _ => 0
  | k + 1, x => odd k x + (if Tn k x % 2 = 1 then 1 else 0)

/-- The remainder of the affine identity, by the recursion `c (k+1) = 3 c k + 2^k` after an
odd step and `c (k+1) = c k` after an even one. -/
def rem : Nat → Nat → Nat
  | 0, _ => 0
  | k + 1, x => if Tn k x % 2 = 1 then 3 * rem k x + 2 ^ k else rem k x

/-- The stronger form of B1, with the remainder named. -/
theorem affine (k x : Nat) :
    2 ^ k * Tn k x = 3 ^ (odd k x) * x + rem k x ∧ (odd k x > 0 → rem k x > 0) := by
  induction k with
  | zero => simp [Tn, odd, rem]
  | succ k ih =>
    obtain ⟨ih1, ih2⟩ := ih
    have hT : Tn (k + 1) x = T (Tn k x) := rfl
    cases Nat.mod_two_eq_zero_or_one (Tn k x) with
    | inl h0 =>
      have hT' : Tn (k + 1) x = Tn k x / 2 := by rw [hT, T, if_pos h0]
      have ho : odd (k + 1) x = odd k x := by simp [odd, h0]
      have hr : rem (k + 1) x = rem k x := by simp [rem, h0]
      rw [hT', ho, hr]
      refine ⟨?_, ih2⟩
      have h2 : 2 * (Tn k x / 2) = Tn k x := by omega
      calc 2 ^ (k + 1) * (Tn k x / 2) = 2 ^ k * (2 * (Tn k x / 2)) := by
            rw [Nat.pow_succ, Nat.mul_assoc]
        _ = 3 ^ (odd k x) * x + rem k x := by rw [h2, ih1]
    | inr h1 =>
      have hne : ¬ (Tn k x % 2 = 0) := by omega
      have hT' : Tn (k + 1) x = (3 * Tn k x + 1) / 2 := by rw [hT, T, if_neg hne]
      have ho : odd (k + 1) x = odd k x + 1 := by simp [odd, h1]
      have hr : rem (k + 1) x = 3 * rem k x + 2 ^ k := by simp [rem, h1]
      rw [hT', ho, hr]
      refine ⟨?_, fun _ => by have := Nat.two_pow_pos k; omega⟩
      have h2 : 2 * ((3 * Tn k x + 1) / 2) = 3 * Tn k x + 1 := by omega
      have e1 : 2 ^ k * (3 * Tn k x) = 3 * (2 ^ k * Tn k x) := Nat.mul_left_comm _ _ _
      have e2 : 3 ^ (odd k x + 1) * x = 3 * (3 ^ (odd k x) * x) := by
        rw [Nat.pow_succ, Nat.mul_comm (3 ^ odd k x) 3, Nat.mul_assoc]
      calc 2 ^ (k + 1) * ((3 * Tn k x + 1) / 2)
          = 2 ^ k * (2 * ((3 * Tn k x + 1) / 2)) := by rw [Nat.pow_succ, Nat.mul_assoc]
        _ = 2 ^ k * (3 * Tn k x) + 2 ^ k := by rw [h2, Nat.mul_add, Nat.mul_one]
        _ = 3 ^ (odd k x + 1) * x + (3 * rem k x + 2 ^ k) := by rw [e1, e2, ih1]; omega

theorem B1 (k x : Nat) : ∃ c : Nat,
    2 ^ k * Tn k x = 3 ^ (odd k x) * x + c ∧
    (odd k x > 0 → c > 0) :=
  ⟨rem k x, affine k x⟩

theorem B2 {k x : Nat} (hx : x > 0) (hk : k > 0)
    (hcyc : Tn k x = x) (ho : odd k x > 0) :
    3 ^ (odd k x) < 2 ^ k := by
  have _ := hk   -- hx and hk are not needed: odd k x > 0 forces k > 0, and the identity of B1 forces x > 0
  have _ := hx
  obtain ⟨c, h, hc⟩ := B1 k x
  have hc' := hc ho
  rw [hcyc] at h
  have hlt : 3 ^ (odd k x) * x < 2 ^ k * x := by omega
  exact Nat.lt_of_mul_lt_mul_right hlt

/-- `reaches f x` is true when the orbit of `x` under `T` meets `1` within `f` steps. -/
def reaches : Nat → Nat → Bool
  | 0, x => x == 1
  | f + 1, x => x == 1 || reaches f (T x)

theorem Tn_succ' (k x : Nat) : Tn (k + 1) x = Tn k (T x) := by
  induction k generalizing x with
  | zero => rfl
  | succ k ih => show T (Tn (k + 1) x) = T (Tn k (T x)); rw [ih]

theorem reaches_sound (f x : Nat) (h : reaches f x = true) : ∃ j, j ≤ f ∧ Tn j x = 1 := by
  induction f generalizing x with
  | zero =>
    simp only [reaches, beq_iff_eq] at h
    exact ⟨0, Nat.le_refl 0, h⟩
  | succ f ih =>
    simp only [reaches, Bool.or_eq_true, beq_iff_eq] at h
    cases h with
    | inl h1 => exact ⟨0, Nat.zero_le _, h1⟩
    | inr h2 =>
      obtain ⟨j, hj, hT⟩ := ih (T x) h2
      exact ⟨j + 1, Nat.succ_le_succ hj, by rw [Tn_succ']; exact hT⟩

/-- `checkAll n` is true when every `x` with `1 ≤ x ≤ n` reaches `1` within `400` steps. -/
def checkAll : Nat → Bool
  | 0 => true
  | n + 1 => reaches 400 (n + 1) && checkAll n

theorem checkAll_sound (n : Nat) (h : checkAll n = true) :
    ∀ x, 1 ≤ x → x ≤ n → reaches 400 x = true := by
  induction n with
  | zero => intro x h1 h2; exact absurd h2 (by omega)
  | succ n ih =>
    intro x h1 h2
    simp only [checkAll, Bool.and_eq_true] at h
    cases Nat.lt_or_ge x (n + 1) with
    | inl hlt => exact ih h.2 x h1 (by omega)
    | inr hge => have hx : x = n + 1 := by omega
                 rw [hx]; exact h.1

theorem B3 : ∀ x, 1 ≤ x → x ≤ 10000 →
    ∃ j, j ≤ 400 ∧ Tn j x = 1 := by
  intro x h1 h2
  exact reaches_sound 400 x (checkAll_sound 10000 (by native_decide) x h1 h2)

theorem B4 : Tn 1 1 = 2 ∧ Tn 2 1 = 1 ∧
    ∀ k, Tn k 1 = 1 ∨ Tn k 1 = 2 := by
  refine ⟨rfl, rfl, ?_⟩
  intro k
  induction k with
  | zero => exact Or.inl rfl
  | succ k ih =>
    cases ih with
    | inl h => exact Or.inr (show T (Tn k 1) = 2 by rw [h]; rfl)
    | inr h => exact Or.inl (show T (Tn k 1) = 1 by rw [h]; rfl)
