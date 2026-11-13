import Mathilda

/-!
  Tests for the Mathilda tactics.  Run from `lean/` with
  `MATHILDA_BIN=../Mathilda lake env lean test/Tactics.lean`
  (see `tests/scripts/lean_tactics.sh`).  `#guard_msgs` makes any unexpected
  message an error, so a clean exit means every check passed.
-/

/-! ## `mathilda_linear_combination` (ideal membership) -/

theorem lc_circle (x y : Rat) (h1 : x^2 + y^2 = 1) (h2 : x = y) : 2*y^2 = 1 := by
  mathilda_linear_combination

theorem lc_int (a b c : Int) (h1 : a + b = c) (h2 : a - b = 0) : 2*a = c := by
  mathilda_linear_combination

theorem lc_cast (x : Int) (n : Nat) (h : x = n) : x^2 = (n : Int)^2 := by
  mathilda_linear_combination

/-! ## `mathilda_sos` (sums of squares) -/

theorem sos_amgm (x y : Rat) : 2*x*y ≤ x^2 + y^2 := by
  mathilda_sos

theorem sos_quartic (x y : Int) : x^4 - 2*x^2 + 1 + y^2 ≥ 0 := by
  mathilda_sos

theorem sos_hyp (x : Rat) (h : 0 ≤ x) : 0 ≤ x^3 := by
  mathilda_sos

/-! ## `mathilda_prime` (Pratt certificates) -/

theorem prime_274177 : ∃ c, Mathilda.PrattCertifies 274177 c := by
  mathilda_prime

theorem prime_1000003 : ∃ c, Mathilda.PrattCertifies 1000003 c := by
  mathilda_prime

/-! ## Only the standard axioms (no `sorryAx`, no `native_decide`) -/

/-- info: 'lc_circle' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms lc_circle
/-- info: 'lc_int' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms lc_int
/-- info: 'lc_cast' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms lc_cast
/-- info: 'sos_amgm' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms sos_amgm
/-- info: 'sos_quartic' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms sos_quartic
/-- info: 'sos_hyp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms sos_hyp
/-- info: 'prime_274177' does not depend on any axioms -/
#guard_msgs in #print axioms prime_274177
/-- info: 'prime_1000003' does not depend on any axioms -/
#guard_msgs in #print axioms prime_1000003

/-! ## `#mathilda` (ask the CAS; proves nothing) -/

/-- info: (x - 1) * (1 + x) * (1 - x + x ^ 2 - x ^ 3 + x ^ 4) * (1 + x + x ^ 2 + x ^ 3 + x ^ 4) -/
#guard_msgs in #mathilda "Factor[x^10-1]"

/-! ## Negative: false goals are not proved -/

/--
error: mathilda: Mathilda found no certificate: the goal is not a polynomial combination of the hypotheses (no cofactors up to degree 4)
-/
#guard_msgs (error) in
example (x y : Rat) (h1 : x^2 + y^2 = 1) (h2 : x = y) : 2*y^2 = 3 := by
  mathilda_linear_combination

/-- error: mathilda: Mathilda found no certificate: no sum-of-squares certificate found -/
#guard_msgs (error) in
example (x : Rat) : x ≤ 1 := by
  mathilda_sos

/-- error: mathilda: Mathilda found no certificate: 274179 is not prime -/
#guard_msgs (error) in
example : ∃ c, Mathilda.PrattCertifies 274179 c := by
  mathilda_prime
