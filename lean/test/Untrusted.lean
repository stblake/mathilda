import Mathilda

/-!
  The oracle is untrusted.  Run with `MATHILDA_BIN=test/bad_oracle.sh`: a fake
  Mathilda that answers every query with a certificate for a *different* (true)
  statement, or -- for goals carrying a marker constant 1001, 1002, ... -- with
  a hostile or malformed answer.  Each tactic must fail; none of these false
  goals is proved, and nothing in an answer is run.  `#guard_msgs` pins every
  rejection, so a clean exit means all of them were caught.
-/

/-! ## Well-formed certificates for the wrong statement: Lean rejects them -/

/--
error: mathilda_linear_combination: Lean rejected the certificate returned by Mathilda (cofactors [1,
 -y - x]); Mathilda is an untrusted oracle and its answer does not check
-/
#guard_msgs (error) in
example (x y : Rat) (h1 : x^2 + y^2 = 1) (h2 : x = y) : 2*y^2 = 3 := by
  mathilda_linear_combination

/--
error: mathilda_sos: Lean rejected the certificate returned by Mathilda (squares [(1, x -
   y, 0)]); Mathilda is an untrusted oracle and its answer does not check
-/
#guard_msgs (error) in
example (x y : Rat) : 3*x*y ≤ x^2 + y^2 := by
  mathilda_sos

/--
error: mathilda_prime: Lean rejected the certificate returned by Mathilda (Pratt certificate [(3, 2, [(2, 1)]), (7, 3, [(2, 1), (3, 1)]), (17, 3, [(2, 4)]), (274177, 5, [(2, 8), (3, 2), (7, 1), (17, 1)])]); Mathilda is an untrusted oracle and its answer does not check
-/
#guard_msgs (error) in
example : ∃ c, Mathilda.PrattCertifies 274179 c := by
  mathilda_prime

/-! ## Hostile and malformed answers: never elaborated, a parse error -/

-- The reviewer's exploit: Lean code (`by_elab`) that adds a kernel-unchecked
-- theorem and writes a file.  It is not a certificate, so it is refused unread.
/--
error: mathilda_sos: Lean could not parse the certificate returned by Mathilda (the identifier `by_elab` is not a variable v1, v2, …); Mathilda is an untrusted oracle and its answer is not used
-/
#guard_msgs (error) in
example (x y : Rat) : y ≤ x + 1001 := by
  mathilda_sos

-- The by_elab answer never ran: no theorem `bogus.evil` was added ...
/-- info: false -/
#guard_msgs in
#eval show Lean.CoreM Bool from return (← Lean.getEnv).contains `bogus.evil

-- ... and its file write did not happen.
/-- info: false -/
#guard_msgs in
#eval System.FilePath.pathExists "/tmp/mathilda_untrusted_pwned"

/--
error: mathilda_sos: Lean could not parse the certificate returned by Mathilda (the identifier `sorry` is not a variable v1, v2, …); Mathilda is an untrusted oracle and its answer is not used
-/
#guard_msgs (error) in
example (x : Rat) : x ≤ x + 1002 - 1 := by mathilda_sos   -- `sorry`

/--
error: mathilda_sos: Lean could not parse the certificate returned by Mathilda (the identifier `x` is not a variable v1, v2, …); Mathilda is an untrusted oracle and its answer is not used
-/
#guard_msgs (error) in
example (x : Rat) : x ≤ x + 1003 - 1 := by mathilda_sos   -- identifier `x`

/--
error: mathilda_sos: Lean could not parse the certificate returned by Mathilda (exponent 100000 exceeds 10000 at `100000`); Mathilda is an untrusted oracle and its answer is not used
-/
#guard_msgs (error) in
example (x : Rat) : x ≤ x + 1004 - 1 := by mathilda_sos   -- exponent 100000

/--
error: mathilda_sos: Lean could not parse the certificate returned by Mathilda (unexpected character '$'); Mathilda is an untrusted oracle and its answer is not used
-/
#guard_msgs (error) in
example (x : Rat) : x ≤ x + 1005 - 1 := by mathilda_sos   -- malformed tokens

/--
error: mathilda_sos: Lean could not parse the certificate returned by Mathilda (a numeral has more than 2000 digits); Mathilda is an untrusted oracle and its answer is not used
-/
#guard_msgs (error) in
example (x : Rat) : x ≤ x + 1006 - 1 := by mathilda_sos   -- 2500-digit numeral

/--
error: mathilda_sos: Lean could not parse the certificate returned by Mathilda (v7 names no variable of the goal); Mathilda is an untrusted oracle and its answer is not used
-/
#guard_msgs (error) in
example (x : Rat) : x ≤ x + 1007 - 1 := by mathilda_sos   -- v7: no such variable

/--
error: mathilda_prime: Lean could not parse the certificate returned by Mathilda (the identifier `by` is not a variable v1, v2, …); Mathilda is an untrusted oracle and its answer is not used
-/
#guard_msgs (error) in
example : ∃ c, Mathilda.PrattCertifies 1009 c := by mathilda_prime   -- `by decide`
