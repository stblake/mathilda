import Mathilda

/-!
  The oracle is untrusted.  Run with `MATHILDA_BIN=test/bad_oracle.sh`: a fake
  Mathilda that answers every query with a certificate for a *different* (true)
  statement.  Each tactic must fail; none of these false goals is proved.
  `#guard_msgs` pins the rejection, so a clean exit means all three were caught.
-/

/--
error: mathilda_linear_combination: Lean rejected the certificate returned by Mathilda (cofactors [1,
 -y - x]); Mathilda is an untrusted oracle and its answer does not check
-/
#guard_msgs (error) in
example (x y : Rat) (h1 : x^2 + y^2 = 1) (h2 : x = y) : 2*y^2 = 3 := by
  mathilda_linear_combination

/--
error: mathilda_sos: Lean rejected the certificate returned by Mathilda (squares [(1, x - y,
   0)]); Mathilda is an untrusted oracle and its answer does not check
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
