import Mathilda

/-!
  Run with `MATHILDA_BIN=/nonexistent/Mathilda`: a missing binary is a clear
  error, never a proof.
-/

/--
error: mathilda: cannot run the Mathilda executable '/nonexistent/Mathilda' (no such file). Set the environment variable MATHILDA_BIN, or `set_option mathilda.path "/path/to/Mathilda"`, or put Mathilda on PATH.
-/
#guard_msgs (error) in
example (x : Rat) (h : x = 1) : x^2 = 1 := by
  mathilda_linear_combination
