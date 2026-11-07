(* Reference B for the finite proxy rows of SIAM100 App. D problem 20 (exp 195):
   Robbins' truncated Abel system, solved in EXACT rational arithmetic by
   LinearSolve, then the root of slog_n(x) = -1/2 by FindRoot at 100 digits.
   (Reference A: siam100_refs/d20_ref.py, own fraction Gaussian elimination +
   mpmath.findroot.)  Also prints h_n(0) for larger n (the sequence that
   approaches the Kneser value).  Run: wolframscript -file d20_ref.wl *)
robbins[n_] := Module[{A, b, v},
  A = Table[k^j/j! - KroneckerDelta[j, k], {j, 0, n - 1}, {k, 1, n}];
  b = UnitVector[n, 1];
  v = LinearSolve[A, b];
  x /. FindRoot[-1 + v . x^Range[n] == -1/2, {x, 1/2}, WorkingPrecision -> 100,
     PrecisionGoal -> 90]];
Do[Print[n, " ", N[robbins[n], 70]], {n, {20, 60}}];
Do[Print[n, " ", N[robbins[n], 20]], {n, {80, 100, 120, 140, 160}}];
