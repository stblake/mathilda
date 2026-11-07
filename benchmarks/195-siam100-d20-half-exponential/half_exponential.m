(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 195 -- SIAM 100-Digit Challenge, Appendix D, Problem 20
   ("half-exponential").

     h(z) is the solution of exp(z) = h(h(z)) that is analytic near the
     origin and increasing for real z.  What is h(0)?  (And the radius of
     convergence of its Maclaurin series -- not asked here.)

   A RESEARCH PROBLEM as posed: the book's remark says further conditions are
   needed for uniqueness and none was known.  The commonly used solution is
   Kneser's (1950) real-analytic half-iterate h = sexp(slog(z) + 1/2), for which
   a uniqueness criterion was later given by Paulsen & Cowgill (Adv. Comput.
   Math. 2017).  slog(0) = -1, so h(0) = sexp(-1/2), the x with
   slog(x) = -1/2.  Only 7 digits are established here (truth.tsv):
   h(0) = 0.4985633.  Robbins' truncated Abel system settles at
   0.49856328509 (n = 140, 160), while the commonly quoted Kneser value is
   0.4985632879: they differ in the 9th digit, unresolved here.

   No system has a built-in tetration / Kneser half-iterate, so the problem
   rows d20.machine and d20.wp50 are SKIP everywhere (N/A).

   THE PROXY rows ask a well-defined finite version: Robbins' "intuitive
   slog".  slog_n(x) = -1 + Sum_{k=1}^n v_k x^k, with v solving the n x n
   linear system that makes slog_n(e^x) - slog_n(x) = 1 hold through x^(n-1):
        Sum_k v_k (k^j/j! - [j == k]) = [j == 0],   j = 0..n-1,
   and h_n(0) = the root of slog_n(x) = -1/2 near 1/2.  h_n(0) -> h(0) slowly
   (n = 60: 1e-7 off).  Idiomatic statement: LinearSolve + FindRoot.
     d20.machine robbins n=20   machine precision (cond(A) ~ 1e8: a perfect
                                double computation keeps ~15 digits)
     d20.wp50 robbins n=60      exact LinearSolve, FindRoot WorkingPrecision 60
   Python: numpy.linalg.solve + scipy brentq / mpmath lu_solve + findroot.
   Maxima: linsolve (exact) + find_root / bf_find_root.  FriCAS: exact
   solve of the linear system over Fraction Integer, then the built-in
   real-root isolation solve(p = 0, eps) (rational bracket to eps), converted
   to DoubleFloat / Float.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"LinearSolve", "FindRoot"}];

Print["SKIP\td20.machine\tKneserHalfIterate"];
Print["SKIP\td20.wp50\tKneserHalfIterate"];

robbinsH0[n_, opts___] := Module[{A, v},
  A = Table[k^j/j! - KroneckerDelta[j, k], {j, 0, n - 1}, {k, 1, n}];
  v = LinearSolve[If[{opts} === {}, N[A], A], UnitVector[n, 1]];
  x /. FindRoot[-1 + v . x^Range[n] == -1/2, {x, 1/2}, opts]];

vrun["d20.machine robbins n=20", robbinsH0[20], 10, 3];
vrun["d20.wp50 robbins n=60", robbinsH0[60, WorkingPrecision -> 60], 48];
