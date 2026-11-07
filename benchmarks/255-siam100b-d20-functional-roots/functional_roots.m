(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 255 -- SIAM 100-Digit Challenge, SEALED HOLDOUT B, family D20
   (Appendix D Problem 20, "half-exponential": h(h(z)) = exp(z)).

   Ten new finite linear-algebra proxies of functional square roots and Abel
   functions, the capability exp 195 tests: build an exactly specified
   (ill-conditioned) matrix, solve it / take a matrix function, then read off
   a value, possibly by a polynomial root.  Truths: truth.tsv
   (siam100_refs/bd20_ref.py == siam100_refs/bd20_ref.gp).

   Robbins' truncated Abel system for base b, L = Log[b]:
       slog_n(x) = -1 + Sum_{k=1}^n v_k x^k,
       Sum_k v_k ((k L)^j/j! - [j == k]) = [j == 0],  j = 0..n-1,
   and the root x of slog_n(x) = level near x0.
     bd20.1   b = E,      n = 10, level -1/2   (h_n(0), easy)
     bd20.2   b = 2,      n = 20, level -1/2
     bd20.3   b = Sqrt2,  n = 30, level -1/2   (b < e^(1/e): real fixed point 2)
     bd20.4   b = 10,     n = 15, level -1/2   (badly scaled: entries to 1e14)
     bd20.5   b = E,      n = 30, level -3/4   (a quarter-iterate value)
     bd20.6   b = E,      n = 24, v_1 = slog_n'(0)  (LinearSolve alone)
     bd20.10  b = E,      n = 45, level -1/2   (harder; a double carries ~8
                                                digits, need 7 on the machine row)
   Carleman square roots: C_jk = [x^k] g(x)^j (j, k = 1..n, upper triangular);
   the first row of MatrixPower[C, 1/2] is the degree-n Taylor polynomial of the
   half-iterate h (h(h(x)) = g(x)), evaluated at x0:
     bd20.7   g = 2x + x^2,       n = 20, x0 = 1/2  (h = (1+x)^Sqrt2 - 1)
     bd20.8   g = Exp[x] - 1,     n = 15, x0 = 1/2  (unipotent: parabolic
                                                     fixed point, not diagonalisable)
     bd20.9   g = 2x/(1 - x^2),   n = 25, x0 = 3/10 (h = Tan[Sqrt2 ArcTan[x]])
   Rows: "bd20.k machine" (machine numbers: N of the matrix) and
   "bd20.k wp50" (base E: the exact rational matrix, as in exp 195; otherwise
   N[matrix, 60]; FindRoot / MatrixPower at WorkingPrecision 60).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"LinearSolve", "FindRoot", "MatrixPower", "StirlingS2"}];

robbinsA[n_, L_] := Table[(k L)^j/j! - KroneckerDelta[j, k], {j, 0, n - 1}, {k, 1, n}];

(* wp: MachinePrecision, or 60 (exact matrix when L is rational) *)
robbinsV[n_, L_, wp_] := Module[{A = robbinsA[n, L]},
  LinearSolve[Which[wp === MachinePrecision, N[A], Element[L, Rationals], A, True, N[A, wp]],
    UnitVector[n, 1]]];

robbinsRoot[n_, L_, level_, x0_, wp_] := Module[{v = robbinsV[n, L, wp], x},
  x /. FindRoot[-1 + v . x^Range[n] == level, {x, x0}, WorkingPrecision -> wp]];

carleman[coef_, n_] := Table[coef[j, k], {j, 1, n}, {k, 1, n}];
c7[j_, k_] := If[j <= k <= 2 j, Binomial[j, k - j] 2^(2 j - k), 0];
c8[j_, k_] := If[k >= j, j! StirlingS2[k, j]/k!, 0];
c9[j_, k_] := If[k >= j && EvenQ[k - j], 2^j Binomial[j + (k - j)/2 - 1, (k - j)/2], 0];

halfIterate[coef_, n_, x0_, wp_] := Module[{C = carleman[coef, n], S},
  S = MatrixPower[If[wp === MachinePrecision, N[C], N[C, wp]], 1/2];
  First[S] . x0^Range[n]];

Do[With[{lab = p[[1]], f = p[[2]], dm = p[[3]]},
    vrun[lab <> " machine", f[MachinePrecision], dm, 3, 60];
    vrun[lab <> " wp50", f[60], 48, 1, 90]],
  {p, {
    {"bd20.1", robbinsRoot[10, 1, -1/2, 1/2, #] &, 10},
    {"bd20.2", robbinsRoot[20, Log[2], -1/2, 54/100, #] &, 10},
    {"bd20.3", robbinsRoot[30, Log[2]/2, -1/2, 63/100, #] &, 10},
    {"bd20.4", robbinsRoot[15, Log[10], -1/2, 39/100, #] &, 10},
    {"bd20.5", robbinsRoot[30, 1, -3/4, 27/100, #] &, 10},
    {"bd20.6", First[robbinsV[24, 1, #]] &, 10},
    {"bd20.7", halfIterate[c7, 20, 1/2, #] &, 10},
    {"bd20.8", halfIterate[c8, 15, 1/2, #] &, 10},
    {"bd20.9", halfIterate[c9, 25, 3/10, #] &, 10},
    {"bd20.10", robbinsRoot[45, 1, -1/2, 1/2, #] &, 7}}}];
