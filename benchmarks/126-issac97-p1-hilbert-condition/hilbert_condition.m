(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 126 -- ISSAC '97 Systems Challenge, Problem 1
   ==========================================================================

       cd benchmarks/126-issac97-p1-hilbert-condition
       ../../Mathilda -file hilbert_condition.m

   PROBLEM.  The 2-norm condition number of the 256 x 256 Hilbert matrix,
   kappa = lambda_max / lambda_min ~ 1.7659543288e389.  No double can hold the
   answer (and lambda_min ~ 1e-389), so the problem is arbitrary precision by
   construction: ~3.5 n digits are lost, the notebook carries 800.

   METHODS (all three of the 1997 solution notebook, verbatim but for Random ->
   a deterministic start vector):
     m1  Eigenvalues[N[Hilbert[256], 800]], ratio of the extreme eigenvalues;
     m2  power iteration on H and on the closed-form inverse Hilbert matrix
         (the Rayleigh quotient of each, iterated to a fixed point);
     m3  LUDecomposition at 500 digits, then 8 power steps on H and 8 inverse
         steps through LUBackSubstitution (`require`d: Mathilda lacks it), and
         the same estimator through LinearSolve[h] (a LinearSolveFunction).
   MACHINE VARIANT.  The same three methods at n = 10 (kappa = 1.6e13), the
   largest size a double still resolves to ~3 digits.
   ARBITRARY-PRECISION SCALING.  m1 at n = 64 / 200 digits and n = 256 / 800.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Eigenvalues", "LUDecomposition", "LUBackSubstitution", "LinearSolve", "FixedPoint"}];
(* Captured BEFORE any definition below mentions the symbol: parsing a
   definition creates the symbol, after which Names[] would report it. *)
haveLUB = haveQ["LUBackSubstitution"];

Hilbert[n_] := Table[1/(i + j - 1), {i, n}, {j, n}];
InverseHilbert[n_] := Table[(-1)^(i + j) (i + j - 1) Binomial[n + i - 1, n - j] *
    Binomial[n + j - 1, n - i] Binomial[i + j - 2, i - 1]^2, {i, n}, {j, n}];

(* ---- m1: all eigenvalues ------------------------------------------------ *)
m1[n_, prec_] := Module[{e = Eigenvalues[N[Hilbert[n], prec]]}, Max[e]/Min[e]];

run["p1 m1 eigenvalues machine n=10", m1[10, MachinePrecision], sig[#, 3] &, 3];
run["p1 m1 eigenvalues wp200 n=64", m1[64, 200], sig[#, 10] &];
run["p1 m1 eigenvalues wp800 n=256", m1[256, 800], sig[#, 10] &, 1, 900];

(* ---- m2: power method on H and its exact inverse ------------------------ *)
normalize[v_List] := v/Sqrt[v . v];
m2[n_, prec_] := Module[{h, ih, rq, fp},
  {h, ih} = N[{Hilbert[n], InverseHilbert[n]}, prec];
  rq[u_List, v_List] := (u . h . u v . ih . v)/(u . u v . v);
  fp = FixedPoint[Function[x, ({{#1, #2}, rq[#1, #2]} &) @@
                     {normalize[h . x[[1, 1]]], normalize[ih . x[[1, 2]]]}],
         {Table[1, {2}, {n}], 2},
         SameTest -> (Abs[1 - #1[[2]]/#2[[2]]] < 10^-6 &)];
  fp[[2]]];

run["p1 m2 power machine n=10", m2[10, MachinePrecision], sig[#, 3] &, 3];
run["p1 m2 power wp800 n=256", m2[256, 800], sig[#, 5] &, 1, 900];

(* ---- m3: LU, then power / inverse-power estimates of the two norms ------ *)
(* Deterministic start vectors replace the notebook's Random[Real, {1, 2}]:
   the estimate's 3rd digit depends on the start, and the check must not. *)
estimateNorm2[mat_?MatrixQ] := Module[{l = Length[mat], x, y},
  y = N[Table[1 + k/l, {k, l}], Precision[mat[[1, 1]]]];
  Do[x = y; y = mat . x, {8}];
  Sqrt[y . Conjugate[y]/x . Conjugate[x]]];
estimateInvNorm2[solve_, l_, prec_] := Module[{x, y},
  x = N[Table[1 + k/l, {k, l}], prec];
  Do[y = x; x = solve[y], {8}];
  Sqrt[x . Conjugate[x]/y . Conjugate[y]]];

m3lu[n_, prec_] := Module[{h = N[Hilbert[n], prec], lud},
  lud = LUDecomposition[h];
  estimateNorm2[h] estimateInvNorm2[LUBackSubstitution[lud, #] &, n, prec]];
m3ls[n_, prec_] := Module[{h = N[Hilbert[n], prec], f},
  f = LinearSolve[h];
  estimateNorm2[h] estimateInvNorm2[f, n, prec]];

run["p1 m3 LUDecomposition only wp500 n=256",
    Length[LUDecomposition[N[Hilbert[256], 500]]], # &];
If[haveLUB,
  run["p1 m3 LUBackSubstitution wp500 n=256", m3lu[256, 500], sig[#, 2] &],
  Print["SKIP\tp1 m3 LUBackSubstitution wp500 n=256\tLUBackSubstitution"]];
run["p1 m3 LinearSolve machine n=10", m3ls[10, MachinePrecision], sig[#, 1] &, 3];
run["p1 m3 LinearSolve wp500 n=256", m3ls[256, 500], sig[#, 2] &];
