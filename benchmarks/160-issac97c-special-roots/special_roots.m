(* bench-timeout: 1200 *)
(* ==========================================================================
   Experiment 160 -- ISSAC '97 holdout C (SEALED), family c5:
   zeros of special functions and special polynomials (counterpart of
   Problem 5)
   ==========================================================================

       cd benchmarks/160-issac97c-special-roots
       ../../Mathilda -file special_roots.m

   Ten zeros: two Lambert W branches (a real one on the -1 branch, a complex
   one on branch 2), the first zero of the Struve function H_0, the DOUBLE
   roots of a terminating 3F2 (by Clausen's formula it is the square of a
   2F1), the first critical-line zero of the Dirichlet L-function of the
   character mod 4, the digamma zero in (-3, -2), the first zero of the
   Kelvin function ber, the parameter where EllipticK = 2 EllipticE, the
   smallest root of the degree-40 partial sum of the exponential series (a
   Szego-curve root, ill-conditioned in a double), and the zero of Zeta' on
   the negative axis.  FindRoot targets are wrapped f[x_?NumericQ]; machine
   rows use default options, wp rows WorkingPrecision -> 40 (FindRoot) or
   WorkingPrecision -> 60 (NSolve on the polynomials).

   TRUE ANSWERS (A = mpmath findroot / polyroots on mpmath's own function,
   60-80 digits; B = PARI/GP at 80-120 digits: lambertw with branch, solve,
   lfunzeros, ellK/ellE, polroots -- or, for Struve and Kelvin, which PARI
   lacks, mpmath findroot on the defining power series summed directly).
   A and B agree to >= 35 digits.  Digits: machine 8 (5 for the double root
   and for the exp partial-sum root, where a double carries 6-7), wp 20 (15
   for the double root at 60 digits).

   c5.1  w E^w = -1/5 on the W_{-1} branch: w = -2.54264135777352642429380615666.
   c5.2  w E^w = 1 on branch 2 (start -2 + 11 I):
         w = -2.40158510486800288417413977468 + 10.7762995161150708984971033464 I.
   c5.3  StruveH[0, x] = 0, first positive zero (start 4.3):
         x = 4.33323782040642167053239927093.
   c5.4  3F2(-30, 31, 1/2; 1, 1; x) = [2F1(-15, 31/2; 1; x)]^2 (Clausen,
         checked coefficient by coefficient in exact arithmetic): every root is
         double.  Smallest: 0.00620338140950646343366027049836.
   c5.5  L(1/2 + I t, chi_-4) = 4^-s (HurwitzZeta[s, 1/4] - HurwitzZeta[s, 3/4])
         = 0, first zero: t = 6.02094890469759665490251152161.
   c5.6  PolyGamma[0, x] = 0 in (-3, -2): x = -2.61072086844414465000153771572.
   c5.7  KelvinBer[0, x] = 0, first zero: x = 2.84891782079513985274970473873.
   c5.8  EllipticK[m] = 2 EllipticE[m] (parameter m):
         m = 0.826114765984970336177373238601.
   c5.9  sum_{k=0}^{40} z^k/k! = 0, root of least modulus with Im > 0:
         z = -12.0759406304122110376023780407 + 0.718482571514124542165114922967 I.
   c5.10 Zeta'[s] = 0 on (-3, -2): s = -2.71726282920457410157058066168.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"FindRoot", "NSolve", "ProductLog", "StruveH", "HurwitzZeta", "PolyGamma", "KelvinBer",
         "EllipticK", "EllipticE", "Zeta"}];

cpx[d_][v_] := If[NumericQ[v], {sig[Re[v], d], sig[Im[v], d]}, "failed"];
fr[d_][r_] := Module[{v = ruleVal[r, t]}, If[NumericQ[v], sig[Re[v], d], "failed"]];
frc[d_][r_] := Module[{v = ruleVal[r, t]}, If[NumericQ[v], cpx[d][v], "failed"]];

g1[t_?NumericQ] := t Exp[t] + 1/5;
run["c5.1 FindRoot machine", FindRoot[g1[t], {t, -3}], fr[8], 3];
run["c5.1 FindRoot wp40", FindRoot[g1[t], {t, -3}, WorkingPrecision -> 40], fr[20]];

g2[t_?NumericQ] := t Exp[t] - 1;
run["c5.2 FindRoot machine", FindRoot[g2[t], {t, -2 + 11 I}], frc[8], 3];
run["c5.2 FindRoot wp40", FindRoot[g2[t], {t, -2 + 11 I}, WorkingPrecision -> 40], frc[20]];

g3[t_?NumericQ] := StruveH[0, t];
run["c5.3 FindRoot machine", FindRoot[g3[t], {t, 43/10}], fr[8], 3];
run["c5.3 FindRoot wp40", FindRoot[g3[t], {t, 43/10}, WorkingPrecision -> 40], fr[20]];

p4 = Sum[Pochhammer[-30, k] Pochhammer[31, k] Pochhammer[1/2, k]/(k!)^3 x^k, {k, 0, 30}];
minAbs[d_][s_] := If[ListQ[s] && Length[s] > 0, sig[Min[Abs[x /. s]], d], "failed"];
run["c5.4 NSolve machine", NSolve[p4 == 0, x], minAbs[5], 3];
run["c5.4 NSolve wp60", NSolve[p4 == 0, x, WorkingPrecision -> 60], minAbs[15]];

g5[t_?NumericQ] := (HurwitzZeta[1/2 + I t, 1/4] - HurwitzZeta[1/2 + I t, 3/4])/4^(1/2 + I t);
run["c5.5 FindRoot machine", FindRoot[g5[t], {t, 6}], fr[8], 3];
run["c5.5 FindRoot wp40", FindRoot[g5[t], {t, 6}, WorkingPrecision -> 40], fr[20]];

g6[t_?NumericQ] := PolyGamma[0, t];
run["c5.6 FindRoot machine", FindRoot[g6[t], {t, -26/10}], fr[8], 3];
run["c5.6 FindRoot wp40", FindRoot[g6[t], {t, -26/10}, WorkingPrecision -> 40], fr[20]];

g7[t_?NumericQ] := KelvinBer[0, t];
run["c5.7 FindRoot machine", FindRoot[g7[t], {t, 28/10}], fr[8], 3];
run["c5.7 FindRoot wp40", FindRoot[g7[t], {t, 28/10}, WorkingPrecision -> 40], fr[20]];

g8[t_?NumericQ] := EllipticK[t] - 2 EllipticE[t];
run["c5.8 FindRoot machine", FindRoot[g8[t], {t, 8/10}], fr[8], 3];
run["c5.8 FindRoot wp40", FindRoot[g8[t], {t, 8/10}, WorkingPrecision -> 40], fr[20]];

p9 = Sum[x^k/k!, {k, 0, 40}];
least[d_][s_] := If[ListQ[s] && Length[s] > 0,
  cpx[d][First[SortBy[Select[x /. s, Im[#] > 0 &], Abs]]], "failed"];
run["c5.9 NSolve machine", NSolve[p9 == 0, x], least[5], 3];
run["c5.9 NSolve wp40", NSolve[p9 == 0, x, WorkingPrecision -> 40], least[20]];

g10[t_?NumericQ] := Zeta'[t];
run["c5.10 FindRoot machine", FindRoot[g10[t], {t, -27/10}], fr[8], 3];
run["c5.10 FindRoot wp40", FindRoot[g10[t], {t, -27/10}, WorkingPrecision -> 40], fr[20]];
