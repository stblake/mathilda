(* bench-timeout: 2400 *)
(* ==========================================================================
   Experiment 159 -- ISSAC '97 holdout C (SEALED), family c4:
   exact coefficients and residues (counterpart of Problem 4)
   ==========================================================================

       cd benchmarks/159-issac97c-coefficients-residues
       ../../Mathilda -file coefficients_residues.m

   Ten coefficients and residues: overpartitions (the eta quotient
   eta(2t)/eta(t)^2), large Schroeder and central Delannoy numbers (algebraic
   generating functions with a square-root singularity at 3 - 2 Sqrt[2]),
   sums of eight squares (theta_3^8), sets of lists (Exp[x/(1-x)], an
   ESSENTIAL singularity at x = 1), involutions (Exp[x + x^2/2]), and the
   residues of Gamma^2 at a double pole, of Zeta^3 at s = 1 (a triple pole:
   needs the Stieltjes constant gamma_1), of 1/(z^15 Sin[z] Sinh[z]) at a
   17th-order pole, and of a rational function at a 30th-order pole.

   THREE ROWS PER PROBLEM.
     exact    SeriesCoefficient / Residue (checked by 15 leading digits and
              the value mod 10^9+7, or the exact rational's numerator and
              denominator mod 10^9+7, or 20 digits of a symbolic constant).
     machine  the Cauchy integral (1/2 Pi) int_0^(2 Pi) f(r e^(it)) (r e^(it))^-n dt
              (or the residue circle) by NIntegrate on a circle of the stated
              radius r, near the saddle point so the integrand's peak exceeds
              the answer by at most ~2e4.  Where the coefficient overflows a
              double (S(1000) ~ 1e761) the machine row is a REDUCED instance
              (n = 150 or 60 or 100), stated per problem.
     wp40     the same integral by Method -> "Trapezoidal" (exponentially
              convergent on a periodic analytic integrand) at 40 digits.

   TRUE ANSWERS.  A: exact integer/rational arithmetic in Python (series
   multiplication, the P-recursive recurrence, a binomial sum, Fractions);
   B: PARI/GP power series (polcoeff of the truncated series), a second
   closed form, or the mpmath trapezoid rule on the stated contour at 30-50
   digits (which also measured the cancellation that sets the digit counts:
   machine 6-7, wp 15, 20 for the residues, 11 for the 11-digit r_8).

   c4.1  pbar(300) = [q^300] prod_k (1+q^k)/(1-q^k) = 175155936072122455560
         (A series product; B PARI).  Contour r = 913/1000 ~ E^(-Pi/(2 Sqrt[300])).
   c4.2  S(n) = [x^n] (1 - x - Sqrt[1 - 6x + x^2])/(2x), large Schroeder.
         A: (n+1) S_n = 3(2n-1) S_(n-1) - (n-2) S_(n-2); B: sum_k C(n,k)
         C(n+k,k)/(k+1) and PARI.  S(1000) = 9.10543459812736067919e760
         (mod 1e9+7: 453374912); machine reduced n = 150:
         2.97634500072195718374e111; contour r = 17/100.
   c4.3  D(n) = [x^n] 1/Sqrt[1 - 6x + x^2], central Delannoy (B: sum_k C(n,k)
         C(n+k,k) and PARI).  exact D(2000) = 1.62223628660764197568e1529
         (mod 476434592); machine n = 150: 3.17858146370208953999533e113;
         wp n = 500: 1.52760689379174916007449381286e381; r = 17/100.
   c4.4  Res_{z=-10} Gamma[z]^2 = 2 PolyGamma[0, 11]/(10!)^2 = 2 (H_10 - gamma)/(10!)^2
         = 3.57186931397414815388235508842e-13; circle |z + 10| = 1/2.
   c4.5  Res_{s=1} Zeta[s]^3 = 3 EulerGamma^2 - 3 StieltjesGamma[1]
         (cube the Laurent series 1/u + gamma - gamma_1 u)
         = 1.21798130787418619753688753669; circle |s - 1| = 1/2.
   c4.6  r_8(1000) = [q^1000] theta_3(q)^8 = 16 sum_{d | n} (-1)^(n+d) d^3
         = 18365675328 (A divisor formula; B series).  Contour r = 99/100,
         theta_3 summed over |k| <= 100.
   c4.7  a(n) = n! [x^n] Exp[x/(1-x)] (A: a(n) = (2n-1) a(n-1) - (n-1)(n-2)
         a(n-2); B: sum_k C(n-1,k-1) n!/k!).  a(300) = 7.99540643814440858750e626
         (mod 27361761); machine reduced n = 60: 3.47714782526715271828e86,
         r = 88/100; wp n = 300, r = 94/100.
   c4.8  Res_{z=0} 1/(z^15 Sin[z] Sinh[z]) = [z^16] z^2/(Sin[z] Sinh[z])
         = 13739/1136785104000 = 1.20858374653719952333225e-8 (A Fraction
         series inversion; B mpmath taylor); circle |z| = 5/2 (poles at Pi, I Pi).
   c4.9  I(n) = n! [x^n] Exp[x + x^2/2], involutions (A: I(n) = I(n-1) +
         (n-1) I(n-2); B: sum_k n!/((n-2k)! k! 2^k)).  I(500) =
         4.23107565308608549950e575 (mod 536729466); machine reduced n = 100:
         2.40533474383334789536e82, r = 10; wp n = 500, r = 22.
   c4.10 Res_{z=1} z^100/((z-1)^30 (z+1)^20) = [u^29] (1+u)^100 (2+u)^-20
         = 7693711204564513690663977952285/35184372088832
         = 218668424297576214.3335026 (A binomial convolution in Fractions;
         B PARI); circle |z - 1| = 1.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"SeriesCoefficient", "Residue", "NIntegrate", "Zeta", "StieltjesGamma"}];

exactChk[v_] := If[IntegerQ[v], {sig[v, 15], Mod[v, 1000000007]}, "unevaluated"];
ratChk[v_] := If[IntegerQ[v] || Head[v] === Rational,
  {sig[v, 15], Mod[Numerator[v], 1000000007], Mod[Denominator[v], 1000000007]}, "unevaluated"];
rsig[d_][v_] := If[!NumericQ[v], "unevaluated",
  Module[{w = If[Precision[v] === Infinity, N[v, d + 10], v]},
    If[Abs[Im[w]] <= 10^-d Abs[w], sig[Re[w], d], "complex"]]];

(* (1/2 Pi) int_0^(2 Pi) h dt, machine and trapezoid at 40 digits; real part. *)
SetAttributes[{cm, cw}, HoldAll];
cm[h_] := Re[NIntegrate[h, {t, 0, 2 Pi}]]/(2 Pi);
cw[h_] := Re[NIntegrate[h, {t, 0, 2 Pi}, Method -> "Trapezoidal", WorkingPrecision -> 40,
                        PrecisionGoal -> 20, MaxRecursion -> 16]]/(2 Pi);

(* integrands: coefficient [x^n] f is f(q) q^-n at q = r e^(it); a residue at
   z0 is f(z) (z - z0) at z = z0 + r e^(it). *)
h1[t_?NumericQ, n_, r_] := Module[{q = r Exp[I t]}, Product[(1 + q^k)/(1 - q^k), {k, 1, n}]/q^n];
h2[t_?NumericQ, n_, r_] := Module[{q = r Exp[I t]}, (1 - q - Sqrt[1 - 6 q + q^2])/(2 q)/q^n];
h3[t_?NumericQ, n_, r_] := Module[{q = r Exp[I t]}, 1/Sqrt[1 - 6 q + q^2]/q^n];
h4[t_?NumericQ, r_] := Module[{w = r Exp[I t]}, Gamma[-10 + w]^2 w];
h5[t_?NumericQ, r_] := Module[{w = r Exp[I t]}, Zeta[1 + w]^3 w];
h6[t_?NumericQ, n_, r_] := Module[{q = r Exp[I t]}, Sum[q^(k^2), {k, -100, 100}]^8/q^n];
h7[t_?NumericQ, n_, r_] := Module[{q = r Exp[I t]}, Exp[q/(1 - q)]/q^n];
h8[t_?NumericQ, r_] := Module[{w = r Exp[I t]}, 1/(w^14 Sin[w] Sinh[w])];
h9[t_?NumericQ, n_, r_] := Module[{q = r Exp[I t]}, Exp[q + q^2/2]/q^n];
h10[t_?NumericQ, r_] := Module[{w = r Exp[I t]}, (1 + w)^100/(w^29 (2 + w)^20)];

run["c4.1 SeriesCoefficient exact",
    SeriesCoefficient[Product[(1 + q^k)/(1 - q^k), {k, 1, 300}], {q, 0, 300}], exactChk];
run["c4.1 contour machine", cm[h1[t, 300, 913/1000]], sig[#, 6] &];
run["c4.1 contour Trapezoidal wp40", cw[h1[t, 300, 913/1000]], sig[#, 15] &, 1, 600];

run["c4.2 SeriesCoefficient exact",
    SeriesCoefficient[(1 - x - Sqrt[1 - 6 x + x^2])/(2 x), {x, 0, 1000}], exactChk];
run["c4.2 contour machine n=150", cm[h2[t, 150, 17/100]], sig[#, 6] &];
run["c4.2 contour Trapezoidal wp40", cw[h2[t, 1000, 17/100]], sig[#, 15] &];

run["c4.3 SeriesCoefficient exact", SeriesCoefficient[1/Sqrt[1 - 6 x + x^2], {x, 0, 2000}], exactChk];
run["c4.3 contour machine n=150", cm[h3[t, 150, 17/100]], sig[#, 6] &];
run["c4.3 contour Trapezoidal wp40 n=500", cw[h3[t, 500, 17/100]], sig[#, 15] &];

run["c4.4 Residue exact", Residue[Gamma[z]^2, {z, -10}], rsig[20]];
run["c4.4 contour machine", cm[h4[t, 1/2]], sig[#, 7] &];
run["c4.4 contour Trapezoidal wp40", cw[h4[t, 1/2]], sig[#, 20] &];

run["c4.5 Residue exact", Residue[Zeta[s]^3, {s, 1}], rsig[20]];
run["c4.5 contour machine", cm[h5[t, 1/2]], sig[#, 7] &];
run["c4.5 contour Trapezoidal wp40", cw[h5[t, 1/2]], sig[#, 20] &];

run["c4.6 SeriesCoefficient exact",
    SeriesCoefficient[Sum[q^(k^2), {k, -40, 40}]^8, {q, 0, 1000}], exactChk];
run["c4.6 contour machine", cm[h6[t, 1000, 99/100]], sig[#, 6] &];
run["c4.6 contour Trapezoidal wp40", cw[h6[t, 1000, 99/100]], sig[#, 11] &, 1, 600];

run["c4.7 SeriesCoefficient exact", 300! SeriesCoefficient[Exp[x/(1 - x)], {x, 0, 300}], exactChk];
run["c4.7 contour machine n=60", 60! cm[h7[t, 60, 88/100]], sig[#, 6] &];
run["c4.7 contour Trapezoidal wp40", 300! cw[h7[t, 300, 94/100]], sig[#, 15] &];

run["c4.8 Residue exact", Residue[1/(z^15 Sin[z] Sinh[z]), {z, 0}], ratChk];
run["c4.8 contour machine", cm[h8[t, 5/2]], sig[#, 7] &];
run["c4.8 contour Trapezoidal wp40", cw[h8[t, 5/2]], sig[#, 20] &];

run["c4.9 SeriesCoefficient exact", 500! SeriesCoefficient[Exp[x + x^2/2], {x, 0, 500}], exactChk];
run["c4.9 contour machine n=100", 100! cm[h9[t, 100, 10]], sig[#, 6] &];
run["c4.9 contour Trapezoidal wp40", 500! cw[h9[t, 500, 22]], sig[#, 15] &, 1, 600];

run["c4.10 Residue exact", Residue[z^100/((z - 1)^30 (z + 1)^20), {z, 1}], ratChk];
run["c4.10 contour machine", cm[h10[t, 1]], sig[#, 6] &];
run["c4.10 contour Trapezoidal wp40", cw[h10[t, 1]], sig[#, 20] &];
