(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 140 -- ISSAC '97 holdout, family 5: roots of special functions
   ==========================================================================

       cd benchmarks/140-issac97h-special-roots
       ../../Mathilda -file special_roots.m

   HOLDOUT.  Ten problems probing what experiment 130 (the largest zero of
   LaguerreL[1000, x], where |L| ~ 1e867) probes -- extreme zeros of
   high-degree orthogonal polynomials and zeros of special functions, some
   where the function's values leave the double range -- on different
   families.  See benchmarks/ISSAC97_HOLDOUT.md.

   Each problem: FindRoot at machine precision (checked to 9 digits; where
   the full problem's values exceed 1e308 near the root the machine case is
   a reduced degree, stated below) and FindRoot with WorkingPrecision -> 50
   (checked to 20 digits, below the default goal of WP/2), both from the
   same starting point (an asymptotic estimate rounded to a few digits).
   Zeros next to 1 are checked as 1 - x, so the check sees the digits that
   matter (machine: 6 digits of 1 - x, which a double resolves from x).

   PROBLEMS AND TRUE VALUES (two methods each, agreeing to all digits shown).
   Orthogonal polynomials: a 60-digit Newton iteration on the three-term
   recurrence (mpmath), and a Sturm-sequence count -- the number of sign
   changes in p_0(x), ..., p_n(x) is the number of zeros of p_n above x --
   evaluated at r -/+ 10^-30, which brackets the claimed zero and proves it
   is the extreme one (for h5.2 and h5.3 also a 150-step bisection on that
   count).  Special functions: a library zero routine and/or a published
   table, plus the same Newton iteration:
   h5.1  largest zero of HermiteH[500, x] = 31.05074638009002022024205371
         (|H_500| ~ 1e1100 there; machine: HermiteH[100, x], 13.40648733814491013849801536)
   h5.2  1 - (largest zero of LegendreP[1000, x]) = 2.888701924489430123709748e-6
   h5.3  1 - (largest zero of JacobiP[300, 3, 1/2, x]) = 2.227861165971667686640295e-4
   h5.4  1 - (largest zero of ChebyshevT[2001, x]) = 1 - Cos[Pi/4002]
         = 3.081169277385704946343859e-7                            (closed form)
   h5.5  100th zero of BesselJ[0, x] = 313.3742660775278447196902451  (library, McMahon)
   h5.6  first zero of BesselJ[100, x] = 108.8361658984097743630979920 (library)
   h5.7  largest zero of LaguerreL[800, 10, x] = 3167.337473886621268542856
         (machine: LaguerreL[80, 10, x], 316.1495015360838567219625)
   h5.8  smallest zero of LaguerreL[2000, x] = 7.227175802102313269819995e-4
   h5.9  50th zero of AiryAi = -38.02100867725525443313246829   (library zero routine)
   h5.10 10th zero of the Riemann zeta function on the critical line, as the
         zero of RiemannSiegelZ[t] = 49.77383247767230218191678468 (published,
         Odlyzko's tables; and a 50-digit zeta-zero routine)
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"FindRoot", "HermiteH", "LegendreP", "JacobiP", "ChebyshevT", "BesselJ",
         "LaguerreL", "AiryAi", "RiemannSiegelZ"}];

(* Each function is wrapped so that it is only ever evaluated at a NUMBER:
   with a symbolic x, LegendreP[1000, x] or ChebyshevT[2001, x] expands into
   its power form, whose coefficients (up to ~1e600) cancel catastrophically
   at any working precision a root finder would choose.  The wrappers keep the
   evaluation on the recurrence / special-function path, as experiment 130
   does with f[x_?NumberQ].  FindRoot then uses a numerical derivative. *)
hF[n_][x_?NumericQ] := HermiteH[n, x];
pF[n_][x_?NumericQ] := LegendreP[n, x];
jF[n_, a_, b_][x_?NumericQ] := JacobiP[n, a, b, x];
tF[n_][x_?NumericQ] := ChebyshevT[n, x];
bF[nu_][x_?NumericQ] := BesselJ[nu, x];
lF[n_, a_][x_?NumericQ] := LaguerreL[n, a, x];
aiF[x_?NumericQ] := AiryAi[x];
zF[x_?NumericQ] := RiemannSiegelZ[x];

fm[f_, x0_] := ruleVal[FindRoot[f[x], {x, x0}], x];
fw[f_, x0_] := ruleVal[FindRoot[f[x], {x, x0}, WorkingPrecision -> 50, MaxIterations -> 200], x];
one[v_] := If[NumericQ[v], 1 - v, "failed"];

run["h5.1 FindRoot machine H100", fm[hF[100], 134/10], sig[#, 9] &, 3];
run["h5.1 FindRoot wp50 H500",    fw[hF[500], 3105/100], sig[#, 20] &];
run["h5.2 FindRoot machine",      one[fm[pF[1000], 999997/1000000]], sig[#, 6] &, 3];
run["h5.2 FindRoot wp50",         one[fw[pF[1000], 999997/1000000]], sig[#, 20] &];
run["h5.3 FindRoot machine",      one[fm[jF[300, 3, 1/2], 9998/10000]], sig[#, 6] &, 3];
run["h5.3 FindRoot wp50",         one[fw[jF[300, 3, 1/2], 9998/10000]], sig[#, 20] &];
run["h5.4 FindRoot machine",      one[fm[tF[2001], 9999997/10000000]], sig[#, 6] &, 3];
run["h5.4 FindRoot wp50",         one[fw[tF[2001], 9999997/10000000]], sig[#, 20] &];
run["h5.5 FindRoot machine",      fm[bF[0], 3134/10], sig[#, 9] &, 3];
run["h5.5 FindRoot wp50",         fw[bF[0], 3134/10], sig[#, 20] &];
run["h5.6 FindRoot machine",      fm[bF[100], 1088/10], sig[#, 9] &, 3];
run["h5.6 FindRoot wp50",         fw[bF[100], 1088/10], sig[#, 20] &];
run["h5.7 FindRoot machine L80",  fm[lF[80, 10], 3161/10], sig[#, 9] &, 3];
run["h5.7 FindRoot wp50 L800",    fw[lF[800, 10], 31673/10], sig[#, 20] &];
run["h5.8 FindRoot machine",      fm[lF[2000, 0], 723/1000000], sig[#, 9] &, 3];
run["h5.8 FindRoot wp50",         fw[lF[2000, 0], 723/1000000], sig[#, 20] &];
run["h5.9 FindRoot machine",      fm[aiF, -3802/100], sig[#, 9] &, 3];
run["h5.9 FindRoot wp50",         fw[aiF, -3802/100], sig[#, 20] &];
run["h5.10 FindRoot machine",     fm[zF, 4977/100], sig[#, 9] &, 3];
run["h5.10 FindRoot wp50",        fw[zF, 4977/100], sig[#, 20] &];
