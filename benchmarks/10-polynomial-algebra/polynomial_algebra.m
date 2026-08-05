(* ==========================================================================
   Experiment 10 -- Polynomial algebra: GCD, resultants, canonical forms
   ==========================================================================
   WHAT IT MEASURES.  The plumbing every other symbolic subsystem sits on:
   src/poly/subresultants.c, ratcanon.c, src/rat.c, src/parfrac.c, src/expand.c.

   WHY IT IS A SEPARATE EXPERIMENT.  Integrate, Simplify, Solve and Factor all
   call these.  When one of those rows is slow, this experiment answers whether
   the cost is in the high-level algorithm or in the arithmetic underneath it --
   the "fix the primitive, then re-profile the composition" lesson that the
   existing corpus records hitting six separate times.

   Checks are DEGREES and evaluated VALUES: both ordering-independent.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"PolynomialGCD", "PolynomialLCM", "Resultant", "Discriminant",
         "Cancel", "Exponent", "CoefficientList", "PolynomialQuotient",
         "PolynomialRemainder", "Subresultants"}];

(* 1. GCD of two dense degree-60 polynomials sharing a degree-20 factor. *)
sh = factorableUPoly[10, x];
a1 = Expand[sh denseUPoly[40, x]];
b1 = Expand[sh denseUPoly[38, x]];
bench["PolynomialGCD, shared deg-20 factor", PolynomialGCD[a1, b1];];
check["PolynomialGCD, shared deg-20 factor", Exponent[PolynomialGCD[a1, b1], x]];

(* 2. Coprime inputs: the worst case, since the remainder sequence runs to the
   bottom before proving the GCD is 1. *)
bench["PolynomialGCD, coprime deg 40",
  PolynomialGCD[denseUPoly[40, x], denseUPoly[39, x] + 1];];
check["PolynomialGCD, coprime deg 40",
  Exponent[PolynomialGCD[denseUPoly[40, x], denseUPoly[39, x] + 1], x]];

(* 3. Discriminant of a degree-20 polynomial: a resultant with its derivative. *)
bench["Discriminant of deg 20", Discriminant[denseUPoly[20, x], x];];
check["Discriminant of deg 20",
  Mod[Discriminant[denseUPoly[20, x], x], 1000003]];

(* 4. Cancel on a large ratio with a shared factor. *)
bench["Cancel deg-60 over deg-58", Cancel[a1/b1];];
check["Cancel deg-60 over deg-58", Round[10^6 N[Cancel[a1/b1] /. x -> 11/10]]];

(* 5. Division with remainder. *)
bench["PolynomialQuotient deg 60 / deg 20",
  PolynomialQuotient[a1, sh, x];];
check["PolynomialQuotient deg 60 / deg 20",
  Exponent[PolynomialQuotient[a1, sh, x], x]];

(* 6. Expand a large power: the raw term-generation cost. *)
bench["Expand (1+x)^400", Expand[(1 + x)^400];];
check["Expand (1+x)^400", Coefficient[Expand[(1 + x)^400], x, 200]];
