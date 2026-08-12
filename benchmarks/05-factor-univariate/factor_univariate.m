(* ==========================================================================
   Experiment 05 -- Univariate polynomial factorisation
   ==========================================================================
   WHAT IT MEASURES.  src/poly/facpoly.c and zupoly.c: Zassenhaus / Cantor-
   Zassenhaus over Z, and the FLINT bridge (src/poly/flint_bridge.c) when
   USE_FLINT is compiled in.

   READ THE BUILD WARNING IN REPORT.md BEFORE ACTING ON THESE ROWS.  Without
   FLINT >= 3.0 the makefile silently substitutes "the classical fallback", so
   these rows measure a different algorithm than an FLINT build does.  That is
   not a Mathilda defect; it is a build configuration, and the report says which
   one it measured.

   VALUE CHECKS ARE FACTOR COUNTS.  Factor's output ordering is not canonical
   across systems, but the NUMBER of irreducible factors (with multiplicity) is.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"Factor", "FactorList", "Expand", "IrreduciblePolynomialQ",
         "PolynomialGCD", "Cyclotomic"}];

(* 1. A product of 8 distinct quadratics: the answer is known by construction. *)
f1 = Expand[factorableUPoly[8, x]];
bench["Factor product of 8 quadratics", Factor[f1];];
check["Factor product of 8 quadratics", Length[FactorList[f1]]];

(* 2. The same shape at 16 factors -- degree 32. *)
f2 = Expand[factorableUPoly[16, x]];
bench["Factor product of 16 quadratics", Factor[f2];];
check["Factor product of 16 quadratics", Length[FactorList[f2]]];

(* 3. x^n - 1: factors into cyclotomics.  Many small factors, high degree. *)
f3 = x^120 - 1;
bench["Factor x^120 - 1", Factor[f3];];
check["Factor x^120 - 1", Length[FactorList[f3]]];

(* 4. A sparse high-degree polynomial: cheap to store, expensive to factor. *)
f4 = sparseUPoly[60, x];
bench["Factor sparse degree 60", Factor[f4];];
check["Factor sparse degree 60", Length[FactorList[f4]]];

(* 5. A dense degree-60 polynomial with mild coefficients, probably irreducible:
   the WORST case, because the algorithm must exhaust its search to prove it. *)
f5 = denseUPoly[60, x];
bench["Factor dense degree 60", Factor[f5];];
check["Factor dense degree 60", Length[FactorList[f5]]];

(* 6. GCD of two large polynomials sharing a known factor. *)
g1 = Expand[factorableUPoly[10, x]];
g2 = Expand[factorableUPoly[6, x] (x^4 + x + 1)];
bench["PolynomialGCD of two deg-20 polys", PolynomialGCD[g1, g2];];
check["PolynomialGCD of two deg-20 polys", Exponent[PolynomialGCD[g1, g2], x]];
