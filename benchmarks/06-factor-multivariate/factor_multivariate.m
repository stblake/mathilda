(* ==========================================================================
   Experiment 06 -- Multivariate factorisation and expansion
   ==========================================================================
   WHAT IT MEASURES.  src/poly/mvfactor.c, mvfactor3.c and mpoly.c: Hensel
   lifting in several variables, plus the Expand path that feeds it.

   Also gated on USE_FLINT -- see the build-warning block in REPORT.md.

   Checks are TERM COUNTS and FACTOR COUNTS, both ordering-independent.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"Factor", "FactorList", "Expand", "Exponent", "Coefficient"}];

(* 1. Expand a product of 6 linear forms in 3 variables. *)
bench["Expand (x+y+z+i), 6 factors", multiPoly[6, {x, y, z}];];
check["Expand (x+y+z+i), 6 factors", Length[multiPoly[6, {x, y, z}]]];

(* 2. Factor it back: the round trip is the real test. *)
m6 = multiPoly[6, {x, y, z}];
bench["Factor 3-var, 6 factors", Factor[m6];];
check["Factor 3-var, 6 factors", Length[FactorList[m6]]];

(* 3. Five variables, four factors: the variable count is what hurts. *)
m5 = multiPoly[4, {x, y, z, u, v}];
bench["Factor 5-var, 4 factors", Factor[m5];];
check["Factor 5-var, 4 factors", Length[FactorList[m5]]];

(* 4. A symmetric form that does NOT factor: the prove-irreducible cost. *)
irr = Expand[x^3 + y^3 + z^3 - 3 x y z + 1];
bench["Factor irreducible 3-var cubic", Factor[irr];];
check["Factor irreducible 3-var cubic", Length[FactorList[irr]]];

(* 5. A power of a sum: (x+y+z)^12 expands to a large multinomial. *)
bench["Expand (x+y+z)^12", Expand[(x + y + z)^12];];
check["Expand (x+y+z)^12", Length[Expand[(x + y + z)^12]]];

(* 6. Factor a difference of high powers in two variables. *)
d2 = x^24 - y^24;
bench["Factor x^24 - y^24", Factor[d2];];
check["Factor x^24 - y^24", Length[FactorList[d2]]];
