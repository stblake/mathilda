(* ==========================================================================
   Experiment 07 -- Groebner bases
   ==========================================================================
   WHAT IT MEASURES.  src/poly/groebner.c, groebnerbasis.c, groebnerwalk.c:
   Buchberger with pair selection, plus the walk between orderings.

   cyclic-n is the standard stress family: cyclic-5 is easy, cyclic-6 is where
   naive pair selection starts to lose, cyclic-7 is a genuine test.  The suite
   stops at 6 deliberately -- a row that cannot finish inside the timeout tells
   us less than a row that finishes and can be compared.

   Checks are BASIS LENGTHS.  A reduced Groebner basis for a fixed ordering is
   unique, so its cardinality is a real invariant here (unlike Factor's ordering).
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"GroebnerBasis", "MonomialOrder", "Eliminate", "Resultant"}];

c3 = cyclicSystem[{x, y, z}];
bench["GroebnerBasis cyclic-3", GroebnerBasis[c3, {x, y, z}];];
check["GroebnerBasis cyclic-3", Length[GroebnerBasis[c3, {x, y, z}]]];

c4 = cyclicSystem[{x, y, z, u}];
bench["GroebnerBasis cyclic-4", GroebnerBasis[c4, {x, y, z, u}];];
check["GroebnerBasis cyclic-4", Length[GroebnerBasis[c4, {x, y, z, u}]]];

c5 = cyclicSystem[{x, y, z, u, v}];
bench["GroebnerBasis cyclic-5", GroebnerBasis[c5, {x, y, z, u, v}], 1];
check["GroebnerBasis cyclic-5", Length[GroebnerBasis[c5, {x, y, z, u, v}]]];

(* A small dense system in 3 variables, unrelated to the cyclic family. *)
s1 = {x^2 + y + z - 1, x + y^2 + z - 1, x + y + z^2 - 1};
bench["GroebnerBasis 3-var quadratic system", GroebnerBasis[s1, {x, y, z}];];
check["GroebnerBasis 3-var quadratic system",
  Length[GroebnerBasis[s1, {x, y, z}]]];

(* Elimination: project a system onto one variable. *)
bench["Eliminate 2 vars from 3", Eliminate[s1, {y, z}];];
check["Eliminate 2 vars from 3",
  Exponent[Total[Flatten[{Eliminate[s1, {y, z}]}] /. Equal -> Subtract], x]];

(* Resultant: the two-polynomial special case of elimination. *)
bench["Resultant of two deg-8 polys",
  Resultant[denseUPoly[8, x] + y, denseUPoly[6, x] - y, x];];
check["Resultant of two deg-8 polys",
  Exponent[Resultant[denseUPoly[8, x] + y, denseUPoly[6, x] - y, x], y]];
