(* ==========================================================================
   Experiment 131 -- ISSAC '97 Systems Challenge, Problem 6
   ==========================================================================

       cd benchmarks/131-issac97-p6-groebner
       ../../Mathilda -file groebner.m

   PROBLEM.  The lexicographic Groebner basis (x > y > z > w) of four dense
   quadrics in four variables.  The basis is in shape position: a degree-16
   univariate polynomial in w and three polynomials linear in x, y, z with
   degree-15 w-coefficients (160 s in 1997).

   METHODS (the notebook's two):
     exact    GroebnerBasis over Q (the exact problem);
     machine  GroebnerBasis[..., CoefficientDomain -> InexactNumbers], the
              notebook's floating-point basis ("150 times faster").
   NOT IN THE NOTEBOOK.  NSolve of the system (the 16 solutions), the
   standard downstream use of a shape basis.
   CHECKS.  The exponent vectors of the basis and the 16-digit leading part
   of the univariate polynomial's constant/leading coefficient ratio.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"GroebnerBasis", "NSolve", "Exponent"}];

eqs = {8 w^2 + 5 w x - 4 w y + 2 w z + 3 w + 5 x^2 + 2 x y - 7 x z - 7 x + 7 y^2 - 8 y z - 7 y + 7 z^2 - 8 z + 8,
       3 w^2 - 5 w x - 3 w y - 6 w z + 9 w + 4 x^2 + 2 x y - 2 x z + 7 x + 9 y^2 + 6 y z + 5 y + 7 z^2 + 7 z + 5,
       -2 w^2 + 9 w x + 9 w y - 7 w z - 4 w + 8 x^2 + 9 x y - 3 x z + 8 x + 6 y^2 - 7 y z + 4 y - 6 z^2 + 8 z + 2,
       7 w^2 + 5 w x + 3 w y - 5 w z - 5 w + 2 x^2 + 9 x y - 7 x z + 4 x - 4 y^2 - 5 y z + 6 y - 4 z^2 - 9 z + 2};
vars = {x, y, z, w};

(* The univariate member is the one free of x, y, z. *)
uni[gb_] := SelectFirst[gb, FreeQ[#, x | y | z] &];
gbChk[gb_, d_] := If[ListQ[gb] && Length[gb] == 4,
  {Sort[(Exponent[#, vars] &) /@ gb],
   sig[Coefficient[uni[gb], w, 0]/Coefficient[uni[gb], w, 16], d]}, "unevaluated"];

run["p6 GroebnerBasis exact lex", GroebnerBasis[eqs, vars], gbChk[#, 16] &, 3];
run["p6 GroebnerBasis InexactNumbers machine",
    GroebnerBasis[eqs, vars, CoefficientDomain -> InexactNumbers], gbChk[#, 6] &, 3];
run["p6 NSolve machine", NSolve[eqs, vars],
    If[ListQ[#], {Length[#], sig[Max[Re[w /. #]], 8]}, "unevaluated"] &, 3];
