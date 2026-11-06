(* bench-timeout: 1200 *)
(* ==========================================================================
   Experiment 129 -- ISSAC '97 Systems Challenge, Problem 4
   ==========================================================================

       cd benchmarks/129-issac97-p4-polynomial-coefficient
       ../../Mathilda -file polynomial_coefficient.m

   PROBLEM.  The coefficient of x^3000 in
       p(x) = (x + 1)^2000 (x^2 + x + 1)^1000 (x^4 + x^3 + x^2 + x + 1)^500,
   an exact 1427-digit integer, 3.97394226558003...e1426.  p has degree 6000,
   so x^3000 is the middle coefficient.

   METHODS (all four of the notebook's):
     m1  the residue of p/x^3001 at 0 as a contour integral over the unit
         circle: NIntegrate at machine precision (whose integrand peaks at
         |p| = 3^1000 5^500 2^2000 ~ 1e1430 -- outside the double range), and
         verbatim with Method -> "Trapezoidal", WorkingPrecision -> 400
         (exact for a periodic integrand once the rule has > 6000 points);
     m2  Residue[p/x^3001, {x, 0}];
     m3  Series[p, {x, 0, 3000}] then SeriesCoefficient;
     m4  p = q^2 with q = (x+1)^1000 (x^2+x+1)^500 (x^4+...+1)^250: the middle
         coefficient of q^2 is the dot product of q's coefficient list with its
         reverse -- an O(deg) step after one Expand of a degree-3000 power.
   EXACT BASELINE.  Coefficient[Expand[p], x, 3000] (not in the notebook).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Expand", "Coefficient", "CoefficientList", "Series", "SeriesCoefficient", "Residue", "NIntegrate"}];

p = (x + 1)^2000 (x^2 + x + 1)^1000 (x^4 + x^3 + x^2 + x + 1)^500;
q = (x + 1)^1000 (x^2 + x + 1)^500 (x^4 + x^3 + x^2 + x + 1)^250;
exactChk[v_] := If[IntegerQ[v], {sig[v, 15], Mod[v, 10^9]}, "unevaluated"];

run["p4 exact Coefficient[Expand]", Coefficient[Expand[p], x, 3000], exactChk, 3];
run["p4 m4 exact CoefficientList square", (#.Reverse[#] &)[CoefficientList[Expand[q], x]],
    exactChk, 3];
run["p4 m3 exact Series to order 3000",
    SeriesCoefficient[Series[p, {x, 0, 3000}], 3000], exactChk];
run["p4 m2 exact Residue", Residue[p/x^3001, {x, 0}], exactChk];

contour = Evaluate[(p/x^3001 /. x -> E^(I phi)) I E^(I phi)];
run["p4 m1 contour NIntegrate machine",
    Im[NIntegrate[contour, {phi, 0, 2 Pi}]]/(2 Pi), sig[#, 6] &];
run["p4 m1 contour Trapezoidal wp400",
    Im[NIntegrate[contour, {phi, 0, 2 Pi}, Method -> "Trapezoidal",
                  WorkingPrecision -> 400, PrecisionGoal -> 16, MaxRecursion -> 10]]/(2 Pi),
    sig[#, 14] &];
