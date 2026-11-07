(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 182 -- SIAM 100-Digit Challenge, Appendix D Problem 7
   ("Area of a region"; problem 9 on the book's web page).

     f(x,y) = exp(-(y + x^3)^2),  g(x,y) = y^2/32 + exp(sin y).
     What is the area of the region of the x-y plane where f > g?
                                       Answer: 2.1150298809461465858...

   Statement chosen (documented, identical in every system): the 1-D form.
   NIntegrate[Boole[f > g], {x,..}, {y,..}] over a bounding box is the most
   literal WL spelling, but a discontinuous 2-D integrand gives a handful of
   digits at best and no system has a 50-digit route for it.  Since
   f > g  <=>  (y + x^3)^2 < s(y)^2,  s(y) = Sqrt[-Log[g(y)]]  (needs g(y) < 1),
   each horizontal line meets the region in ONE x-interval, of length
       w(y) = Surd[-y + s(y), 3] - Surd[-y - s(y), 3],
   and g(y) < 1 on two "islands" of y: [y1, 0] and [y2, y3].  The case is
     - FindRoot for the island ends y1 ~ -2.85, y2 ~ 3.77, y3 ~ 4.45
       (g(y) = 1) and for the interior kink yc ~ -0.83 where -y - s(y) = 0
       (a cube root of zero: vertical tangent of w), starting values read off
       a plot of g;  y = 0 is an exact island end (g(0) = 1);
     - NIntegrate of w over {y, y1, yc, 0} and {y, y2, y3}.
   w has square-root endpoint singularities at y1, y2, y3 and a |y|^(1/6) one
   at 0.  machine: default goals;  wp50: WorkingPrecision -> 60 in both the
   FindRoots and the NIntegrates, PrecisionGoal -> 50.
   True value and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate", "FindRoot", "Surd"}];

gA[y_] := y^2/32 + Exp[Sin[y]];
wA[y_] := Surd[-y + Sqrt[-Log[gA[y]]], 3] - Surd[-y - Sqrt[-Log[gA[y]]], 3];

areaM[] := Module[{y1, yc, y2, y3},
  y1 = y /. FindRoot[gA[y] == 1, {y, -2.8}];
  yc = y /. FindRoot[-Log[gA[y]] == y^2, {y, -0.8}];
  y2 = y /. FindRoot[gA[y] == 1, {y, 3.8}];
  y3 = y /. FindRoot[gA[y] == 1, {y, 4.45}];
  NIntegrate[wA[y], {y, y1, yc, 0}] + NIntegrate[wA[y], {y, y2, y3}]];

area50[] := Module[{y1, yc, y2, y3},
  y1 = y /. FindRoot[gA[y] == 1, {y, -28/10}, WorkingPrecision -> 60];
  yc = y /. FindRoot[-Log[gA[y]] == y^2, {y, -8/10}, WorkingPrecision -> 60];
  y2 = y /. FindRoot[gA[y] == 1, {y, 38/10}, WorkingPrecision -> 60];
  y3 = y /. FindRoot[gA[y] == 1, {y, 445/100}, WorkingPrecision -> 60];
  NIntegrate[wA[y], {y, y1, yc, 0}, WorkingPrecision -> 60, PrecisionGoal -> 50] +
  NIntegrate[wA[y], {y, y2, y3}, WorkingPrecision -> 60, PrecisionGoal -> 50]];

vrun["d7.machine", areaM[], 10, 3];
vrun["d7.wp50", area50[], 48];
