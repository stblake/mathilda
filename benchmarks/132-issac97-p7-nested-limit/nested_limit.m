(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 132 -- ISSAC '97 Systems Challenge, Problem 7
   ==========================================================================

       cd benchmarks/132-issac97-p7-nested-limit
       ../../Mathilda -file nested_limit.m

   PROBLEM.  Limit[f[g[x]]/g[f[x]], x -> 0] for two signed sums f, g of twelve
   three-fold compositions of Sin/Tan/Sinh/Tanh each.  Massive cancellation:
   f(x) ~ (7769/3274425) x^15 and g(x) ~ (4/63) x^11, so the answer is the
   exact rational ratio of leading terms, 3060352455.6720...
   (= (7769/3274425)(4/63)^15 / ((4/63)(7769/3274425)^11)).

   METHODS (the notebook's two, plus the direct one):
     m1  Series of f and g to order 25, then the ratio of leading terms
         (exact);
     m2  numerically, approaching 0 along 10^-i at 200 i digits until the
         ratio is stationary to 12 digits (cancellation of ~180 digits per
         step -- this is arbitrary-precision evaluation of 24 nested
         elementary functions);
     direct  Limit[f[g[x]]/g[f[x]], x -> 0] (capped).
   No machine variant: at double precision f and g are pure rounding noise
   for every x where the ratio is near its limit.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Series", "Normal", "Limit"}];

f[x_] = Tan[Tanh[Sin[x]]] + Tanh[Sin[Tan[x]]] + Sin[Tan[Tanh[x]]] - Tan[Sin[Tanh[x]]] -
  Sin[Tanh[Tan[x]]] - Tanh[Tan[Sin[x]]] - Tan[Sinh[Tanh[x]]] - Sinh[Tanh[Tan[x]]] -
  Tanh[Tan[Sinh[x]]] + Tan[Tanh[Sinh[x]]] + Tanh[Sinh[Tan[x]]] + Sinh[Tan[Tanh[x]]];
g[x_] = Sinh[Tanh[Sin[x]]] + Tanh[Sin[Sinh[x]]] + Sin[Sinh[Tanh[x]]] - Sinh[Sin[Tanh[x]]] -
  Sin[Tanh[Sinh[x]]] - Tanh[Sinh[Sin[x]]] - Tan[Sinh[Sin[x]]] - Sinh[Sin[Tan[x]]] -
  Sin[Tan[Sinh[x]]] + Tan[Sin[Sinh[x]]] + Sin[Sinh[Tan[x]]] + Sinh[Tan[Sin[x]]];

lead[s_] := First[SortBy[List @@ Expand[s], Exponent[#, x] &]];
serChk[s_] := If[FreeQ[s, Series | SeriesData], {Exponent[lead[s], x], lead[s]/x^Exponent[lead[s], x]},
  "unevaluated"];

(* As in the notebook, the two series are computed once (fSeries, gSeries) and
   the ratio row reuses them -- so it times only the substitution, never a
   system's cache of a Series it has just computed. *)
run["p7 m1 Series f order 25", fSeries = Normal[Series[f[x], {x, 0, 25}]], serChk];
run["p7 m1 Series g order 25", gSeries = Normal[Series[g[x], {x, 0, 25}]], serChk];

m1[] := Module[{fs = lead[fSeries], gs = lead[gSeries]}, (fs /. x -> gs)/(gs /. x -> fs)];
run["p7 m1 Series ratio of leading terms", m1[], sig[#, 12] &, 3];

m2[] := Module[{oldRatio = 10^9, newRatio, i = 2},
  While[newRatio = N[f[g[N[10^(-i), 200 i]]]/g[f[N[10^(-i), 200 i]]], 12];
        Abs[(newRatio - oldRatio)/newRatio] > 10^(-12), i = i + 1; oldRatio = newRatio];
  newRatio];
run["p7 m2 numeric approach wp200i", m2[], sig[#, 11] &, 3];

run["p7 direct Limit", Limit[f[g[x]]/g[f[x]], x -> 0], sig[#, 12] &];
