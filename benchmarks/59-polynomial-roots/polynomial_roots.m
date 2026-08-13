(* ==========================================================================
   Experiment 59 -- High-degree polynomial roots (group E)
   ==========================================================================
   WHAT IT MEASURES.  NSolve / NRoots on dense univariate polynomials --
   src/numerical_roots/nroots.c, nroots_aberth.c (Aberth-Ehrlich), nroots_jt.c
   (Jenkins-Traub).  Experiment 14 finds a single scalar root; this finds ALL
   roots of degree 20 -> 100, the many-root regime where the companion-matrix /
   simultaneous-iteration methods live.

   Baseline is numpy.roots (companion-matrix eigenvalues via LAPACK), so both
   sides compute the same roots and the order-invariant check Total[|root|]
   agrees to 4 places even at degree 100.  The dense coefficients come from
   denseUPoly (data.m), shared with the Python side; roots-of-unity and the
   Wilkinson polynomial add a unit-circle case and an ill-conditioned case.

   Checks: sum of |root| (order-invariant, exercises the root-finder, unlike the
   Vieta-trivial sum of roots), or the root count.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"NSolve", "NRoots"}];

p20 = denseUPoly[20, x];
p50 = denseUPoly[50, x];
p100 = denseUPoly[100, x];
p40 = denseUPoly[40, x];
pRootUnity = x^64 - 1;
pWilk = Product[(x - i), {i, 1, 15}];

bench["NSolve dense degree 20", NSolve[p20 == 0, x];];
check["NSolve dense degree 20", Round[10^4 Total[Abs[x /. NSolve[p20 == 0, x]]]]];

bench["NSolve dense degree 50", NSolve[p50 == 0, x];];
check["NSolve dense degree 50", Round[10^4 Total[Abs[x /. NSolve[p50 == 0, x]]]]];

bench["NSolve dense degree 100", NSolve[p100 == 0, x];];
check["NSolve dense degree 100", Round[10^4 Total[Abs[x /. NSolve[p100 == 0, x]]]]];

bench["NSolve roots of unity degree 64", NSolve[pRootUnity == 0, x];];
check["NSolve roots of unity degree 64",
      Round[10^4 Total[Abs[x /. NSolve[pRootUnity == 0, x]]]]];

bench["NSolve Wilkinson degree 15", NSolve[pWilk == 0, x];];
check["NSolve Wilkinson degree 15", Round[10^3 Total[Abs[x /. NSolve[pWilk == 0, x]]]]];

bench["NRoots dense degree 40", NRoots[p40 == 0, x];];
check["NRoots dense degree 40", Length[NRoots[p40 == 0, x]]];
