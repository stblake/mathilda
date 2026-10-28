(* ==========================================================================
   Experiment 121 -- Together, Cancel and Apart at stress sizes
   ==========================================================================

       cd benchmarks/121-rational-functions
       ../../Mathilda -file rational_functions.m
       python3 rational_functions.py
       maxima --very-quiet -b rational_functions.mac
       gp -q rational_functions.gp

   WHAT IT MEASURES.  The rational-function engine (src/poly/ratfun.c) at sizes
   where the algorithm, not the per-call overhead, decides the time: many
   summands over distinct denominators, a nested continued fraction, gcds of
   degree 200 and 500 hidden in degree-400 and degree-1000 quotients, dense and
   sparse multivariate gcds, and partial fractions with many factors, high
   multiplicity, irreducible quadratics and five parameters.  The small cases
   live in 01 (Together of 8, Apart of degree 8), 03 and 10 (two Cancels).

   THE BASELINES.  Maxima (ratsimp, partfrac), sympy (cancel, apart), PARI/GP
   (its native rational functions reduce on construction; it has no general
   partial fractions, so the Apart cases are SKIP there).  sympy's `together`
   does not cancel, so its Together column is `cancel`, which combines and
   reduces in one step -- the same contract.

   THE CHECKS ARE EXACT INTEGERS.  Every value is the result at a rational
   point, so 10^6 v is an exact rational and Round needs no floating point.
   Together and Cancel also pin the denominator's degree in x, which a result
   left unreduced would get wrong:
       check = 10^40 deg_x(Denominator) + Round[10^6 v]      (|10^6 v| < 10^39)
   Apart checks the value only, because systems split the polynomial part and
   the content differently.  The point is x = 11/10 except for the 40-factor
   Apart, whose value there is ~1e-30 (it uses x = 401/10).
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"Together", "Cancel", "Apart", "Denominator", "Exponent"}];

pts = {x -> 11/10, y -> 1/3, z -> 2/7, u -> 3/11, v -> 5/13, w -> 7/17,
       a -> 1/5, b -> 2/9, c -> 3/13, d -> 5/19, f -> 7/23};
chk[t_] := 10^40 Exponent[Denominator[t], x] + Round[10^6 (t /. pts)];
chkA[t_] := Round[10^6 (t /. pts)];
chkA[t_, x0_] := Round[10^6 (t /. x -> x0)];

(* ---- Together ------------------------------------------------------------ *)

tog50 = Sum[(i^2 + 1)/(x - i), {i, 1, 50}];
bench["Together, 50 linear summands", Together[tog50];];
check["Together, 50 linear summands", chk[Together[tog50]]];

tog200 = Sum[(i + 1)/(x - i), {i, 1, 200}];
bench["Together, 200 linear summands", Together[tog200];];
check["Together, 200 linear summands", chk[Together[tog200]]];

togMV = Sum[(x + i y)/((x + i y + z) (y - i z + 1)), {i, 1, 12}];
bench["Together, 12 trivariate summands", Together[togMV];];
check["Together, 12 trivariate summands", chk[Together[togMV]]];

nest32 = Nest[x + 1/(1 + #) &, x, 32];
bench["Together, nested depth 32", Together[nest32];];
check["Together, nested depth 32", chk[Together[nest32]]];

(* ---- Cancel -------------------------------------------------------------- *)

sh200 = factorableUPoly[100, x];
c400n = Expand[sh200 denseUPoly[200, x]];
c400d = Expand[sh200 (denseUPoly[199, x] + 1)];
bench["Cancel deg-400, deg-200 gcd", Cancel[c400n/c400d];];
check["Cancel deg-400, deg-200 gcd", chk[Cancel[c400n/c400d]]];

sh500 = factorableUPoly[250, x];
c1000n = Expand[sh500 denseUPoly[500, x]];
c1000d = Expand[sh500 (denseUPoly[499, x] + 1)];
bench["Cancel deg-1000, deg-500 gcd", Cancel[c1000n/c1000d];];
check["Cancel deg-1000, deg-500 gcd", chk[Cancel[c1000n/c1000d]]];

g3 = (1 + x + y + z)^6;
c3n = Expand[g3 (x + 2 y - z + 3)^4];
c3d = Expand[g3 (x - y + 2 z + 1)^4];
bench["Cancel dense trivariate, deg-6 gcd", Cancel[c3n/c3d];];
check["Cancel dense trivariate, deg-6 gcd", chk[Cancel[c3n/c3d]]];

g6 = 1 + x y^3 + z^2 u v^4 + w^5 + x w;
c6n = Expand[g6^2 (x^2 + z v + 2) (y u + w^2 + 1)];
c6d = Expand[g6^2 (x z + u^3 + 3) (v w + y^2 + 1)];
bench["Cancel sparse 6-variable", Cancel[c6n/c6d];];
check["Cancel sparse 6-variable", chk[Cancel[c6n/c6d]]];

(* ---- Apart --------------------------------------------------------------- *)

ap40 = ratFn[40, x];
bench["Apart, 40 linear factors", Apart[ap40, x];];
check["Apart, 40 linear factors", chkA[Apart[ap40, x], 401/10]];

apMult = (x^3 + 2 x + 5)/((x - 1)^10 (x + 2)^10);
bench["Apart, multiplicity 10 and 10", Apart[apMult, x];];
check["Apart, multiplicity 10 and 10", chkA[Apart[apMult, x]]];

apQuad = denseUPoly[9, x]/Product[x^2 + i, {i, 1, 6}];
bench["Apart, 6 irreducible quadratics", Apart[apQuad, x];];
check["Apart, 6 irreducible quadratics", chkA[Apart[apQuad, x]]];

apPar = 1/((x + a) (x + b) (x + c) (x + d) (x + f));
bench["Apart, 5 parameters", Apart[apPar, x];];
check["Apart, 5 parameters", chkA[Apart[apPar, x]]];
