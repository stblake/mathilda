(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 122 -- polynomial gcds, lcms and ideals: PolynomialGCD,
   PolynomialLCM, PolynomialReduce, Eliminate, CylindricalDecomposition and
   GroebnerBasis at stress sizes
   ==========================================================================

       cd benchmarks/122-polynomial-ideals
       ../../Mathilda    -file polynomial_ideals.m
       maxima --very-quiet -r 'batchload("polynomial_ideals.mac")$'
       gp -q polynomial_ideals.gp

   WHAT IT MEASURES.  The FLINT-backed polynomial gcd (src/poly/polygcd_flint.c,
   src/poly/poly.c), the Groebner engines (src/poly/groebner.c, gbmod.c and
   gb_flint.c) and their consumers -- multivariate division, elimination -- and
   the cylindrical algebraic decomposition (src/solve/reduce_cad.c).  The small
   cases live in 05, 07 and 10; this one is sized so the algorithm, not the
   per-call overhead, decides the time.  The cases run light to heavy, so a
   timeout loses only the tail.

   THE BASELINES.  Maxima (gcd/lcm; the grobner package: poly_reduced_grobner,
   poly_normal_form, poly_elimination_ideal) and PARI/GP (gcd/lcm only: it has no
   Groebner bases, multivariate division or elimination, so those cases are SKIP
   there).  Neither has a cylindrical decomposition, so those cases are SKIP
   in both.

   THE CHECKS ARE EXACT INTEGERS AND IGNORE NORMALISATION.  Systems disagree on
   the content, sign and monic-ness of a gcd or a basis element, but not on the
   ratio of its values at two rational points, which is invariant under scaling:
       poly:  10^40 deg_x(g) + Round[10^6 g(p1)/g(p2)]
       list:  10^40 Length[G] + Round[10^6 Sum_g g(p1)/g(p2)]
   A reduced Groebner basis for a fixed order is unique, and the remainder of a
   division by one is unique up to scale, so both are pinned.  Under Modulus the
   ratio is taken mod p at two integer points.  A decomposition is checked by
   membership on a fixed grid of rational points (off every boundary), folded
   into a bitmask.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"PolynomialGCD", "PolynomialLCM", "PolynomialReduce", "Eliminate",
         "CylindricalDecomposition", "GroebnerBasis", "MonomialOrder"}];

pt1[vs_] := Thread[vs -> Table[(2 i + 1)/(3 i + 5), {i, Length[vs]}]];
pt2[vs_] := Thread[vs -> Table[(3 i + 2)/(5 i + 7), {i, Length[vs]}]];
q1[vs_] := Thread[vs -> Table[i + 2, {i, Length[vs]}]];
q2[vs_] := Thread[vs -> Table[2 i + 5, {i, Length[vs]}]];
ratio[g_, vs_] := (g /. pt1[vs])/(g /. pt2[vs]);
ratioMod[g_, vs_, p_] := Mod[(g /. q1[vs]) PowerMod[g /. q2[vs], -1, p], p];
polyChk[g_, vs_] := 10^40 Exponent[g, First[vs]] + Round[10^6 ratio[g, vs]];
polyChkMod[g_, vs_, p_] := 10^12 Exponent[g, First[vs]] + ratioMod[g, vs, p];
listChk[G_, vs_] := 10^40 Length[G] + Round[10^6 Total[ratio[#, vs] & /@ G]];
listChkMod[G_, vs_, p_] := 10^12 Length[G] + Mod[Total[ratioMod[#, vs, p] & /@ G], p];
elimPolys[r_] := (Subtract @@ #) & /@ Flatten[{r /. And -> List}];
cadChk[cad_, vs_, grid_] :=
  Total[MapIndexed[Boole[TrueQ[cad /. Thread[vs -> #1]]] 2^(First[#2] - 1) &, grid]];
axis[n_, s_] := Table[s ((2 k - n - 1)/n + 1/97), {k, n}];

katU[vs_, i_] := If[Abs[i] < Length[vs], vs[[Abs[i] + 1]], 0];
katsuraSystem[vs_] := Module[{n = Length[vs] - 1},
  Append[Table[Expand[Sum[katU[vs, l] katU[vs, m - l], {l, -n, n}] - katU[vs, m]],
               {m, 0, n - 1}],
         Sum[katU[vs, l], {l, -n, n}] - 1]];

xs = {x}; v3 = {x, y, z}; v6 = {x, y, z, u, v, w};
p = 32003;

(* ---- PolynomialGCD ------------------------------------------------------- *)

sh100 = factorableUPoly[100, x];
g400a = Expand[sh100 denseUPoly[200, x]];
g400b = Expand[sh100 (denseUPoly[199, x] + 1)];
bench["PolynomialGCD, deg-400, deg-200 gcd", PolynomialGCD[g400a, g400b];];
check["PolynomialGCD, deg-400, deg-200 gcd", polyChk[PolynomialGCD[g400a, g400b], xs]];

sh250 = factorableUPoly[250, x];
g1000a = Expand[sh250 denseUPoly[500, x]];
g1000b = Expand[sh250 (denseUPoly[499, x] + 1)];
bench["PolynomialGCD, deg-1000, deg-500 gcd", PolynomialGCD[g1000a, g1000b];];
check["PolynomialGCD, deg-1000, deg-500 gcd", polyChk[PolynomialGCD[g1000a, g1000b], xs]];

cp1000a = Expand[denseUPoly[1000, x]];
cp1000b = Expand[denseUPoly[999, x] + x^500 + 2];
bench["PolynomialGCD, deg-1000 coprime", PolynomialGCD[cp1000a, cp1000b];];
check["PolynomialGCD, deg-1000 coprime", polyChk[PolynomialGCD[cp1000a, cp1000b], xs]];

shBig = Expand[(10^50 + 3) x^100 + denseUPoly[99, x]];
big200a = Expand[shBig (denseUPoly[100, x] + 7^60)];
big200b = Expand[shBig (denseUPoly[101, x] - 11^50)];
bench["PolynomialGCD, deg-200, 60-digit coefficients", PolynomialGCD[big200a, big200b];];
check["PolynomialGCD, deg-200, 60-digit coefficients",
  polyChk[PolynomialGCD[big200a, big200b], xs]];

g2 = Expand[(x^2 + 3 x y + y^3 + 1)^5];
bi40a = Expand[g2 (x - y^2 + 3)^8];
bi40b = Expand[g2 (x^3 + y + 5)^6];
bench["PolynomialGCD, dense bivariate, deg-15 gcd", PolynomialGCD[bi40a, bi40b];];
check["PolynomialGCD, dense bivariate, deg-15 gcd", polyChk[PolynomialGCD[bi40a, bi40b], v3]];

g3 = (1 + x + y + z)^6;
tri6a = Expand[g3 (x + 2 y - z + 3)^4];
tri6b = Expand[g3 (x - y + 2 z + 1)^4];
bench["PolynomialGCD, dense trivariate, deg-6 gcd", PolynomialGCD[tri6a, tri6b];];
check["PolynomialGCD, dense trivariate, deg-6 gcd", polyChk[PolynomialGCD[tri6a, tri6b], v3]];

g6 = 1 + x y^3 + z^2 u v^4 + w^5 + x w;
sp6a = Expand[g6^2 (x^2 + z v + 2) (y u + w^2 + 1)];
sp6b = Expand[g6^2 (x z + u^3 + 3) (v w + y^2 + 1)];
bench["PolynomialGCD, sparse 6-variable", PolynomialGCD[sp6a, sp6b];];
check["PolynomialGCD, sparse 6-variable", polyChk[PolynomialGCD[sp6a, sp6b], v6]];

ps10 = Table[Expand[(x^2 + y^2 + x y + 1)^2 (x + y + i)^3 (x - i y + 1)^2], {i, 1, 10}];
bench["PolynomialGCD, 10 bivariate arguments", PolynomialGCD @@ ps10;];
check["PolynomialGCD, 10 bivariate arguments", polyChk[PolynomialGCD @@ ps10, v3]];

bench["PolynomialGCD, Modulus 32003, deg-1000",
  PolynomialGCD[g1000a, g1000b, Modulus -> 32003];];
check["PolynomialGCD, Modulus 32003, deg-1000",
  polyChkMod[PolynomialGCD[g1000a, g1000b, Modulus -> 32003], xs, p]];

bench["PolynomialGCD, Modulus 32003, trivariate",
  PolynomialGCD[tri6a, tri6b, Modulus -> 32003];];
check["PolynomialGCD, Modulus 32003, trivariate",
  polyChkMod[PolynomialGCD[tri6a, tri6b, Modulus -> 32003], v3, p]];

(* ---- PolynomialLCM ------------------------------------------------------- *)

bench["PolynomialLCM, deg-400, deg-200 gcd", PolynomialLCM[g400a, g400b];];
check["PolynomialLCM, deg-400, deg-200 gcd", polyChk[PolynomialLCM[g400a, g400b], xs]];

bench["PolynomialLCM, deg-1000, deg-500 gcd", PolynomialLCM[g1000a, g1000b];];
check["PolynomialLCM, deg-1000, deg-500 gcd", polyChk[PolynomialLCM[g1000a, g1000b], xs]];

bench["PolynomialLCM, dense trivariate, deg-6 gcd", PolynomialLCM[tri6a, tri6b];];
check["PolynomialLCM, dense trivariate, deg-6 gcd", polyChk[PolynomialLCM[tri6a, tri6b], v3]];

bench["PolynomialLCM, sparse 6-variable", PolynomialLCM[sp6a, sp6b];];
check["PolynomialLCM, sparse 6-variable", polyChk[PolynomialLCM[sp6a, sp6b], v6]];

bench["PolynomialLCM, 10 bivariate arguments", PolynomialLCM @@ ps10;];
check["PolynomialLCM, 10 bivariate arguments", polyChk[PolynomialLCM @@ ps10, v3]];

bench["PolynomialLCM, Modulus 32003, trivariate",
  PolynomialLCM[tri6a, tri6b, Modulus -> 32003];];
check["PolynomialLCM, Modulus 32003, trivariate",
  polyChkMod[PolynomialLCM[tri6a, tri6b, Modulus -> 32003], v3, p]];

(* ---- PolynomialReduce ---------------------------------------------------- *)

(* The divisors are reduced Groebner bases (computed untimed), so the remainder
   is unique up to scale in every system. *)

c5v = {a1, a2, a3, a4, a5};
cyc5 = cyclicSystem[c5v];
cyc5drl = GroebnerBasis[cyc5, c5v, MonomialOrder -> DegreeReverseLexicographic];
fc5 = Expand[(a1 + 2 a2 + 3 a3 + 4 a4 + 5 a5 + 1)^7];
bench["PolynomialReduce, by the cyclic-5 DRL basis",
  PolynomialReduce[fc5, cyc5drl, c5v, MonomialOrder -> DegreeReverseLexicographic];];
check["PolynomialReduce, by the cyclic-5 DRL basis",
  polyChk[Last[PolynomialReduce[fc5, cyc5drl, c5v,
    MonomialOrder -> DegreeReverseLexicographic]], c5v]];

k4v = {k0, k1, k2, k3, k4};
kat4 = katsuraSystem[k4v];
kat4drl = GroebnerBasis[kat4, k4v, MonomialOrder -> DegreeReverseLexicographic];
fk4 = Expand[(k0 - 2 k1 + 3 k2 - k3 + 2 k4 + 1)^9];
bench["PolynomialReduce, by the katsura-4 DRL basis",
  PolynomialReduce[fk4, kat4drl, k4v, MonomialOrder -> DegreeReverseLexicographic];];
check["PolynomialReduce, by the katsura-4 DRL basis",
  polyChk[Last[PolynomialReduce[fk4, kat4drl, k4v,
    MonomialOrder -> DegreeReverseLexicographic]], k4v]];

sq3 = {x^2 + y + z - 1, x + y^2 + z - 1, x + y + z^2 - 1};
sq3lex = GroebnerBasis[sq3, v3];
f3 = Expand[(x + 2 y + 3 z + 1)^14];
bench["PolynomialReduce, by a 3-variable Lex basis", PolynomialReduce[f3, sq3lex, v3];];
check["PolynomialReduce, by a 3-variable Lex basis",
  polyChk[Last[PolynomialReduce[f3, sq3lex, v3]], v3]];

sq3dl = GroebnerBasis[sq3, v3, MonomialOrder -> DegreeLexicographic];
bench["PolynomialReduce, by a 3-variable DegLex basis",
  PolynomialReduce[f3, sq3dl, v3, MonomialOrder -> DegreeLexicographic];];
check["PolynomialReduce, by a 3-variable DegLex basis",
  polyChk[Last[PolynomialReduce[f3, sq3dl, v3, MonomialOrder -> DegreeLexicographic]], v3]];

(* ---- Eliminate ----------------------------------------------------------- *)

(* Every case has a principal elimination ideal, so the eliminant is one
   polynomial, unique up to scale. *)

bench["Eliminate, 2 variables from 3 quadratics", Eliminate[Thread[sq3 == 0], {y, z}];];
check["Eliminate, 2 variables from 3 quadratics",
  listChk[elimPolys[Eliminate[Thread[sq3 == 0], {y, z}]], xs]];

crv = {x == t^3 - 3 t, y == t^5 + t^2 - 1};
bench["Eliminate, plane curve of degree 5", Eliminate[crv, t];];
check["Eliminate, plane curve of degree 5", listChk[elimPolys[Eliminate[crv, t]], {x, y}]];

srf = {x == s t, y == s^2 - t, z == s + t^3};
bench["Eliminate, parametric surface", Eliminate[srf, {s, t}];];
check["Eliminate, parametric surface", listChk[elimPolys[Eliminate[srf, {s, t}]], v3]];

chn = {x == y^2 + z, y == z^2 + w, z == w^2 + 1, w^3 + x w - 2 == 0};
bench["Eliminate, 3 variables from a 4-variable chain", Eliminate[chn, {y, z, w}];];
check["Eliminate, 3 variables from a 4-variable chain",
  listChk[elimPolys[Eliminate[chn, {y, z, w}]], xs]];

(* ---- CylindricalDecomposition -------------------------------------------- *)

grid2[s_] := Tuples[{axis[7, s], axis[7, s]}];
grid3[s_] := Tuples[{axis[5, s], axis[5, s], axis[5, s]}];

cDisc = x^2 + y^2 < 1;
bench["CylindricalDecomposition, unit disc", CylindricalDecomposition[cDisc, {x, y}];];
check["CylindricalDecomposition, unit disc",
  cadChk[CylindricalDecomposition[cDisc, {x, y}], {x, y}, grid2[6/5]]];

cEll = x^2/4 + y^2 < 1 && x^2 + y^2/4 < 1;
bench["CylindricalDecomposition, two ellipses", CylindricalDecomposition[cEll, {x, y}];];
check["CylindricalDecomposition, two ellipses",
  cadChk[CylindricalDecomposition[cEll, {x, y}], {x, y}, grid2[11/5]]];

cQrt = x^4 + y^4 - 3 x y + x < 1;
bench["CylindricalDecomposition, quartic curve", CylindricalDecomposition[cQrt, {x, y}];];
check["CylindricalDecomposition, quartic curve",
  cadChk[CylindricalDecomposition[cQrt, {x, y}], {x, y}, grid2[8/5]]];

cBall = x^2 + y^2 + z^2 < 1 && x^2 + y^2 < 1/4;
bench["CylindricalDecomposition, ball and cylinder",
  CylindricalDecomposition[cBall, v3];];
check["CylindricalDecomposition, ball and cylinder",
  cadChk[CylindricalDecomposition[cBall, v3], v3, grid3[6/5]]];

cUmb = x^2 < y^2 z && x^2 + y^2 + z^2 < 4;
bench["CylindricalDecomposition, Whitney umbrella in a ball",
  CylindricalDecomposition[cUmb, v3];];
check["CylindricalDecomposition, Whitney umbrella in a ball",
  cadChk[CylindricalDecomposition[cUmb, v3], v3, grid3[11/5]]];

(* ---- GroebnerBasis ------------------------------------------------------- *)

drl = MonomialOrder -> DegreeReverseLexicographic;

bench["GroebnerBasis, 3-variable quadratics, Lex", GroebnerBasis[sq3, v3];];
check["GroebnerBasis, 3-variable quadratics, Lex", listChk[GroebnerBasis[sq3, v3], v3]];

bench["GroebnerBasis, cyclic-5, DRL", GroebnerBasis[cyc5, c5v, drl];];
check["GroebnerBasis, cyclic-5, DRL", listChk[GroebnerBasis[cyc5, c5v, drl], c5v]];

bench["GroebnerBasis, katsura-4, DRL", GroebnerBasis[kat4, k4v, drl];];
check["GroebnerBasis, katsura-4, DRL", listChk[GroebnerBasis[kat4, k4v, drl], k4v]];

bench["GroebnerBasis, elimination form, plane curve",
  GroebnerBasis[{x - t^3 + 3 t, y - t^5 - t^2 + 1}, {x, y}, {t}];];
check["GroebnerBasis, elimination form, plane curve",
  listChk[GroebnerBasis[{x - t^3 + 3 t, y - t^5 - t^2 + 1}, {x, y}, {t}], {x, y}]];

bench["GroebnerBasis, cyclic-5, Modulus 32003",
  GroebnerBasis[cyc5, c5v, Modulus -> 32003];];
check["GroebnerBasis, cyclic-5, Modulus 32003",
  listChkMod[GroebnerBasis[cyc5, c5v, Modulus -> 32003], c5v, p]];

k5v = {k0, k1, k2, k3, k4, k5};
kat5 = katsuraSystem[k5v];
bench["GroebnerBasis, katsura-5, DRL", GroebnerBasis[kat5, k5v, drl];];
check["GroebnerBasis, katsura-5, DRL", listChk[GroebnerBasis[kat5, k5v, drl], k5v]];

bench["GroebnerBasis, cyclic-5, Lex", GroebnerBasis[cyc5, c5v];];
check["GroebnerBasis, cyclic-5, Lex", listChk[GroebnerBasis[cyc5, c5v], c5v]];

c6v = {b1, b2, b3, b4, b5, b6};
cyc6 = cyclicSystem[c6v];
bench["GroebnerBasis, cyclic-6, DRL, Modulus 32003",
  GroebnerBasis[cyc6, c6v, drl, Modulus -> 32003];];
check["GroebnerBasis, cyclic-6, DRL, Modulus 32003",
  listChkMod[GroebnerBasis[cyc6, c6v, drl, Modulus -> 32003], c6v, p]];

bench["GroebnerBasis, cyclic-6, DRL", GroebnerBasis[cyc6, c6v, drl];];
check["GroebnerBasis, cyclic-6, DRL", listChk[GroebnerBasis[cyc6, c6v, drl], c6v]];

bench["GroebnerBasis, katsura-4, Lex", GroebnerBasis[kat4, k4v];];
check["GroebnerBasis, katsura-4, Lex", listChk[GroebnerBasis[kat4, k4v], k4v]];

k6v = {k0, k1, k2, k3, k4, k5, k6};
kat6 = katsuraSystem[k6v];
bench["GroebnerBasis, katsura-6, DRL", GroebnerBasis[kat6, k6v, drl];];
check["GroebnerBasis, katsura-6, DRL", listChk[GroebnerBasis[kat6, k6v, drl], k6v]];

(* ---- GroebnerBasis over Q: the standard families ------------------------- *)

(* The families the Groebner literature measures engines on, in DRL and in
   Lex, zero-dimensional (katsura, cyclic, eco, noon, reimer, the dense random
   systems of the ISSAC '97 holdout, experiment 141) and positive-dimensional
   (the plane-curve elimination above; cyclic-n has positive-dimensional
   components only from n = 8 on).  Light to heavy within each order, the
   heaviest last; a case whose single run takes seconds is timed once. *)

ecoSystem[x_List] := Module[{n = Length[x]},
  Append[Table[Expand[(x[[k]] + Sum[x[[i]] x[[i + k]], {i, 1, n - k - 1}]) x[[n]] - k], {k, 1, n - 1}],
         Sum[x[[i]], {i, 1, n - 1}] + 1]];
noonSystem[x_List] := Table[Expand[10 x[[i]] (Total[x^2] - x[[i]]^2) - 11 x[[i]] + 10], {i, Length[x]}];
reimerSystem[x_List] := Table[Sum[(-1)^(i + 1) 2 x[[i]]^(j + 1), {i, Length[x]}] - 1, {j, 1, Length[x]}];

c7v = {d1, d2, d3, d4, d5, d6, d7}; cyc7 = cyclicSystem[c7v];
k7v = {k0, k1, k2, k3, k4, k5, k6, k7}; kat7 = katsuraSystem[k7v];
k8v = {k0, k1, k2, k3, k4, k5, k6, k7, k8}; kat8 = katsuraSystem[k8v];
e6v = {e1, e2, e3, e4, e5, e6}; eco6 = ecoSystem[e6v];
e7v = {e1, e2, e3, e4, e5, e6, e7}; eco7 = ecoSystem[e7v];
e8v = {e1, e2, e3, e4, e5, e6, e7, e8}; eco8 = ecoSystem[e8v];
n4v = {n1, n2, n3, n4}; noon4 = noonSystem[n4v];
n5v = {n1, n2, n3, n4, n5}; noon5 = noonSystem[n5v];
n6v = {n1, n2, n3, n4, n5, n6}; noon6 = noonSystem[n6v];
r4v = {r1, r2, r3, r4}; reim4 = reimerSystem[r4v];
r5v = {r1, r2, r3, r4, r5}; reim5 = reimerSystem[r5v];
(* experiment 141's h6.1, h6.3, h6.8, h6.10: 4 and 5 dense quadrics with
   1-digit coefficients, 4 with 6-digit ones, 2 quadrics + 2 cubics *)
dv4 = {x, y, z, w}; dv5 = {v, x, y, z, w};
dq4 = {6 + 5*w - w^2 + 2*x - 8*w*x - 3*x^2 - 9*y + 8*w*y + 4*x*y + 9*y^2 - 7*z - 9*w*z + 2*x*z - 9*y*z - 8*z^2, -7 - 5*w - 6*w^2 - 2*x - 4*w*x - 6*x^2 + 9*y + 8*w*y + x*y - y^2 + 8*z - 4*w*z - 2*x*z + 4*y*z - 3*z^2, -5 + 9*w - 7*w^2 - 2*x + 8*w*x - 5*x^2 - 4*y - 7*w*y - 6*x*y + 8*y^2 - 9*z + 5*w*z - 7*x*z - y*z - 8*z^2, 6 + 7*w - 5*w^2 + 6*x - 7*w*x + 7*x^2 - 9*y - 9*w*y - 3*x*y + 9*y^2 + 9*z - 9*w*z + 2*x*z - 2*y*z - 5*z^2};
dq5 = {-3 - 4*v + 7*v^2 - w - 3*v*w + 3*w^2 - 3*x - 9*v*x - 5*w*x + 7*x^2 - 8*y - 5*v*y + 6*w*y + 4*x*y - 3*y^2 - 9*z + 5*v*z - 7*w*z - 9*x*z - 4*y*z - 4*z^2, -4 + 4*v - v^2 - w - v*w - 5*w^2 - 5*x - 7*v*x - w*x - 5*x^2 - 6*y + 9*v*y - 6*w*y - 6*x*y + 9*y^2 - 5*z + 6*v*z + 9*w*z - 8*x*z - 5*y*z - 5*z^2, 8 + 8*v - 6*v^2 + 3*w + 4*v*w + 3*w^2 + 5*x - 4*v*x + 3*w*x - 9*x^2 + 4*y - 2*v*y + 7*w*y + 6*x*y + 8*y^2 + 3*z - 5*v*z + 9*w*z + 3*x*z - 5*y*z + 6*z^2, 9 - 3*v - 2*v^2 + 2*w + 9*v*w + 5*w^2 - 2*x - 2*v*x - 4*w*x + 5*x^2 + 7*y + 5*v*y - 3*w*y + 6*x*y - 3*y^2 - 7*z - 3*v*z - 2*w*z + 6*x*z - 4*y*z + 5*z^2, -9 - v + 3*v^2 - w - 4*v*w + 6*w^2 + 5*x - 5*v*x - 6*w*x + 4*x^2 + 8*y - 6*v*y + w*y - 6*x*y + 7*y^2 - 9*z + 9*v*z + 7*w*z + 8*x*z + 5*y*z + 6*z^2};
dq6d = {854771 - 564557*w - 899147*w^2 - 615329*x - 308540*w*x + 851798*x^2 + 974711*y + 313139*w*y + 222097*x*y - 880766*y^2 - 410696*z - 632850*w*z - 936705*x*z + 914522*y*z + 926361*z^2, -72493 - 960291*w + 809274*w^2 + 656550*x - 65296*w*x - 932720*x^2 - 294045*y - 410727*w*y - 205850*x*y + 888575*y^2 + 802615*z + 708191*w*z + 321448*x*z - 365950*y*z - 614411*z^2, 68559 + 748792*w - 328623*w^2 - 114720*x + 293330*w*x - 223946*x^2 - 443363*y - 191475*w*y - 21700*x*y + 300870*y^2 - 991793*z + 747752*w*z - 948354*x*z - 300782*y*z + 306202*z^2, -469472 - 21909*w + 482639*w^2 + 248427*x + 66724*w*x - 134523*x^2 + 248804*y - 646416*w*y + 612064*x*y + 863308*y^2 - 136206*z - 389466*w*z + 809411*x*z + 574646*y*z + 478494*z^2};
dq23 = {-8 + 3*w + 2*w^2 - 2*x + 5*w*x - 2*x^2 + 9*y - 2*w*y - 6*x*y - y^2 + 7*z + 5*w*z + 4*x*z - 2*y*z + 4*z^2, 2 + 8*w - 5*w^2 + 8*x + 5*w*x + 8*x^2 + y + 9*w*y - 7*x*y - 7*y^2 - 2*z - 8*w*z - 8*x*z + 9*y*z - 5*z^2, 7 - 2*w + 5*w^2 + 5*w^3 - 9*x + 6*w*x - 3*w^2*x + 2*x^2 + 3*w*x^2 - x^3 - 8*y - 6*w*y + 2*w^2*y - 3*x*y - 5*w*x*y + 5*x^2*y + 5*y^2 - 6*w*y^2 - x*y^2 + 6*y^3 + 2*z - 5*w*z - 6*w^2*z - 3*x*z - 8*w*x*z - 8*x^2*z - 6*y*z + 7*w*y*z + x*y*z + 7*y^2*z - z^2 - 6*w*z^2 - 2*x*z^2 - 5*y*z^2 + 3*z^3, 5 - 6*w + 6*w^2 - 2*w^3 - 3*x + 2*w*x - 4*w^2*x - 6*x^2 + 3*w*x^2 - 7*x^3 - 4*y + 4*w*y - 6*w^2*y + 7*x*y - w*x*y + 7*x^2*y + 5*y^2 - 4*w*y^2 - 2*x*y^2 - 3*y^3 + z - 3*w*z - 4*w^2*z + 4*x*z + 4*w*x*z + 4*x^2*z - 6*y*z + 7*w*y*z + 7*x*y*z + 2*y^2*z + 9*z^2 + 3*w*z^2 - 7*x*z^2 - 8*y*z^2 - 7*z^3};

famDRL = {
  {"noon-4", noon4, n4v}, {"eco-6", eco6, e6v}, {"reimer-4", reim4, r4v},
  {"dense 4 quadrics", dq4, dv4}, {"dense 4 quadrics, 6-digit", dq6d, dv4},
  {"dense 5 quadrics", dq5, dv5}, {"dense 2 quadrics + 2 cubics", dq23, dv4},
  {"eco-7", eco7, e7v}, {"noon-5", noon5, n5v}, {"katsura-7", kat7, k7v},
  {"eco-8", eco8, e8v}, {"reimer-5", reim5, r5v}, {"noon-6", noon6, n6v}};
famLex = {
  {"katsura-5", kat5, k5v}, {"reimer-4", reim4, r4v}, {"noon-4", noon4, n4v},
  {"eco-6", eco6, e6v}, {"dense 4 quadrics", dq4, dv4}, {"dense 4 quadrics, 6-digit", dq6d, dv4},
  {"dense 5 quadrics", dq5, dv5}, {"dense 2 quadrics + 2 cubics", dq23, dv4},
  {"eco-7", eco7, e7v}, {"cyclic-6", cyc6, c6v}, {"katsura-6", kat6, k6v},
  {"noon-5", noon5, n5v}, {"reimer-5", reim5, r5v}};

(* the check reads the basis of the last timed run: one computation fewer *)
Do[With[{lab = "GroebnerBasis, " <> f[[1]] <> ", DRL", s = f[[2]], vs = f[[3]]},
   bench[lab, gbRes = GroebnerBasis[s, vs, drl];];
   check[lab, listChk[gbRes, vs]]], {f, famDRL}];
Do[With[{lab = "GroebnerBasis, " <> f[[1]] <> ", Lex", s = f[[2]], vs = f[[3]]},
   bench[lab, gbRes = GroebnerBasis[s, vs];];
   check[lab, listChk[gbRes, vs]]], {f, famLex}];

(* seconds per run: timed once, after the rest *)
benchOnce["GroebnerBasis, katsura-8, DRL", gbRes = GroebnerBasis[kat8, k8v, drl];];
check["GroebnerBasis, katsura-8, DRL", listChk[gbRes, k8v]];
benchOnce["GroebnerBasis, katsura-7, Lex", gbRes = GroebnerBasis[kat7, k7v];];
check["GroebnerBasis, katsura-7, Lex", listChk[gbRes, k7v]];
benchOnce["GroebnerBasis, cyclic-7, DRL", gbRes = GroebnerBasis[cyc7, c7v, drl];];
check["GroebnerBasis, cyclic-7, DRL", listChk[gbRes, c7v]];
