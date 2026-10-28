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
