(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 151 -- ISSAC '97 holdout B, family 6: polynomial systems and
   Groebner bases (the counterpart of Problem 6, four dense quadrics)
   ==========================================================================

       cd benchmarks/151-issac97b-polynomial-systems
       ../../Mathilda -file polynomial_systems.m

   The original is one generic dense system in shape position.  These ten
   differ in structure: a structured benchmark system (Katsura-5, 6
   unknowns), three dense cubics, a Lagrange-multiplier system, the
   Arnborg-Lazard system, an inconsistent overdetermined system, an
   elimination (implicitisation) problem, 20-digit coefficients, a system
   whose eliminant has degree 4 but 24 solutions (not in shape position), a
   non-radical system (Cox-Little-O'Shea), and the 6-unknown Trinks system.
   For each: the exact lex GroebnerBasis (checked by its length, the degree of
   its univariate member, and that member's trailing-nonzero/leading
   coefficient ratio, invariants of the reduced basis), and NSolve at machine
   precision and at 40 digits (checked by the number of distinct solutions
   and the largest real part of the last unknown).

   True answers: sympy's Buchberger / FGLM lex bases (an implementation
   independent of both benchmarked systems), the eliminant's roots by mpmath
   polyroots at 60 digits; second methods per problem:
   b6.1  Katsura-5 (u0..u5): 32 solutions (the known 2^n count), eliminant
         degree 32 in u5, max Re u5 = 1/3.
   b6.2  three dense cubics in x, y, z: 27 = 3^3 (Bezout, generic)
         solutions, eliminant degree 27, max Re z = 1.00377595654634977900.
   b6.3  critical points of f = x^3 + 2y^3 + 3z^3 + xyz on
         x^2 + 2y^2 + 3z^2 = 1 (unknowns x, y, z, l): 14 solutions; max Re l
         = 3/2, i.e. max f = 2 l/3 = 1 (Euler's identity), confirmed by direct
         numerical maximisation of f on the ellipsoid (scipy).
   b6.4  Arnborg-Lazard: 20 solutions (literature count), max Re z
         = 9.41270710281753576346.
   b6.5  x^2+y^2+z^2 = 3, xy+yz-2x+1, x^2-yz+3z-2, xz+y^2-x-y-1: inconsistent,
         reduced basis {1} (sympy), NSolve -> no solutions.
   b6.6  implicitise x = t^3 - 3t, y = t^4 - t^2 + t: the eliminant
         -x^4 - 3x^3 - 9x^2 y + x^2 - 12xy + 11x + y^3 - 12y^2 + 33y (9 terms,
         degree 4 in x, 3 in y; equal to the resultant Res_t, sympy); the
         curve meets x^2 + y^2 = 4 in 8 points, max Re y
         = 9.16731117205546875599 (mpmath roots of the degree-8 t-polynomial).
   b6.7  two sextics in x, y with 20-digit coefficients: 36 solutions,
         eliminant degree 36 in y, which divides Res_x (degree 36, squarefree);
         max Re y = 1.33135484579998119546.
   b6.8  power sums p1 = 1, p2 = 3, p3 = 4 and xyzw = -1: the unknowns are the
         roots, in some order, of t^4 - t^3 - t^2 - 1 = (t + 1)(t^3 - 2t^2 + t - 1)
         (Newton's identities, by hand), so 4! = 24 solutions and an
         eliminant of degree 4; max Re w = the real root of t^3 - 2t^2 + t - 1
         = 1.75487766624669276005.
   b6.9  x^2 - 2xz + 5, xy^2 + yz^3, 3y^2 - 8z^3 (CLO, Ideals, Varieties and
         Algorithms, ch. 3): eliminant z^9-degree with a square factor
         (square-free part degree 7), max Re z = 10.4195973828311373737.
   b6.10 Trinks (W, p, z, t, s, b): 10 solutions, eliminant degree 10 in b,
         max Re b = -0.06729047613313375077.

   CHECKS.  Exact: {Length[basis], eliminant degree, sig[ratio, 15]} (b6.5:
   the basis itself; b6.6: {Length, degree in x, degree in y, number of terms,
   leading-coefficient ratio}).  NSolve: {distinct solutions, sig[max Re, d]}
   with d = 8 at machine precision (4 for the 20-digit coefficients of
   b6.7) and 18 at 40 digits.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"GroebnerBasis", "NSolve", "Exponent", "CoefficientList"}];

(* The eliminant: the basis member free of every unknown but the last. *)
uniOf[gb_, rest_] := SelectFirst[gb, FreeQ[#, Alternatives @@ rest] &];
ratio[u_, v_] := Module[{c = DeleteCases[CoefficientList[u, v], 0]}, First[c]/Last[c]];
gbChk[vars_][gb_] := If[ListQ[gb] && FreeQ[gb, GroebnerBasis],
  Module[{u = uniOf[gb, Most[vars]]},
    If[u === Missing["NotFound"] || Head[u] === Missing, "no eliminant",
       {Length[gb], Exponent[u, Last[vars]], sig[ratio[u, Last[vars]], 15]}]], "unevaluated"];
(* distinct solutions (rounded to 1e-6) and the largest real part of the last unknown *)
nsChk[vars_, d_][sols_] := If[ListQ[sols] && FreeQ[sols, NSolve],
  If[sols === {}, {0},
    Module[{pts = vars /. sols},
      {Length[Union[Round[10^6 Join[Re[#], Im[#]]] & /@ pts]], sig[Max[Re[pts[[All, -1]]]], d]}]],
  "unevaluated"];

triple[tag_, eqs_, vars_, dm_] := (
  run[tag <> " GroebnerBasis exact", GroebnerBasis[eqs, vars], gbChk[vars], 3];
  run[tag <> " NSolve wp40", NSolve[eqs, vars, WorkingPrecision -> 40], nsChk[vars, 18]];
  run[tag <> " NSolve machine", NSolve[eqs, vars], nsChk[vars, dm], 3]);
SetAttributes[triple, HoldAll];

(* ---- b6.1 Katsura-5 --------------------------------------------------- *)
katsura = {u0^2 - u0 + 2 u1^2 + 2 u2^2 + 2 u3^2 + 2 u4^2 + 2 u5^2,
  2 u0 u1 + 2 u1 u2 - u1 + 2 u2 u3 + 2 u3 u4 + 2 u4 u5,
  2 u0 u2 + u1^2 + 2 u1 u3 + 2 u2 u4 - u2 + 2 u3 u5,
  2 u0 u3 + 2 u1 u2 + 2 u1 u4 + 2 u2 u5 - u3,
  2 u0 u4 + 2 u1 u3 + 2 u1 u5 + u2^2 - u4,
  u0 + 2 u1 + 2 u2 + 2 u3 + 2 u4 + 2 u5 - 1};
triple["b6.1", katsura, {u0, u1, u2, u3, u4, u5}, 8];

(* ---- b6.2 three dense cubics ------------------------------------------ *)
cubics = {x^3 + 2 x^2 y - 3 x y z + y^3 - 2 z^3 + 4 x z - y^2 + 5 z - 3,
  2 x^3 - y^3 + x y^2 + 3 y z^2 - z^3 + x^2 - 2 y z + x - 1,
  -x^3 + 3 x z^2 + y^2 z + 2 z^3 - 4 x y + y^2 - 2 z + 2};
triple["b6.2", cubics, {x, y, z}, 8];

(* ---- b6.3 Lagrange critical points ------------------------------------ *)
fL = x^3 + 2 y^3 + 3 z^3 + x y z; gL = x^2 + 2 y^2 + 3 z^2 - 1;
lagr = Append[Table[D[fL, v] - l D[gL, v], {v, {x, y, z}}], gL];
triple["b6.3", lagr, {x, y, z, l}, 8];

(* ---- b6.4 Arnborg-Lazard ---------------------------------------------- *)
arnborg = {x^2 y z + x y^2 z + x y z^2 + x y z + x y + x z + y z,
  x^2 y^2 z + x y^2 z^2 + x^2 y z + x y z + y z + x + z,
  x^2 y^2 z^2 + x^2 y^2 z + x y^2 z + x y z + x z + z + 1};
triple["b6.4", arnborg, {x, y, z}, 8];

(* ---- b6.5 inconsistent ------------------------------------------------ *)
incons = {x^2 + y^2 + z^2 - 3, x y + y z - 2 x + 1, x^2 - y z + 3 z - 2, x z + y^2 - x - y - 1};
run["b6.5 GroebnerBasis exact", GroebnerBasis[incons, {x, y, z}], If[ListQ[#], #, "unevaluated"] &, 3];
run["b6.5 NSolve wp40", NSolve[incons, {x, y, z}, WorkingPrecision -> 40], nsChk[{x, y, z}, 18]];
run["b6.5 NSolve machine", NSolve[incons, {x, y, z}], nsChk[{x, y, z}, 8], 3];

(* ---- b6.6 implicitisation + intersection ------------------------------ *)
implChk[gb_] := If[ListQ[gb] && Length[gb] >= 1 && FreeQ[gb, GroebnerBasis | t],
  Module[{g = Expand[First[gb]]},
    {Length[gb], Exponent[g, x], Exponent[g, y], Length[g], Coefficient[g, y, 3]/Coefficient[g, x, 4]}], "unevaluated"];
run["b6.6 GroebnerBasis eliminate exact", GroebnerBasis[{x - (t^3 - 3 t), y - (t^4 - t^2 + t)}, {x, y}, {t}], implChk, 3];
curve = {x - (t^3 - 3 t), y - (t^4 - t^2 + t), x^2 + y^2 - 4};
run["b6.6 NSolve wp40", NSolve[curve, {t, x, y}, WorkingPrecision -> 40], nsChk[{t, x, y}, 18]];
run["b6.6 NSolve machine", NSolve[curve, {t, x, y}], nsChk[{t, x, y}, 8], 3];

(* ---- b6.7 20-digit coefficients --------------------------------------- *)
bigc = {123456789012345678901 x^6 - 98765432109876543210 x^3 y^2 + 31415926535897932384 y^6 -
    27182818284590452353 x y + 11111111111111111111,
  x^5 y - 16180339887498948482 x^2 y^4 + 14142135623730950488 y^3 - 17320508075688772935 x^4 +
    22360679774997896964 y - 1};
triple["b6.7", bigc, {x, y}, 4];

(* ---- b6.8 power sums (24 solutions, eliminant degree 4) --------------- *)
psums = {x + y + z + w - 1, x^2 + y^2 + z^2 + w^2 - 3, x^3 + y^3 + z^3 + w^3 - 4, x y z w + 1};
triple["b6.8", psums, {x, y, z, w}, 8];

(* ---- b6.9 non-radical (CLO) ------------------------------------------- *)
clo = {x^2 - 2 x z + 5, x y^2 + y z^3, 3 y^2 - 8 z^3};
triple["b6.9", clo, {x, y, z}, 8];

(* ---- b6.10 Trinks ----------------------------------------------------- *)
trinks = {45 p + 35 s - 165 b - 36, 35 p + 40 zz + 25 tt - 27 s, 15 ww + 25 p s + 30 zz - 18 tt - 165 b^2,
  -9 ww + 15 p tt + 20 zz s, ww p + 2 zz tt - 11 b^3, 99 ww - 11 b s + 3 b^2};
triple["b6.10", trinks, {ww, p, zz, tt, s, b}, 8];
