(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 161 -- ISSAC '97 holdout C (SEALED), family c6:
   polynomial systems (counterpart of Problem 6)
   ==========================================================================

       cd benchmarks/161-issac97c-polynomial-systems
       ../../Mathilda -file polynomial_systems.m

   Ten zero-dimensional systems from applications: totally mixed Nash
   equilibria of a 3-player and of a 4-player 2-strategy game (multilinear
   indifference conditions), a mass-action chemical equilibrium with
   conservation laws, the forward kinematics of a planar 3-RPR parallel robot
   (assembly modes; rotation as (c, s) with c^2 + s^2 = 1), the steady states
   of a three-gene repressilator, two dense quartic curves (Bezout 16), the
   Lagrange system for max x + 2y + 3z on x^4 + y^4 + z^4 = 1, the Noonburg
   neural-network system (NOT in shape position: 21 solutions over 14 values of
   the last variable), the feet of the normals from (3, 3) to a cubic curve,
   and the steady states of two diffusively coupled Brusselator cells.

   THREE ROWS PER SYSTEM.
     exact    GroebnerBasis[eqs, vars] (lex, vars as listed): the number of
              basis elements, the degree of the univariate member in the last
              variable, and 15 digits of its constant/leading coefficient
              ratio (normalisation independent).
     machine  NSolve[eqs, vars]: the number of solutions and 8 digits of the
              largest real part of the last variable over all solutions.
     wp40     the same with WorkingPrecision -> 40, 20 digits.

   TRUE ANSWERS (A = sympy's lex Groebner basis over Q, the eliminant's roots
   at 80 digits by mpmath polyroots, every solution back-substituted into the
   original equations: residual <= 2e-74; B = Maxima's poly_reduced_grobner,
   an independent Buchberger implementation, whose basis has the same length
   and whose univariate member equals A's up to a scalar).  The solution count
   is the eliminant degree where the basis is in shape position (all but
   c6.7; every eliminant is squarefree); for c6.7 it is the dimension of the
   quotient ring (21 standard monomials), radical because the eliminant of
   the generic form x + 3y - 2z has degree 21 and is squarefree.

        #   vars          len deg  c0/cN                          #sol  max Re(last var)
   c6.1  {p1,p2,p3}       3   2   1/3                              2   2/3
   c6.2  {a,b,c,d}        4   4   9/4                              4   1.65840319545635272232263624243
   c6.3  {x,y,c,s}        4   4   -0.531615234375                  4   0.984693406500657076632830920783
   c6.4  {u,v,w}          3   9   -130                             9   1.51598022769282058968467151318
   c6.5  {x,y}            2  16   3103                            16   1.79759338472934391991989819789
   c6.6  {x,y,z,lam}      4  36   -1.575841326647290296589e-4     36   1.17206247758678450467230451977
   c6.7  {x,y,z}          4  14   0.1375                          21   1.68372096585234179107896172109
   c6.8  {x,y}            2   7   5                                7   3.23552961063904388209241601064
   c6.9  {x1,y1,x2,y2}    4   5   -357.65625                       5   6.27469440017220052664355112728
   c6.10 {p1,p2,p3,p4}    4   9   -0.07032260594940866061513273    9   3.93580934808760382950002583463
   c6.6 also has a closed form: the maximum of x + 2y + 3z is 4 lam, and by
   Hoelder it is (1 + 2^(4/3) + 3^(4/3))^(3/4), so max lam = 1.172062477586...
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"GroebnerBasis", "NSolve", "Exponent", "Coefficient"}];

sys[1] = {{1 - 3 p2 - 2 p3 + 5 p2 p3, -2 + 3 p1 + 4 p3 - 6 p1 p3, 1 - 2 p1 - 3 p2 + 5 p1 p2}, {p1, p2, p3}};
sys[2] = {{c - 3 a b, d - 5 c^2, a + c + 2 d - 2, b + c + 2 d - 3}, {a, b, c, d}};
sys[3] = {{x^2 + y^2 - 50, (x + 5 c - 10)^2 + (y + 5 s)^2 - 40,
           (x + 2 c - 4 s - 4)^2 + (y + 2 s + 4 c - 8)^2 - 45, c^2 + s^2 - 1}, {x, y, c, s}};
sys[4] = {{u (1 + w^2) - 5, v (1 + u^2) - 5, w (1 + v^2) - 5}, {u, v, w}};
sys[5] = {{x^4 - 3 x^2 y + y^3 - 2 x y + x - 5, y^4 + 2 x^3 - x y^2 + 3 y - 7}, {x, y}};
sys[6] = {{4 lam x^3 - 1, 4 lam y^3 - 2, 4 lam z^3 - 3, x^4 + y^4 + z^4 - 1}, {x, y, z, lam}};
sys[7] = {{x (y^2 + z^2) - 11/10 x + 1, y (x^2 + z^2) - 11/10 y + 1, z (x^2 + y^2) - 11/10 z + 1}, {x, y, z}};
sys[8] = {{y^2 - x^3 + 2 x - 1, 2 y (x - 3) - (y - 3) (2 - 3 x^2)}, {x, y}};
sys[9] = {{2 - 6 x1 + x1^2 y1 + (x2 - x1), 5 x1 - x1^2 y1 + (y2 - y1),
           2 - 6 x2 + x2^2 y2 + (x1 - x2), 5 x2 - x2^2 y2 + (y1 - y2)}, {x1, y1, x2, y2}};
sys[10] = {{2 p2 p3 p4 - 7 p2 p3 + 8 p2 p4 - 5 p2 - 6 p3 p4 + 3 p3 - 8 p4 + 1,
            4 p1 p3 p4 - 8 p1 p3 - 7 p1 p4 - 8 p1 + 4 p3 p4 + 7 p3 - 3 p4 + 9,
            -6 p1 p2 p4 + 4 p1 p2 - 8 p1 p4 - 2 p1 + 9 p2 p4 - 7 p2 + 8 p4 - 7,
            -2 p1 p2 p3 + 9 p1 p2 + 3 p1 p3 + 9 p1 - 8 p2 p3 - 8 p2 + 9 p3 - 2}, {p1, p2, p3, p4}};

uniOf[gb_, vs_] := Select[gb, FreeQ[#, Alternatives @@ Most[vs]] &];
gbRow[gb_, l_, u_] := Module[{n = Exponent[u, l]},
  {Length[gb], n, sig[Coefficient[u, l, 0]/Coefficient[u, l, n], 15]}];
gbChk[vs_][gb_] := Which[
  !ListQ[gb], "unevaluated",
  Length[uniOf[gb, vs]] != 1, "unevaluated",
  True, gbRow[gb, Last[vs], Expand[First[uniOf[gb, vs]]]]];
nsChk[vs_, d_][s_] := If[ListQ[s] && Length[s] > 0,
  {Length[s], sig[Max[Re[Last[vs] /. s]], d]}, "unevaluated"];

Do[
  With[{eqs = sys[k][[1]], vs = sys[k][[2]], lab = "c6." <> ToString[k]},
    run[lab <> " GroebnerBasis exact", GroebnerBasis[eqs, vs], gbChk[vs]];
    run[lab <> " NSolve machine", NSolve[eqs, vs], nsChk[vs, 8]];
    run[lab <> " NSolve wp40", NSolve[eqs, vs, WorkingPrecision -> 40], nsChk[vs, 20]]],
  {k, 10}];
