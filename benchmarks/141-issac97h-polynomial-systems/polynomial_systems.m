(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 141 -- ISSAC '97 holdout, family 6: polynomial systems
   ==========================================================================

       cd benchmarks/141-issac97h-polynomial-systems
       ../../Mathilda -file polynomial_systems.m

   HOLDOUT.  Ten problems probing what experiment 131 (the lex Groebner basis
   of four dense quadrics) probes -- exact Groebner bases of zero-dimensional
   systems and their numerical solutions -- on different systems: new dense
   random-coefficient systems of the same and larger size, mixed degrees,
   six-digit coefficients, a grevlex basis, and two classical benchmarks
   (Katsura-5, cyclic-5).  See benchmarks/ISSAC97_HOLDOUT.md.

   Coefficients of s1-s4 and s7-s10 are uniform random nonzero integers
   (Python random.seed(97); [-9, 9], s8 [-999999, 999999]), fixed here.

   Each problem has three cases:
     exact    GroebnerBasis[sys, vars] (lex; s4 grevlex);
     machine  NSolve[sys, vars];
     wp40     NSolve[sys, vars, WorkingPrecision -> 40].
   CHECKS.  exact: the number of basis elements, the sorted exponent vectors
   (Exponent[g, vars]) and, where the basis has a univariate member in the
   last variable, its constant/leading coefficient ratio to 16 digits.
   numeric: the number of solutions and Max[Re[v]] over them for the check
   variable v (machine 8 digits, wp40 18 digits).

   PROBLEMS AND TRUE VALUES.  The ratio was obtained from the exact basis and
   independently as (-1)^d times the product of the last coordinates of the
   d solutions from a 60-digit NSolve -- the univariate member of a shape
   basis is c Product[(v - v_i)] -- and the two agree to >= 19 digits; the
   solution count is the Bezout number for the dense generic systems and the
   known count (32, 70) for Katsura-5 and cyclic-5; Max Re agrees between a
   60-digit and a machine NSolve to >= 10 digits.
     sys  shape                          solutions  ratio (const/lead)       Max Re[v]
   h6.1  4 quadrics in x,y,z,w            16  deg 16  0.8696578864525617678  w 0.4997973735980418593
   h6.2  3 cubics in x,y,z                27  deg 27  0.01501362839974060569 z 3.134381580032942796
   h6.3  5 quadrics in v,x,y,z,w          32  deg 32  192.0781013427937613   w 3.922013766919696123
   h6.4  4 quadrics, grevlex: 12 basis elements      (no univariate member) 16 w 0.9016989787499557851
   h6.5  Katsura-5 (u1..u5, u0)           32  deg 32  9.976845909679557194e-11 u1 (see table)
   h6.6  cyclic-5 (a..e)                  70  deg 15  -1                     e 2.118033988749894848
   h6.7  degree 7 and 8 in x,y            56  deg 56  1.393462280185198192   y 1.599325999854060862
   h6.8  4 quadrics, 6-digit coefficients 16  deg 16  36.78270028974930448   w 3.359152893765975268
   h6.9  3 quadrics + 1 cubic in x,y,z,w  24  deg 24  7.172513660038497679   w 0.9978641849262317102
   h6.10 2 quadrics + 2 cubics            36  deg 36  -32.06374098385247742  w 3.388175712657317654
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"GroebnerBasis", "NSolve", "Exponent"}];

s1 = {6 + 5*w - w^2 + 2*x - 8*w*x - 3*x^2 - 9*y + 8*w*y + 4*x*y + 9*y^2 - 7*z - 9*w*z + 2*x*z - 9*y*z - 8*z^2, -7 - 5*w - 6*w^2 - 2*x - 4*w*x - 6*x^2 + 9*y + 8*w*y + x*y - y^2 + 8*z - 4*w*z - 2*x*z + 4*y*z - 3*z^2, -5 + 9*w - 7*w^2 - 2*x + 8*w*x - 5*x^2 - 4*y - 7*w*y - 6*x*y + 8*y^2 - 9*z + 5*w*z - 7*x*z - y*z - 8*z^2, 6 + 7*w - 5*w^2 + 6*x - 7*w*x + 7*x^2 - 9*y - 9*w*y - 3*x*y + 9*y^2 + 9*z - 9*w*z + 2*x*z - 2*y*z - 5*z^2};
s2 = {-4 + 4*x + 8*x^2 - 9*x^3 + 4*y + 8*x*y - 8*x^2*y + 9*y^2 + 5*x*y^2 - 8*y^3 + 2*z + 8*x*z + 8*x^2*z - 2*y*z + 4*x*y*z + 7*y^2*z + 8*z^2 - 6*x*z^2 - 9*y*z^2 - 8*z^3, -9 - 3*x + 7*x^2 + 5*x^3 + 6*y + 8*x*y - 6*x^2*y + 9*y^2 + x*y^2 - 2*y^3 - 7*z + 7*x*z + 2*x^2*z - 6*y*z - 8*x*y*z - 4*y^2*z + 9*z^2 - 3*x*z^2 + y*z^2 - 7*z^3, 3 - 7*x - 6*x^2 + 8*x^3 + 7*y + 8*x*y + 5*x^2*y + 3*y^2 + 5*x*y^2 - 4*y^3 - 4*z - x*z - 9*x^2*z - 3*y*z + 7*x*y*z + 8*y^2*z + 2*z^2 - 8*x*z^2 + 7*y*z^2 - 8*z^3};
s3 = {-3 - 4*v + 7*v^2 - w - 3*v*w + 3*w^2 - 3*x - 9*v*x - 5*w*x + 7*x^2 - 8*y - 5*v*y + 6*w*y + 4*x*y - 3*y^2 - 9*z + 5*v*z - 7*w*z - 9*x*z - 4*y*z - 4*z^2, -4 + 4*v - v^2 - w - v*w - 5*w^2 - 5*x - 7*v*x - w*x - 5*x^2 - 6*y + 9*v*y - 6*w*y - 6*x*y + 9*y^2 - 5*z + 6*v*z + 9*w*z - 8*x*z - 5*y*z - 5*z^2, 8 + 8*v - 6*v^2 + 3*w + 4*v*w + 3*w^2 + 5*x - 4*v*x + 3*w*x - 9*x^2 + 4*y - 2*v*y + 7*w*y + 6*x*y + 8*y^2 + 3*z - 5*v*z + 9*w*z + 3*x*z - 5*y*z + 6*z^2, 9 - 3*v - 2*v^2 + 2*w + 9*v*w + 5*w^2 - 2*x - 2*v*x - 4*w*x + 5*x^2 + 7*y + 5*v*y - 3*w*y + 6*x*y - 3*y^2 - 7*z - 3*v*z - 2*w*z + 6*x*z - 4*y*z + 5*z^2, -9 - v + 3*v^2 - w - 4*v*w + 6*w^2 + 5*x - 5*v*x - 6*w*x + 4*x^2 + 8*y - 6*v*y + w*y - 6*x*y + 7*y^2 - 9*z + 9*v*z + 7*w*z + 8*x*z + 5*y*z + 6*z^2};
s4 = {-9 - 3*w + 6*w^2 + 8*x + 2*w*x + 9*x^2 + 4*y + w*y - x*y - 3*y^2 + z + 3*w*z + 7*x*z - 2*y*z - 7*z^2, 7 - 3*w + 7*w^2 - x + 7*w*x + x^2 - 7*y - 8*w*y - x*y + y^2 + 9*z - 3*w*z - 5*x*z - 3*y*z - 6*z^2, -2 + 9*w + 4*w^2 - 2*x - 2*w*x - 2*x^2 + 8*y - 4*w*y - x*y + 7*y^2 + 8*z - 9*w*z - 6*x*z - 4*y*z + 7*z^2, 9 - 8*w + w^2 - 6*x - 6*w*x + 6*x^2 + 2*y - 2*w*y + 8*x*y + 8*y^2 + 7*z + 4*w*z + 4*x*z + 4*y*z - 2*z^2};
s5 = {-1 + u0 + 2*u1 + 2*u2 + 2*u3 + 2*u4 + 2*u5, -u0 + u0^2 + 2*u1^2 + 2*u2^2 + 2*u3^2 + 2*u4^2 + 2*u5^2, -u1 + 2*u0*u1 + 2*u1*u2 + 2*u2*u3 + 2*u3*u4 + 2*u4*u5, u1^2 - u2 + 2*u0*u2 + 2*u1*u3 + 2*u2*u4 + 2*u3*u5, 2*u1*u2 - u3 + 2*u0*u3 + 2*u1*u4 + 2*u2*u5, u2^2 + 2*u1*u3 - u4 + 2*u0*u4 + 2*u1*u5};
s6 = {a + b + c + d + e, a*b + b*c + c*d + a*e + d*e, a*b*c + b*c*d + a*b*e + a*d*e + c*d*e, a*b*c*d + a*b*c*e + a*b*d*e + a*c*d*e + b*c*d*e, -1 + a*b*c*d*e};
s7 = {6 + 8*x - x^2 + 3*x^3 - 8*x^4 + x^5 + 2*x^6 - 3*x^7 - 2*y - 7*x*y - 5*x^2*y - 5*x^3*y - 3*x^4*y + 8*x^5*y - 7*x^6*y - 9*y^2 - 5*x*y^2 - 6*x^2*y^2 + 9*x^3*y^2 + 3*x^4*y^2 - 9*x^5*y^2 + 8*y^3 - 5*x*y^3 - 5*x^2*y^3 + 8*x^3*y^3 - 3*x^4*y^3 - 9*y^4 - 5*x*y^4 - 4*x^2*y^4 - 5*x^3*y^4 - 4*y^5 - 4*x*y^5 - x^2*y^5 + 2*y^6 + 5*x*y^6 - 8*y^7, -7 - 2*x + 7*x^2 + 8*x^3 - 8*x^4 - 4*x^5 - 7*x^6 - 7*x^7 + 2*x^8 + 5*y - 7*x*y + 3*x^2*y - 3*x^3*y - 8*x^4*y + 8*x^5*y - 4*x^6*y - 6*x^7*y + 6*y^2 - 3*x*y^2 + 8*x^2*y^2 - 5*x^3*y^2 - 8*x^4*y^2 + 2*x^5*y^2 - 9*x^6*y^2 - y^3 + 4*x*y^3 - 4*x^2*y^3 + 2*x^3*y^3 + 9*x^4*y^3 - 8*x^5*y^3 - 3*y^4 - 5*x*y^4 + 2*x^2*y^4 + x^3*y^4 - 2*x^4*y^4 - 9*y^5 + 3*x*y^5 - 4*x^2*y^5 - 9*x^3*y^5 + 4*y^6 - 3*x*y^6 - 5*x^2*y^6 + 7*y^7 + 5*x*y^7 + y^8};
s8 = {854771 - 564557*w - 899147*w^2 - 615329*x - 308540*w*x + 851798*x^2 + 974711*y + 313139*w*y + 222097*x*y - 880766*y^2 - 410696*z - 632850*w*z - 936705*x*z + 914522*y*z + 926361*z^2, -72493 - 960291*w + 809274*w^2 + 656550*x - 65296*w*x - 932720*x^2 - 294045*y - 410727*w*y - 205850*x*y + 888575*y^2 + 802615*z + 708191*w*z + 321448*x*z - 365950*y*z - 614411*z^2, 68559 + 748792*w - 328623*w^2 - 114720*x + 293330*w*x - 223946*x^2 - 443363*y - 191475*w*y - 21700*x*y + 300870*y^2 - 991793*z + 747752*w*z - 948354*x*z - 300782*y*z + 306202*z^2, -469472 - 21909*w + 482639*w^2 + 248427*x + 66724*w*x - 134523*x^2 + 248804*y - 646416*w*y + 612064*x*y + 863308*y^2 - 136206*z - 389466*w*z + 809411*x*z + 574646*y*z + 478494*z^2};
s9 = {-7 - 2*w - w^2 + x + 4*w*x + 8*x^2 - 8*y + 2*w*y + 7*x*y + 4*y^2 - 5*z + 5*w*z + 2*x*z + 8*y*z + 8*z^2, 2 + 4*w - 7*w^2 - 7*x + 6*w*x - 8*x^2 + 9*y + 7*w*y + 9*x*y - 2*y^2 - 3*z + 6*w*z - 7*x*z + 5*y*z - z^2, -4 + 7*w - 3*w^2 - 2*x + 9*w*x - 4*x^2 + y - 6*w*y + 9*x*y + 7*y^2 - 6*z + 8*w*z + 2*x*z + 4*y*z - 6*z^2, 5 + 9*w + 2*w^2 + 9*w^3 - 2*x - 8*w*x - 5*w^2*x + 9*x^2 - 2*w*x^2 + 2*x^3 - 3*y - 8*w*y - 6*w^2*y - 2*x*y - 4*w*x*y - x^2*y + 6*y^2 + 5*w*y^2 - 7*x*y^2 + 9*y^3 - 2*z + w*z - 4*w^2*z + 9*x*z + 7*w*x*z + 4*x^2*z + 8*y*z - 3*w*y*z + 3*x*y*z + 4*y^2*z - 4*z^2 + 9*w*z^2 + 3*x*z^2 - 2*y*z^2 - 5*z^3};
s10 = {-8 + 3*w + 2*w^2 - 2*x + 5*w*x - 2*x^2 + 9*y - 2*w*y - 6*x*y - y^2 + 7*z + 5*w*z + 4*x*z - 2*y*z + 4*z^2, 2 + 8*w - 5*w^2 + 8*x + 5*w*x + 8*x^2 + y + 9*w*y - 7*x*y - 7*y^2 - 2*z - 8*w*z - 8*x*z + 9*y*z - 5*z^2, 7 - 2*w + 5*w^2 + 5*w^3 - 9*x + 6*w*x - 3*w^2*x + 2*x^2 + 3*w*x^2 - x^3 - 8*y - 6*w*y + 2*w^2*y - 3*x*y - 5*w*x*y + 5*x^2*y + 5*y^2 - 6*w*y^2 - x*y^2 + 6*y^3 + 2*z - 5*w*z - 6*w^2*z - 3*x*z - 8*w*x*z - 8*x^2*z - 6*y*z + 7*w*y*z + x*y*z + 7*y^2*z - z^2 - 6*w*z^2 - 2*x*z^2 - 5*y*z^2 + 3*z^3, 5 - 6*w + 6*w^2 - 2*w^3 - 3*x + 2*w*x - 4*w^2*x - 6*x^2 + 3*w*x^2 - 7*x^3 - 4*y + 4*w*y - 6*w^2*y + 7*x*y - w*x*y + 7*x^2*y + 5*y^2 - 4*w*y^2 - 2*x*y^2 - 3*y^3 + z - 3*w*z - 4*w^2*z + 4*x*z + 4*w*x*z + 4*x^2*z - 6*y*z + 7*w*y*z + 7*x*y*z + 2*y^2*z + 9*z^2 + 3*w*z^2 - 7*x*z^2 - 8*y*z^2 - 7*z^3};

(* {label, system, variables (lex order, last = eliminated last), monomial order, check variable} *)
problems = {
  {"h6.1", s1, {x, y, z, w}, Lexicographic, w},
  {"h6.2", s2, {x, y, z}, Lexicographic, z},
  {"h6.3", s3, {v, x, y, z, w}, Lexicographic, w},
  {"h6.4", s4, {x, y, z, w}, DegreeReverseLexicographic, w},
  {"h6.5", s5, {u1, u2, u3, u4, u5, u0}, Lexicographic, u1},
  {"h6.6", s6, {a, b, c, d, e}, Lexicographic, e},
  {"h6.7", s7, {x, y}, Lexicographic, y},
  {"h6.8", s8, {x, y, z, w}, Lexicographic, w},
  {"h6.9", s9, {x, y, z, w}, Lexicographic, w},
  {"h6.10", s10, {x, y, z, w}, Lexicographic, w}};

gbChk[gb_, vs_] := If[!ListQ[gb] || !FreeQ[gb, GroebnerBasis], "unevaluated",
  Module[{lv = Last[vs], uni = Select[gb, FreeQ[#, Alternatives @@ Most[vs]] &]},
    {Length[gb], Sort[(Exponent[#, vs] &) /@ gb],
     If[Length[uni] == 1,
        sig[Coefficient[First[uni], lv, 0]/Coefficient[First[uni], lv, Exponent[First[uni], lv]], 16],
        "none"]}]];
solChk[s_, cv_, d_] := If[ListQ[s] && Length[s] > 0 && MatchQ[s, {{__Rule} ..}],
  {Length[s], sig[Max[Re[cv /. s]], d]}, "unevaluated"];

Do[With[{lab = p[[1]], sys = p[[2]], vs = p[[3]], ord = p[[4]], cv = p[[5]]},
   run[lab <> " GroebnerBasis exact " <> If[ord === Lexicographic, "lex", "grevlex"],
       GroebnerBasis[sys, vs, MonomialOrder -> ord], gbChk[#, vs] &, 3];
   run[lab <> " NSolve machine", NSolve[sys, vs], solChk[#, cv, 8] &, 3];
   run[lab <> " NSolve wp40", NSolve[sys, vs, WorkingPrecision -> 40], solChk[#, cv, 18] &]],
  {p, problems}];
