(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 243 -- SIAM 100-Digit Challenge, SEALED HOLDOUT B for Appendix D
   Problem 8 (experiment 183, "a Sobolev constant").

   The capability: NIntegrate of a rational function over ALL of R^N
   (infinite ranges in every coordinate, algebraic decay), N = 1..4, as
   written -- no polar/spherical reduction (that is the reference method,
   siam100_refs/bd8_ref.py).  Ten new problems:

     c_m(R^N) = Sqrt[ (2 Pi)^-N Int_{R^N} 1/Sum_{|k|<=m} x^(2k) dx ]   (k a multi-index)
     bd8.1   c_2(R^1)                         (easy; = (2 Sqrt 3)^(-1/2))
     bd8.2   c_6(R^1)                         (degree-12 denominator)
     bd8.3   c_3(R^2)
     bd8.4   c_4(R^2)
     bd8.5   c_3(R^3)                         (3-D, 20 monomials)
     bd8.6   c_2(R^3)                         (HARDER: tail ~ r^-4 in 3-D, slow decay)
     bd8.7   Int_{R^2} 1/(1 + x^4 + y^4)                 (= Pi EllipticK[1/2])
     bd8.8   Int_{R^3} (1 + x^2 + y^2 + z^2)^-3          (= Pi^2/4)
     bd8.9   Int_{R^4} (1 + x^2 + y^2 + z^2 + w^2)^-3    (HARDER: 4-D; = Pi^2/2)
     bd8.10  Int_{R^2} 1/(1 + x^4 + x^2 y^2 + y^4)       (= Pi EllipticK[1/4])

   machine: default NIntegrate;  wp50: WorkingPrecision -> 60,
   PrecisionGoal -> 50 (expected to be out of reach in 3-D and 4-D).
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate"}];

p1[m_] := Sum[x^(2 i), {i, 0, m}];
p2[m_] := Sum[x^(2 i) y^(2 j), {i, 0, m}, {j, 0, m - i}];
p3[m_] := Sum[x^(2 i) y^(2 j) z^(2 k), {i, 0, m}, {j, 0, m - i}, {k, 0, m - i - j}];

cM[1, m_] := With[{p = p1[m]}, Sqrt[NIntegrate[1/p, {x, -Infinity, Infinity}]/(2 Pi)]];
cM[2, m_] := With[{p = p2[m]}, Sqrt[NIntegrate[1/p, {x, -Infinity, Infinity}, {y, -Infinity, Infinity}]/(2 Pi)^2]];
cM[3, m_] := With[{p = p3[m]}, Sqrt[NIntegrate[1/p, {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, {z, -Infinity, Infinity}]/(2 Pi)^3]];
cW[1, m_] := With[{p = p1[m]}, Sqrt[NIntegrate[1/p, {x, -Infinity, Infinity},
                  WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)]];
cW[2, m_] := With[{p = p2[m]}, Sqrt[NIntegrate[1/p, {x, -Infinity, Infinity}, {y, -Infinity, Infinity},
                  WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)^2]];
cW[3, m_] := With[{p = p3[m]}, Sqrt[NIntegrate[1/p, {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, {z, -Infinity, Infinity},
                  WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)^3]];

Do[With[{n = cs[[1]], m = cs[[2]], k = cs[[3]], rp = If[cs[[1]] == 3, 1, 3]},
  vrun["bd8." <> ToString[k] <> " machine", cM[n, m], 10, rp, 60];
  vrun["bd8." <> ToString[k] <> " wp50", cW[n, m], 48, 1, 120]],
  {cs, {{1, 2, 1}, {1, 6, 2}, {2, 3, 3}, {2, 4, 4}, {3, 3, 5}, {3, 2, 6}}}];

vrun["bd8.7 machine", NIntegrate[1/(1 + x^4 + y^4), {x, -Infinity, Infinity}, {y, -Infinity, Infinity}], 10, 3, 60];
vrun["bd8.7 wp50", NIntegrate[1/(1 + x^4 + y^4), {x, -Infinity, Infinity}, {y, -Infinity, Infinity},
                              WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 120];
vrun["bd8.8 machine", NIntegrate[(1 + x^2 + y^2 + z^2)^-3, {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, {z, -Infinity, Infinity}], 10, 1, 60];
vrun["bd8.8 wp50", NIntegrate[(1 + x^2 + y^2 + z^2)^-3, {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, {z, -Infinity, Infinity},
                              WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 120];
vrun["bd8.9 machine", NIntegrate[(1 + x^2 + y^2 + z^2 + w^2)^-3,
                                 {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, {z, -Infinity, Infinity}, {w, -Infinity, Infinity}], 10, 1, 60];
vrun["bd8.9 wp50", NIntegrate[(1 + x^2 + y^2 + z^2 + w^2)^-3,
                              {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, {z, -Infinity, Infinity}, {w, -Infinity, Infinity},
                              WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 120];
vrun["bd8.10 machine", NIntegrate[1/(1 + x^4 + x^2 y^2 + y^4), {x, -Infinity, Infinity}, {y, -Infinity, Infinity}], 10, 3, 60];
vrun["bd8.10 wp50", NIntegrate[1/(1 + x^4 + x^2 y^2 + y^4), {x, -Infinity, Infinity}, {y, -Infinity, Infinity},
                               WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 120];
