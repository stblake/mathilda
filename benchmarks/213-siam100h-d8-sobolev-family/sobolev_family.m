(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 213 -- SIAM 100-Digit Challenge HOLDOUT A, family D8
   (generalises experiment 183, "a Sobolev constant").

   The capability: an N-dimensional integral over ALL of R^N of the
   reciprocal of a positive polynomial (or a power of one), c^2 = (2 Pi)^-N
   Int_{R^N} 1/P(x) dx, answer c = Sqrt[c^2] -- at sizes where the truth is
   known to 55+ digits (the original's c_10(R^10) is known to 5).  Ten
   instances, varied in kind (k a multi-index, |k| = Sum k_j):

     hd8.1   c_1(R^1)   P = 1 + x^2                      = 1/Sqrt[2]   (easy)
     hd8.2   c_2(R^1)   P = 1 + x^2 + x^4                = (2 Sqrt[3])^(-1/2)
     hd8.3   c_5(R^1)   P = Sum_{k<=5} x^(2k)
     hd8.4   c_10(R^1)  P = Sum_{k<=10} x^(2k)           (degree 20)
     hd8.5   c_3(R^2)   P = Sum_{|k|<=3} x^(2k1) y^(2k2) (10 monomials)
     hd8.6   c_4(R^2)   P = Sum_{|k|<=4} ...             (15 monomials)
     hd8.7   c_2(R^3)   P = Sum_{|k|<=2} ... in 3-D      (harder: a 3-D integral)
     hd8.8   P = (1 + x^2 + y^2)^(3/2)       Bessel potential, = 1/Sqrt[2 Pi]
     hd8.9   P = (1 + x^2 + y^2 + z^2)^(5/2) Bessel potential, = 1/(Pi Sqrt[6])  (3-D)
     hd8.10  P = 1 + x^2 + y^4               anisotropic,
                                             c^2 = Gamma(1/4)^2/(8 Pi^(3/2))

   Idiomatic statement: the N-dimensional NIntegrate over R^N as written,
   default goals (machine) and WorkingPrecision -> 60, PrecisionGoal -> 50
   (wp50).  The polynomials are written out.  No polar/spherical reduction:
   that is the reference method (siam100_refs/hd8_ref.{py,wl}).  Caps: 30 s
   per machine run (warm-up + 3 timed), 60 s per wp50 run.
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate"}];

(* hd8.1  c_1(R^1) *)
vrun["hd8.1 machine", Sqrt[NIntegrate[1/(1 + x^2), {x, -Infinity, Infinity}]/(2 Pi)^1], 10, 3, 30];
vrun["hd8.1 wp50", Sqrt[NIntegrate[1/(1 + x^2), {x, -Infinity, Infinity},
    WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)^1], 48, 1, 60];
(* hd8.2  c_2(R^1) *)
vrun["hd8.2 machine", Sqrt[NIntegrate[1/(1 + x^2 + x^4), {x, -Infinity, Infinity}]/(2 Pi)^1], 10, 3, 30];
vrun["hd8.2 wp50", Sqrt[NIntegrate[1/(1 + x^2 + x^4), {x, -Infinity, Infinity},
    WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)^1], 48, 1, 60];
(* hd8.3  c_5(R^1) *)
vrun["hd8.3 machine", Sqrt[NIntegrate[1/(1 + x^2 + x^4 + x^6 + x^8 + x^10), {x, -Infinity, Infinity}]/(2 Pi)^1], 10, 3, 30];
vrun["hd8.3 wp50", Sqrt[NIntegrate[1/(1 + x^2 + x^4 + x^6 + x^8 + x^10), {x, -Infinity, Infinity},
    WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)^1], 48, 1, 60];
(* hd8.4  c_10(R^1) *)
vrun["hd8.4 machine", Sqrt[NIntegrate[1/(1 + x^2 + x^4 + x^6 + x^8 + x^10 + x^12 + x^14 + x^16 + x^18 + x^20), {x, -Infinity, Infinity}]/(2 Pi)^1], 10, 3, 30];
vrun["hd8.4 wp50", Sqrt[NIntegrate[1/(1 + x^2 + x^4 + x^6 + x^8 + x^10 + x^12 + x^14 + x^16 + x^18 + x^20), {x, -Infinity, Infinity},
    WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)^1], 48, 1, 60];
(* hd8.5  c_3(R^2) *)
vrun["hd8.5 machine", Sqrt[NIntegrate[1/(1 + x^2 + y^2 + x^4 + x^2*y^2 + y^4 + x^6 + x^4*y^2 + x^2*y^4 + y^6), {x, -Infinity, Infinity}, {y, -Infinity, Infinity}]/(2 Pi)^2], 10, 3, 30];
vrun["hd8.5 wp50", Sqrt[NIntegrate[1/(1 + x^2 + y^2 + x^4 + x^2*y^2 + y^4 + x^6 + x^4*y^2 + x^2*y^4 + y^6), {x, -Infinity, Infinity}, {y, -Infinity, Infinity},
    WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)^2], 48, 1, 60];
(* hd8.6  c_4(R^2) *)
vrun["hd8.6 machine", Sqrt[NIntegrate[1/(1 + x^2 + y^2 + x^4 + x^2*y^2 + y^4 + x^6 + x^4*y^2 + x^2*y^4 + y^6 + x^8 + x^6*y^2 + x^4*y^4 + x^2*y^6 + y^8), {x, -Infinity, Infinity}, {y, -Infinity, Infinity}]/(2 Pi)^2], 10, 3, 30];
vrun["hd8.6 wp50", Sqrt[NIntegrate[1/(1 + x^2 + y^2 + x^4 + x^2*y^2 + y^4 + x^6 + x^4*y^2 + x^2*y^4 + y^6 + x^8 + x^6*y^2 + x^4*y^4 + x^2*y^6 + y^8), {x, -Infinity, Infinity}, {y, -Infinity, Infinity},
    WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)^2], 48, 1, 60];
(* hd8.7  c_2(R^3) *)
vrun["hd8.7 machine", Sqrt[NIntegrate[1/(1 + x^2 + y^2 + z^2 + x^4 + x^2*y^2 + x^2*z^2 + y^4 + y^2*z^2 + z^4), {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, {z, -Infinity, Infinity}]/(2 Pi)^3], 10, 3, 30];
vrun["hd8.7 wp50", Sqrt[NIntegrate[1/(1 + x^2 + y^2 + z^2 + x^4 + x^2*y^2 + x^2*z^2 + y^4 + y^2*z^2 + z^4), {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, {z, -Infinity, Infinity},
    WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)^3], 48, 1, 60];
(* hd8.8  Bessel potential (1+|x|^2)^(-3/2) on R^2 *)
vrun["hd8.8 machine", Sqrt[NIntegrate[1/((1 + x^2 + y^2)^(3/2)), {x, -Infinity, Infinity}, {y, -Infinity, Infinity}]/(2 Pi)^2], 10, 3, 30];
vrun["hd8.8 wp50", Sqrt[NIntegrate[1/((1 + x^2 + y^2)^(3/2)), {x, -Infinity, Infinity}, {y, -Infinity, Infinity},
    WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)^2], 48, 1, 60];
(* hd8.9  Bessel potential (1+|x|^2)^(-5/2) on R^3 *)
vrun["hd8.9 machine", Sqrt[NIntegrate[1/((1 + x^2 + y^2 + z^2)^(5/2)), {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, {z, -Infinity, Infinity}]/(2 Pi)^3], 10, 3, 30];
vrun["hd8.9 wp50", Sqrt[NIntegrate[1/((1 + x^2 + y^2 + z^2)^(5/2)), {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, {z, -Infinity, Infinity},
    WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)^3], 48, 1, 60];
(* hd8.10  anisotropic 1/(1+x^2+y^4) on R^2 *)
vrun["hd8.10 machine", Sqrt[NIntegrate[1/(1 + x^2 + y^4), {x, -Infinity, Infinity}, {y, -Infinity, Infinity}]/(2 Pi)^2], 10, 3, 30];
vrun["hd8.10 wp50", Sqrt[NIntegrate[1/(1 + x^2 + y^4), {x, -Infinity, Infinity}, {y, -Infinity, Infinity},
    WorkingPrecision -> 60, PrecisionGoal -> 50]/(2 Pi)^2], 48, 1, 60];
