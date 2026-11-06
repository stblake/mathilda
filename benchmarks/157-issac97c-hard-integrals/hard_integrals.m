(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 157 -- ISSAC '97 holdout C (SEALED), family c2:
   hard definite integrals, numerically (counterpart of Problem 2)
   ==========================================================================

       cd benchmarks/157-issac97c-hard-integrals
       ../../Mathilda -file hard_integrals.m

   Ten integrals that defeat a default quadrature in ten different ways:
   slowly decaying Bessel products (conditionally convergent), an infinitely
   oscillating integrable singularity, a Stieltjes transform next to its branch
   cut, an x^(-0.999) endpoint singularity, an accelerating (Airy) phase, an
   unannounced interior inverse-square-root singularity, a log integrand whose
   O(1) values cancel to 1e-12, a 2-D corner singularity, and the Borwein
   integral whose value differs from Pi/2 in the 11th digit.  Every row is
   NIntegrate with default options (machine) or WorkingPrecision -> 40 (wp40);
   the integration range is given as stated, with NO hint of where the
   difficulty sits.

   TRUE ANSWERS (A = the closed form below, derived or quoted, evaluated in
   mpmath at 60 digits; B = an independent mpmath quadrature at 30 digits:
   quadosc, tanh-sinh split at the difficulty, or a Levin-accelerated sum of
   period integrals).  A and B agree to >= 25 digits except c2.10 (17 digits,
   which still resolves the 2.3e-11 deficit; A is the published rational).
   Digits checked: machine 6, wp40 15 (12 for c2.8, which cancels 12 digits).

   c2.1  int_0^oo x J0(x) J0(2x)^2 dx = 1/(2 Pi Delta), Delta = Sqrt[15]/4 the
         area of the (1,2,2) triangle (Sonine-Gegenbauer) = 2/(Pi Sqrt[15])
         = 0.164374518416399948150658425158.  Integrand ~ x^(-1/2).
   c2.2  int_0^1 Sin[1/x]/x dx = int_1^oo Sin[t]/t dt = Pi/2 - SinIntegral[1]
         = 0.624713256427713604289968377817.
   c2.3  int_-1^1 Sqrt[1-x^2]/(x - z) dx = -Pi (z - Sqrt[z-1] Sqrt[z+1]) at
         z = 1/2 + I delta (the semicircle's Stieltjes transform; the pole sits
         delta from the path).  machine delta=1e-3:
         -1.56898252904292563667809530742 + 2.71755987209473961110709894196 I;
         wp40 delta=1e-10:
         -1.57079632661351668280789990658 + 2.72069904603716751055632205414 I.
   c2.4  int_0^oo x^(a-1) E^-x/(x+1) dx = Gamma[a] E Gamma[1-a, 1], a = 1/1000
         = 998.828034469465151587773873665  (B subtracts the endpoint term
         1/a analytically before integrating).
   c2.5  int_0^oo Cos[t^3/3 - 2 t] dt = Pi AiryAi[-2]
         = 0.714421505810163771281836190476.
   c2.6  int_0^2 dx/Sqrt[Abs[x^2 - 2]] = Pi/2 + Log[1 + Sqrt[2]]
         = 2.45216991381443964446393101662  (singular at Sqrt[2], unannounced).
   c2.7  int_0^oo J0(x) J0(2x) dx = (1/Pi) EllipticK[1/4] (Weber-Schafheitlin,
         parameter m = k^2 = 1/4) = 0.536591003574682187526420853985.
   c2.8  int_0^Pi Log[1 - 2 r Cos[t] + r^2] dt = 2 Pi Log[r] for r > 1
         (Jensen).  machine r = 1 + 1e-4: 6.28287116885660785085983559428e-4;
         wp40 r = 1 + 1e-12: 6.28318530717644488427169906772e-12.
   c2.9  int_0^1 int_0^1 -Log[x y]/(1 - x y) dx dy = 2 Zeta[3]
         = 2.40411380631918857079947632302  (Beukers).
   c2.10 int_0^oo prod_{k=0}^{7} Sinc[x/(2k+1)] dx
         = Pi/2 - (6879714958723010531/935615849440640907310521750000) Pi
         = 1.57079632677179604650584089425  (Borwein & Borwein 2001; Pi/2 =
         1.5707963267948966...).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"NIntegrate", "BesselJ"}];

cpx[d_][v_] := If[NumericQ[v], {sig[Re[v], d], sig[Im[v], d]}, "unevaluated"];

run["c2.1 NIntegrate machine", NIntegrate[x BesselJ[0, x] BesselJ[0, 2 x]^2, {x, 0, Infinity}], sig[#, 6] &, 3];
run["c2.1 NIntegrate wp40", NIntegrate[x BesselJ[0, x] BesselJ[0, 2 x]^2, {x, 0, Infinity}, WorkingPrecision -> 40, PrecisionGoal -> 20], sig[#, 15] &];

run["c2.2 NIntegrate machine", NIntegrate[Sin[1/x]/x, {x, 0, 1}], sig[#, 6] &, 3];
run["c2.2 NIntegrate wp40", NIntegrate[Sin[1/x]/x, {x, 0, 1}, WorkingPrecision -> 40, PrecisionGoal -> 20], sig[#, 15] &];

run["c2.3 NIntegrate machine", NIntegrate[Sqrt[1 - x^2]/(x - (1/2 + I/1000)), {x, -1, 1}], cpx[6], 3];
run["c2.3 NIntegrate wp40", NIntegrate[Sqrt[1 - x^2]/(x - (1/2 + I/10^10)), {x, -1, 1}, WorkingPrecision -> 40, PrecisionGoal -> 20], cpx[15]];

run["c2.4 NIntegrate machine", NIntegrate[x^(1/1000 - 1) Exp[-x]/(x + 1), {x, 0, Infinity}], sig[#, 6] &, 3];
run["c2.4 NIntegrate wp40", NIntegrate[x^(1/1000 - 1) Exp[-x]/(x + 1), {x, 0, Infinity}, WorkingPrecision -> 40, PrecisionGoal -> 20], sig[#, 15] &];

run["c2.5 NIntegrate machine", NIntegrate[Cos[t^3/3 - 2 t], {t, 0, Infinity}], sig[#, 6] &, 3];
run["c2.5 NIntegrate wp40", NIntegrate[Cos[t^3/3 - 2 t], {t, 0, Infinity}, WorkingPrecision -> 40, PrecisionGoal -> 20], sig[#, 15] &];

run["c2.6 NIntegrate machine", NIntegrate[1/Sqrt[Abs[x^2 - 2]], {x, 0, 2}], sig[#, 6] &, 3];
run["c2.6 NIntegrate wp40", NIntegrate[1/Sqrt[Abs[x^2 - 2]], {x, 0, 2}, WorkingPrecision -> 40, PrecisionGoal -> 20], sig[#, 15] &];

run["c2.7 NIntegrate machine", NIntegrate[BesselJ[0, x] BesselJ[0, 2 x], {x, 0, Infinity}], sig[#, 6] &, 3];
run["c2.7 NIntegrate wp40", NIntegrate[BesselJ[0, x] BesselJ[0, 2 x], {x, 0, Infinity}, WorkingPrecision -> 40, PrecisionGoal -> 20], sig[#, 15] &];

run["c2.8 NIntegrate machine",
    NIntegrate[Log[1 - 2 (1 + 10^-4) Cos[t] + (1 + 10^-4)^2], {t, 0, Pi}], sig[#, 6] &, 3];
run["c2.8 NIntegrate wp40",
    NIntegrate[Log[1 - 2 (1 + 10^-12) Cos[t] + (1 + 10^-12)^2], {t, 0, Pi}, WorkingPrecision -> 40, PrecisionGoal -> 20], sig[#, 12] &];

run["c2.9 NIntegrate machine", NIntegrate[-Log[x y]/(1 - x y), {x, 0, 1}, {y, 0, 1}], sig[#, 6] &, 3];
run["c2.9 NIntegrate wp40", NIntegrate[-Log[x y]/(1 - x y), {x, 0, 1}, {y, 0, 1}, WorkingPrecision -> 40, PrecisionGoal -> 20], sig[#, 15] &];

borwein = Product[Sinc[x/(2 k + 1)], {k, 0, 7}];
run["c2.10 NIntegrate machine", NIntegrate[Evaluate[borwein], {x, 0, Infinity}], sig[#, 6] &, 3];
run["c2.10 NIntegrate wp40", NIntegrate[Evaluate[borwein], {x, 0, Infinity}, WorkingPrecision -> 40, PrecisionGoal -> 20], sig[#, 15] &];
