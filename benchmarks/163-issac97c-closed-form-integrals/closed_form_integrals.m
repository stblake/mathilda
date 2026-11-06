(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 163 -- ISSAC '97 holdout C (SEALED), family c8:
   definite integrals with closed forms (counterpart of Problem 8)
   ==========================================================================

       cd benchmarks/163-issac97c-closed-form-integrals
       ../../Mathilda -file closed_form_integrals.m

   Ten integrals whose values are Clausen, Glaisher, Apery-type (Zeta[3],
   Dirichlet beta(3)), Binet, Vardi, Weber and Beta-derivative constants.  None
   is a textbook table entry in the form given: each needs a change of
   variables, a series interchange or a special-function identity.

   THREE ROWS PER INTEGRAL: exact Integrate (checked through N of the result,
   20 digits), NIntegrate machine (7 digits: several have log or log-log
   endpoint singularities) and NIntegrate with WorkingPrecision -> 40 (20
   digits).

   TRUE ANSWERS (A = the closed form, derived as stated, evaluated in mpmath at
   50 digits; B = mpmath tanh-sinh quadrature at 45 digits, quadosc for
   c8.5).  A and B agree to >= 1e-31 (>= 1e-41 for all but c8.9).

   c8.1  int_0^(Pi/3) Log[2 Sin[x/2]] dx = -Cl_2(Pi/3) (the Clausen function is
         -int_0^t Log[2 Sin[x/2]]) = -1.01494160640965362502120255427.
   c8.2  int_0^(1/2) LogGamma[x] dx = (5/24) Log[2] + (3/2) Log[Glaisher]
         + (1/4) Log[Pi] = 0.803719849629681710151993020173.
   c8.3  int_0^(Pi/2) x Log[Sin[x]] dx = 7 Zeta[3]/16 - (Pi^2/8) Log[2]
         (Fourier series Log[2 Sin x] = -sum Cos[2kx]/k, integrated termwise)
         = -0.329236162849817068243549448583.
   c8.4  int_0^oo ArcTan[t]/(E^(2 Pi t) - 1) dt = 1/2 - Log[2 Pi]/4 (Binet's
         second formula for LogGamma at z = 1) = 0.0405307333976636291098351317972.
   c8.5  int_0^oo BesselJ[0, x] Sin[x/2]/x dx = ArcSin[1/2] = Pi/6 (the Weber
         integral int J0(x) Sin[a x]/x = ArcSin[a], a <= 1)
         = 0.523598775598298873077107230547.
   c8.6  int_(Pi/4)^(Pi/2) Log[Log[Tan[x]]] dx
         = (Pi/2) Log[Gamma[3/4] Sqrt[2 Pi]/Gamma[1/4]] (Vardi 1988)
         = -0.260442806300988445540093868790.
   c8.7  int_0^1 x LogGamma[x] dx = Log[2 Pi]/4 - Log[Glaisher]
         = 0.210714789568552108342911874627 (A: quadrature; B: PSLQ on
         {v, 1, Log 2, Log Pi, Log A, Zeta[3]/Pi^2} at 40 digits returns
         4 v - Log 2 - Log Pi + 4 Log A = 0).
   c8.8  int_0^1 Log[x] Log[1-x]^2/x dx = -Pi^4/180 (expand Log[1-x]^2 =
         2 sum H_(n-1) x^n/n and use sum H_(n-1)/n^3 = Zeta[4]/4)
         = -0.541161616855569095758001848271.
   c8.9  int_0^1 x^(-1/3) (1-x)^(-2/3) Log[x] dx = d/da B(a, 1/3) at a = 2/3
         = B(2/3, 1/3) (PolyGamma[2/3] - PolyGamma[1]) = Pi^2/3 - Sqrt[3] Pi Log[3]
         = -2.68811867848189624972422298605.
   c8.10 int_0^1 ArcTan[x] Log[x]/x dx = -beta(3) = -Pi^3/32 (termwise:
         int_0^1 x^(2k) Log[x] = -1/(2k+1)^2) = -0.968946146259369380483634845847.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Integrate", "NIntegrate", "LogGamma", "BesselJ"}];

(* The value is real.  An exact result may be written with complex
   intermediates (PolyLog at (-1)^(1/3)) and a machine quadrature of
   Log[Log[Tan[x]]] may pick up a rounding-level imaginary part next to Pi/4;
   either is accepted when |Im| <= 10^-d |value|, and checked on the real part. *)
rsig[d_][v_] := If[!NumericQ[v], "unevaluated",
  Module[{w = If[Precision[v] === Infinity, N[v, d + 10], v]},
    If[Abs[Im[w]] <= 10^-d Abs[w], sig[Re[w], d], "complex"]]];

ig[1] = {Log[2 Sin[x/2]], {x, 0, Pi/3}};
ig[2] = {LogGamma[x], {x, 0, 1/2}};
ig[3] = {x Log[Sin[x]], {x, 0, Pi/2}};
ig[4] = {ArcTan[x]/(Exp[2 Pi x] - 1), {x, 0, Infinity}};
ig[5] = {BesselJ[0, x] Sin[x/2]/x, {x, 0, Infinity}};
ig[6] = {Log[Log[Tan[x]]], {x, Pi/4, Pi/2}};
ig[7] = {x LogGamma[x], {x, 0, 1}};
ig[8] = {Log[x] Log[1 - x]^2/x, {x, 0, 1}};
ig[9] = {x^(-1/3) (1 - x)^(-2/3) Log[x], {x, 0, 1}};
ig[10] = {ArcTan[x] Log[x]/x, {x, 0, 1}};

Do[
  With[{f = ig[k][[1]], r = ig[k][[2]], lab = "c8." <> ToString[k]},
    run[lab <> " Integrate exact", Integrate[f, r], rsig[20]];
    run[lab <> " NIntegrate machine", NIntegrate[f, r], rsig[7], 3];
    run[lab <> " NIntegrate wp40", NIntegrate[f, r, WorkingPrecision -> 40, PrecisionGoal -> 20],
        rsig[20]]],
  {k, 10}];
