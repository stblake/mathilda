(* bench-timeout: 1200 *)
(* ==========================================================================
   Experiment 139 -- ISSAC '97 holdout, family 4: exact coefficient extraction
   ==========================================================================

       cd benchmarks/139-issac97h-coefficients
       ../../Mathilda -file coefficients.m

   HOLDOUT.  Ten problems probing what experiment 129 (the middle coefficient
   of a degree-6000 product of powers) probes -- one exact coefficient of a
   huge polynomial or power series, or a residue at a high-order pole -- on
   different generating functions.  See benchmarks/ISSAC97_HOLDOUT.md.

   Each problem has three cases:
     exact    Coefficient[Expand[p], x, k], SeriesCoefficient or Residue;
     machine  the coefficient as a Cauchy contour integral, default
              NIntegrate, on a REDUCED instance whose integrand stays inside
              the double range (and loses at most ~8 digits to cancellation),
              checked to 5 digits (a machine NIntegrate's default goal);
     wpNN     the same contour integral on the FULL instance with the
              trapezoidal rule (exact for a periodic integrand once the rule
              has more points than the band of frequencies), at a working
              precision covering the cancellation max|integrand| / |answer|.
   The contour is |x| = r for a coefficient and |x - a| = rho for a residue:
       coefficient [x^k] f = (1/(2 Pi r^k)) Int_0^2Pi f(r E^(I t)) E^(-I k t) dt.

   PROBLEMS AND TRUE VALUES.  Every exact value was computed by pure integer
   arithmetic in Python (no CAS), by a second route where one exists, and the
   check is {leading 15 digits, exponent} plus numerator and denominator mod
   10^9:
   h4.1  [x^2000] (1 + x + x^2)^2000 = 1.909454965696130777e952 (953 digits)
         = Sum[Binomial[n, 2k] Binomial[2k, k]] (both routes)
         machine n = 200: 9.168205127076471545e93
   h4.2  [x^2250] (1 + x)^1500 (1 + x + x^2 + x^3)^1000 = 3.984929649324310828e1051
         = Sum[Binomial[1000, j] Binomial[2500, 2250 - 2j]] (convolution and closed sum)
         machine 150/100, k = 225: 7.169778242085354729e103
   h4.3  [x^2200] (x + 2)^1500 (x^3 - x + 1)^700 (x^2 + 3x + 1)^400
         = 2.837757356724861352e992 (convolution)
         machine 150/70/40, k = 220: 6.229683408098898651e97
   h4.4  [x^20000] 1/((1-x)(1-x^2)(1-x^3)(1-x^5)(1-x^7)(1-x^11))
         = 11585914270681804 (coin-change recurrence)
         machine k = 1000: 3876087157; wp contour k = 2000: 119680802813 (at
         k = 20000 the trapezoidal rule's convergence test is fooled by
         aliasing at every recursion depth)
   h4.5  Residue[1/((x - 1)^200 (x + 2)^300), {x, 1}] = -Binomial[498, 199]/3^499
         = -1.002978091496547155e-94 (closed form: [t^199] (3 + t)^-300)
         machine 20/30: -Binomial[48, 19]/3^49 = -4.823184391527019831e-11
   h4.6  [x^1000] (1 - x)^1000 (1 + x)^2000 = -6.656284598533242651e449
         (alternating closed sum and convolution; ~39 digits cancel on |x| = 1)
         machine 100/200, k = 100: 3.902905082656172828e41
   h4.7  [x^2525] Product[1 + x^k, {k, 100}] = 17310240059487250166 x 10^8
         = 1.731024005948725017e27 (partitions into distinct parts <= 100, DP)
         machine Product to 30, k = 232: 8679280
   h4.8  [x^1500] 1/(Sqrt[1 - 4x] (1 - x)) = Sum[Binomial[2k, k], {k, 0, 1500}]
         = 2.389556237309696426e901 (closed sum)
         machine k = 150: 1.251527767776154611e89
   h4.9  [x^401] Tan[x] = 2.889611317106376218e-79 (zigzag / boustrophedon
         numbers T_401/401!)          machine k = 41: 1.158764443279885220e-8
   h4.10 Residue[1/Sin[x]^201, {x, 0}] = Binomial[200, 100]/2^200
         = 0.05634847900925642225 (closed form: x = ArcSin[t] turns the
         residue into [t^200] (1 - t^2)^(-1/2))
         machine Sin^21: Binomial[20, 10]/2^20 = 0.176197052001953125
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Expand", "Coefficient", "SeriesCoefficient", "Residue", "NIntegrate"}];

exactChk[v_] := If[IntegerQ[v] || Head[v] === Rational,
  {sig[v, 15], Mod[Numerator[v], 10^9], Mod[Denominator[v], 10^9]}, "unevaluated"];

(* Cauchy integrals.  cm: default NIntegrate (machine); cw: trapezoidal rule
   at WorkingPrecision wp.  Real part: every answer here is real. *)
cm[f_, k_, r_] := Re[NIntegrate[Evaluate[(f /. x -> r E^(I t)) E^(-I k t)], {t, 0, 2 Pi}]]/(2 Pi r^k);
cw[f_, k_, r_, wp_] := Re[NIntegrate[Evaluate[(f /. x -> r E^(I t)) E^(-I k t)], {t, 0, 2 Pi},
     Method -> "Trapezoidal", WorkingPrecision -> wp, PrecisionGoal -> 16, MaxRecursion -> 14]]/(2 Pi r^k);
rm[f_, a_, rho_] := Re[NIntegrate[Evaluate[(f /. x -> a + rho E^(I t)) rho E^(I t)], {t, 0, 2 Pi}]]/(2 Pi);
rw[f_, a_, rho_, wp_] := Re[NIntegrate[Evaluate[(f /. x -> a + rho E^(I t)) rho E^(I t)], {t, 0, 2 Pi},
     Method -> "Trapezoidal", WorkingPrecision -> wp, PrecisionGoal -> 16, MaxRecursion -> 14]]/(2 Pi);

p1[n_] := (1 + x + x^2)^n;
p2[a_, b_] := (1 + x)^a (1 + x + x^2 + x^3)^b;
p3[a_, b_, c_] := (x + 2)^a (x^3 - x + 1)^b (x^2 + 3 x + 1)^c;
r4 = 1/((1 - x) (1 - x^2) (1 - x^3) (1 - x^5) (1 - x^7) (1 - x^11));
r5[a_, b_] := 1/((x - 1)^a (x + 2)^b);
p6[a_, b_] := (1 - x)^a (1 + x)^b;
p7[m_] := Product[1 + x^k, {k, m}];
s8 = 1/(Sqrt[1 - 4 x] (1 - x));

run["h4.1 exact Coefficient[Expand]", Coefficient[Expand[p1[2000]], x, 2000], exactChk];
run["h4.1 contour NIntegrate machine n=200", cm[p1[200], 200, 1], sig[#, 5] &, 3];
run["h4.1 contour Trapezoidal wp30", cw[p1[2000], 2000, 1, 30], sig[#, 14] &];

run["h4.2 exact Coefficient[Expand]", Coefficient[Expand[p2[1500, 1000]], x, 2250], exactChk];
run["h4.2 contour NIntegrate machine", cm[p2[150, 100], 225, 1], sig[#, 5] &, 3];
run["h4.2 contour Trapezoidal wp30", cw[p2[1500, 1000], 2250, 1, 30], sig[#, 14] &];

run["h4.3 exact Coefficient[Expand]", Coefficient[Expand[p3[1500, 700, 400]], x, 2200], exactChk];
run["h4.3 contour NIntegrate machine", cm[p3[150, 70, 40], 220, 1], sig[#, 5] &, 3];
run["h4.3 contour Trapezoidal wp350", cw[p3[1500, 700, 400], 2200, 1, 350], sig[#, 14] &];

run["h4.4 exact SeriesCoefficient k=20000", SeriesCoefficient[r4, {x, 0, 20000}], exactChk];
run["h4.4 contour NIntegrate machine k=1000", cm[r4, 1000, 99/100], sig[#, 5] &, 3];
run["h4.4 contour Trapezoidal wp40 k=2000", cw[r4, 2000, 99/100, 40], sig[#, 14] &];

run["h4.5 exact Residue", Residue[r5[200, 300], {x, 1}], exactChk];
run["h4.5 contour NIntegrate machine 20/30", rm[r5[20, 30], 1, 1], sig[#, 5] &, 3];
run["h4.5 contour Trapezoidal wp30", rw[r5[200, 300], 1, 1, 30], sig[#, 14] &];

run["h4.6 exact Coefficient[Expand]", Coefficient[Expand[p6[1000, 2000]], x, 1000], exactChk];
run["h4.6 contour NIntegrate machine", cm[p6[100, 200], 100, 1], sig[#, 5] &, 3];
run["h4.6 contour Trapezoidal wp70", cw[p6[1000, 2000], 1000, 1, 70], sig[#, 14] &];

run["h4.7 exact Coefficient[Expand]", Coefficient[Expand[p7[100]], x, 2525], exactChk];
run["h4.7 contour NIntegrate machine m=30", cm[p7[30], 232, 1], sig[#, 5] &, 3];
run["h4.7 contour Trapezoidal wp30", cw[p7[100], 2525, 1, 30], sig[#, 14] &];

run["h4.8 exact SeriesCoefficient k=1500", SeriesCoefficient[s8, {x, 0, 1500}], exactChk];
run["h4.8 contour NIntegrate machine k=150", cm[s8, 150, 24/100], sig[#, 5] &, 3];
run["h4.8 contour Trapezoidal wp60 k=1500", cw[s8, 1500, 24/100, 60], sig[#, 14] &];

run["h4.9 exact SeriesCoefficient Tan k=401", SeriesCoefficient[Tan[x], {x, 0, 401}], exactChk];
run["h4.9 contour NIntegrate machine k=41", cm[Tan[x], 41, 3/2], sig[#, 5] &, 3];
run["h4.9 contour Trapezoidal wp40 k=401", cw[Tan[x], 401, 3/2, 40], sig[#, 14] &];

run["h4.10 exact Residue", Residue[1/Sin[x]^201, {x, 0}], exactChk];
run["h4.10 contour NIntegrate machine Sin^21", rm[1/Sin[x]^21, 0, 1], sig[#, 5] &, 3];
run["h4.10 contour Trapezoidal wp50", rw[1/Sin[x]^201, 0, 1, 50], sig[#, 14] &];
