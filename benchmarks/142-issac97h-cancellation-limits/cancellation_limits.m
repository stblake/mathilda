(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 142 -- ISSAC '97 holdout, family 7: limits with massive cancellation
   ==========================================================================

       cd benchmarks/142-issac97h-cancellation-limits
       ../../Mathilda -file cancellation_limits.m

   HOLDOUT.  Ten problems probing what experiment 132 (Limit[f[g[x]]/g[f[x]],
   x -> 0] for twelve-term sums of nested Sin/Tan/Sinh/Tanh) probes -- limits
   whose leading terms cancel to high order -- on different compositions
   (the inverse functions ArcSin/ArcTan/ArcSinh/ArcTanh, deeper nestings)
   and on limits at Infinity where the cancellation is exponential.  See
   benchmarks/ISSAC97_HOLDOUT.md.

   Each problem has three cases:
     exact    Limit[e, x -> 0 or Infinity], checked exact and to 20 digits;
     machine  the expression at one point x_m in double precision (x_m = 1/10
              or 10), where the cancellation still leaves >= 6 digits.  For
              the compositions h7.2, h7.3, h7.10 a double cannot evaluate
              f[g[x]] at any x where it approximates the limit (the inner
              value is ~1e-9 and f of it cancels completely), so the machine
              case is the REDUCED quotient f[x]/g[x] at x = 1/10;
     wp       N[e /. x -> x_w, 30] at an exact point x_w deep in the
              cancellation (100 to 450 digits cancel), checked to 20 digits:
              the system must raise its internal precision on its own, with
              $MaxExtraPrecision raised to 1000 for these rows (the default
              cap of 50 extra digits cannot reach them in any system, which
              then rightly answers N::meprec rather than 30 digits).

   PROBLEMS AND TRUE VALUES.  Limits: exact, from the leading terms of the
   Taylor expansions (f ~ a x^p, g ~ b x^q gives f[g]/g[f] -> a^(1-q) b^(p-1))
   or, at Infinity, by asymptotic expansion (h7.4 and h7.5 are Gruntz's
   thesis examples), each confirmed by evaluating the expression in mpmath at
   points approaching the limit at 3000 and 6000 digits (the two precisions
   agree to all digits shown).  Point values: mpmath at 200/400 (x_m) and
   3000/6000 (x_w) digits.
     f2 = Tan[Sin[x]] - Sin[Tan[x]] ~ x^7/30,  g2 = Tan[Tanh[x]] - Tanh[Tan[x]] ~ 2 x^7/45;
     f3, g3: the twelve-term sums of experiment 132 with Sin, Tan, Sinh, Tanh
       replaced by ArcSin, ArcTan, ArcSinh, ArcTanh: f3 ~ -16 x^11/189,
       g3 ~ -20 x^11/189;
     f10 = Tan[Sin[Tan[Sin[x]]]] - Sin[Tan[Sin[Tan[x]]]] ~ x^7/15,
       g10 = Sinh[Tanh[x]] - Tanh[Sinh[x]] ~ x^7/30.
   problem                                          limit        x_m value                x_w value
   h7.1  (Sin[Tan]-Tan[Sin])/(ArcSin[ArcTan]-ArcTan[ArcSin])  1   1.016735900392242564   10^-8: 1.0000000000000001666667
   h7.2  f2[g2[x]]/g2[f2[x]]                 4096/729 = 5.618655692684499314   f2/g2 0.7586975708497706116   10^-6: 5.618655692684505410761
   h7.3  f3[g3[x]]/g3[f3[x]]          (5/4)^10 = 9.313225746154785156   f3/g3 0.7999597155712885935   10^-4: 9.313225746154790314691
   h7.4  E^x (E^(1/x - E^-x) - E^(1/x)), x -> oo          -1   -1.105145831114269735   1000: -1.001000500166708341668
   h7.5  (E^(x E^-x/(E^-x + E^(-2x^2/(x+1)))) - E^x)/x  -E^2   -6.150321866214310661   1000: -7.374307488864213536807
   h7.6  (1 + 1/x)^(x^2) E^-x, x -> oo            E^(-1/2)   0.6256388323119630930   10^30: 0.6065306597126334236038
   h7.7  (Tan[Tan[Tan]]-Sin[Sin[Sin]])/(Tan[Tan]-Sin[Sin])  3/2   1.503869110345547682   10^-12: 1.5000000000000000000000004
   h7.8  (Log[Cos[x]] + x^2/2 + x^4/12 + x^6/45)/x^8  -17/2520   -0.006767975352357744219   10^-10: -0.006746031746031746031768
   h7.9  (Sinh[Tan]-Tan[Sinh])/(ArcSinh[ArcTan]-ArcTan[ArcSinh])  1   1.051089324812166719   10^-8: 1.0000000000000005000000
   h7.10 f10[g10[x]]/g10[f10[x]]                 1/64   f10/g10 2.063705497271428363   10^-5: 0.01562499996571180559310
   (h7.6's x_w case evaluates the same function as E^(x^2 Log[1 + 1/x] - x):
   the power form at an exact x = 10^30 is an exact rational to the power 10^60.)
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Limit", "N"}];

f2[x_] = Tan[Sin[x]] - Sin[Tan[x]];  g2[x_] = Tan[Tanh[x]] - Tanh[Tan[x]];
f3[x_] = ArcTan[ArcTanh[ArcSin[x]]] + ArcTanh[ArcSin[ArcTan[x]]] + ArcSin[ArcTan[ArcTanh[x]]] - ArcTan[ArcSin[ArcTanh[x]]] -
  ArcSin[ArcTanh[ArcTan[x]]] - ArcTanh[ArcTan[ArcSin[x]]] - ArcTan[ArcSinh[ArcTanh[x]]] - ArcSinh[ArcTanh[ArcTan[x]]] -
  ArcTanh[ArcTan[ArcSinh[x]]] + ArcTan[ArcTanh[ArcSinh[x]]] + ArcTanh[ArcSinh[ArcTan[x]]] + ArcSinh[ArcTan[ArcTanh[x]]];
g3[x_] = ArcSinh[ArcTanh[ArcSin[x]]] + ArcTanh[ArcSin[ArcSinh[x]]] + ArcSin[ArcSinh[ArcTanh[x]]] - ArcSinh[ArcSin[ArcTanh[x]]] -
  ArcSin[ArcTanh[ArcSinh[x]]] - ArcTanh[ArcSinh[ArcSin[x]]] - ArcTan[ArcSinh[ArcSin[x]]] - ArcSinh[ArcSin[ArcTan[x]]] -
  ArcSin[ArcTan[ArcSinh[x]]] + ArcTan[ArcSin[ArcSinh[x]]] + ArcSin[ArcSinh[ArcTan[x]]] + ArcSinh[ArcTan[ArcSin[x]]];
f10[x_] = Tan[Sin[Tan[Sin[x]]]] - Sin[Tan[Sin[Tan[x]]]];  g10[x_] = Sinh[Tanh[x]] - Tanh[Sinh[x]];

(* {label, expression, limit point, machine expression (or the reduced
   quotient), x_m, wp expression, x_w} *)
problems = {
 {"h7.1", (Sin[Tan[x]] - Tan[Sin[x]])/(ArcSin[ArcTan[x]] - ArcTan[ArcSin[x]]), 0, Automatic, 1/10, Automatic, 10^-8},
 {"h7.2", f2[g2[x]]/g2[f2[x]], 0, f2[x]/g2[x], 1/10, Automatic, 10^-6},
 {"h7.3", f3[g3[x]]/g3[f3[x]], 0, f3[x]/g3[x], 1/10, Automatic, 10^-4},
 {"h7.4", Exp[x] (Exp[1/x - Exp[-x]] - Exp[1/x]), Infinity, Automatic, 10, Automatic, 1000},
 {"h7.5", (Exp[x Exp[-x]/(Exp[-x] + Exp[-2 x^2/(x + 1)])] - Exp[x])/x, Infinity, Automatic, 10, Automatic, 1000},
 {"h7.6", (1 + 1/x)^(x^2) Exp[-x], Infinity, Automatic, 10, Exp[x^2 Log[1 + 1/x] - x], 10^30},
 {"h7.7", (Tan[Tan[Tan[x]]] - Sin[Sin[Sin[x]]])/(Tan[Tan[x]] - Sin[Sin[x]]), 0, Automatic, 1/10, Automatic, 10^-12},
 {"h7.8", (Log[Cos[x]] + x^2/2 + x^4/12 + x^6/45)/x^8, 0, Automatic, 1/10, Automatic, 10^-10},
 {"h7.9", (Sinh[Tan[x]] - Tan[Sinh[x]])/(ArcSinh[ArcTan[x]] - ArcTan[ArcSinh[x]]), 0, Automatic, 1/10, Automatic, 10^-8},
 {"h7.10", f10[g10[x]]/g10[f10[x]], 0, f10[x]/g10[x], 1/10, Automatic, 10^-5}};

exactChk[v_] := Which[!NumericQ[v], "unevaluated", Precision[v] =!= Infinity, "inexact", True, sig[v, 20]];

Do[With[{lab = p[[1]], e = p[[2]], pt = p[[3]], em = If[p[[4]] === Automatic, p[[2]], p[[4]]],
         xm = N[p[[5]]], ew = If[p[[6]] === Automatic, p[[2]], p[[6]]], xw = p[[7]]},
   run[lab <> " Limit exact", Limit[e, x -> pt], exactChk];
   run[lab <> " point machine" <> If[p[[4]] === Automatic, "", " reduced f/g"], em /. x -> xm,
       sig[#, If[lab === "h7.3", 4, 5]] &, 3];
   run[lab <> " point N[.., 30]", Block[{$MaxExtraPrecision = 1000}, N[ew /. x -> xw, 30]], sig[#, 20] &, 3]],
  {p, problems}];
