(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 204 -- SIAM 100-Digit Challenge, HOLDOUT A for Problem 9 (an
   optimal parameter of a parametric integral): ten new argmax / argmin
   problems.  As in experiment 174 the integral is stated as written, by
   NIntegrate, inside the optimiser (a ?NumericQ helper: a black box).

     hc9.1  argmax_{a in [0,5]}  (2 + cos 10a) Int_0^2 x^a sin(a/(2-x)) dx
     hc9.2  argmax_{a in [0,5]}  (2 + sin 10a) Int_0^1 x^a sin(a/(1-x)) dx
     hc9.3  argmin_{a in [0,2]}  Int_0^inf x^a e^-x dx           (easy: Gamma's minimum - 1)
     hc9.4  argmax_{a in [0,4]}  (1/pi) Int_0^pi cos(t - a sin t) dt   (easy: J1's first maximum)
     hc9.5  argmax_{a in [0,4]}  (2 + sin 3a) Int_0^1 x^(-1/2) cos(a log x) dx
            (infinite oscillation and a singularity at x = 0)
     hc9.6  argmax_{b in [0,3]}  Int_0^b sin(x^2) dx             (parameter in the limit; = sqrt(pi))
     hc9.7  argmax_{a in [0,3]}  Int_0^inf x^a e^-x sin x dx      (an interior maximum just above the boundary values)
     hc9.8  argmin_{a in [0,10]} Int_0^inf e^-x cos(a sqrt x) dx
     hc9.9  argmax_{a in [0,5]}  (2 + sin 10a) Int_0^2 x^a cos(a/(2-x)) dx
     hc9.10 argmax_{a in [0,5]}  (2 + sin 10a) Int_0^3 x^a sin(a/(3-x)) dx
   Rows (as in experiment 174):
     "<id> machine"  NMaximize / NMinimize over the interval constraint (the
                     question is a GLOBAL extremum);
     "<id> wp50"     FindMaximum / FindMinimum from a start near the answer at
                     WorkingPrecision 60 (a global search at 60 digits is out
                     of reach; the 50-digit tier asks for the local refinement).
   Machine rows are capped at 120 s (as in experiment 174), wp50 rows at 60 s.
   (A WL NMaximize that overruns its cap can abort the whole vrun statement,
   which then prints no line: scored as no answer.)  Closed forms (Meijer G, Gamma, Bessel,
   Dawson) are NOT used as cases.  True values: truth.tsv (siam100_refs/hc9_ref.py).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate", "NMaximize", "NMinimize", "FindMaximum", "FindMinimum"}];

f1[a_?NumericQ, o___]  := (2 + Cos[10 a]) NIntegrate[x^a Sin[a/(2 - x)], {x, 0, 2}, o];
f2[a_?NumericQ, o___]  := (2 + Sin[10 a]) NIntegrate[x^a Sin[a/(1 - x)], {x, 0, 1}, o];
f3[a_?NumericQ, o___]  := NIntegrate[x^a Exp[-x], {x, 0, Infinity}, o];
f4[a_?NumericQ, o___]  := NIntegrate[Cos[t - a Sin[t]], {t, 0, Pi}, o]/Pi;
f5[a_?NumericQ, o___]  := (2 + Sin[3 a]) NIntegrate[x^(-1/2) Cos[a Log[x]], {x, 0, 1}, o];
f6[b_?NumericQ, o___]  := NIntegrate[Sin[x^2], {x, 0, b}, o];
f7[a_?NumericQ, o___]  := NIntegrate[x^a Exp[-x] Sin[x], {x, 0, Infinity}, o];
f8[a_?NumericQ, o___]  := NIntegrate[Exp[-x] Cos[a Sqrt[x]], {x, 0, Infinity}, o];
f9[a_?NumericQ, o___]  := (2 + Sin[10 a]) NIntegrate[x^a Cos[a/(2 - x)], {x, 0, 2}, o];
f10[a_?NumericQ, o___] := (2 + Sin[10 a]) NIntegrate[x^a Sin[a/(3 - x)], {x, 0, 3}, o];

vrun["hc9.1 machine", a /. Last[NMaximize[{f1[a], 0 <= a <= 5}, a]], 10, 1, 120];
vrun["hc9.1 wp50", a /. Last[FindMaximum[f1[a, WorkingPrecision -> 60], {a, 64/100}, WorkingPrecision -> 60]], 48, 1, 60];
vrun["hc9.2 machine", a /. Last[NMaximize[{f2[a], 0 <= a <= 5}, a]], 10, 1, 120];
vrun["hc9.2 wp50", a /. Last[FindMaximum[f2[a, WorkingPrecision -> 60], {a, 2/10}, WorkingPrecision -> 60]], 48, 1, 60];
vrun["hc9.3 machine", a /. Last[NMinimize[{f3[a], 0 <= a <= 2}, a]], 10, 1, 120];
vrun["hc9.3 wp50", a /. Last[FindMinimum[f3[a, WorkingPrecision -> 60], {a, 46/100}, WorkingPrecision -> 60]], 48, 1, 60];
vrun["hc9.4 machine", a /. Last[NMaximize[{f4[a], 0 <= a <= 4}, a]], 10, 1, 120];
vrun["hc9.4 wp50", a /. Last[FindMaximum[f4[a, WorkingPrecision -> 60], {a, 184/100}, WorkingPrecision -> 60]], 48, 1, 60];
vrun["hc9.5 machine", a /. Last[NMaximize[{f5[a], 0 <= a <= 4}, a]], 10, 1, 120];
vrun["hc9.5 wp50", a /. Last[FindMaximum[f5[a, WorkingPrecision -> 60], {a, 15/100}, WorkingPrecision -> 60]], 48, 1, 60];
vrun["hc9.6 machine", b /. Last[NMaximize[{f6[b], 0 <= b <= 3}, b]], 10, 1, 120];
vrun["hc9.6 wp50", b /. Last[FindMaximum[f6[b, WorkingPrecision -> 60], {b, 177/100}, WorkingPrecision -> 60]], 48, 1, 60];
vrun["hc9.7 machine", a /. Last[NMaximize[{f7[a], 0 <= a <= 3}, a]], 10, 1, 120];
vrun["hc9.7 wp50", a /. Last[FindMaximum[f7[a, WorkingPrecision -> 60], {a, 15/10}, WorkingPrecision -> 60]], 48, 1, 60];
vrun["hc9.8 machine", a /. Last[NMinimize[{f8[a], 0 <= a <= 10}, a]], 10, 1, 120];
vrun["hc9.8 wp50", a /. Last[FindMinimum[f8[a, WorkingPrecision -> 60], {a, 3}, WorkingPrecision -> 60]], 48, 1, 60];
vrun["hc9.9 machine", a /. Last[NMaximize[{f9[a], 0 <= a <= 5}, a]], 10, 1, 120];
vrun["hc9.9 wp50", a /. Last[FindMaximum[f9[a, WorkingPrecision -> 60], {a, 12/100}, WorkingPrecision -> 60]], 48, 1, 60];
vrun["hc9.10 machine", a /. Last[NMaximize[{f10[a], 0 <= a <= 5}, a]], 10, 1, 120];
vrun["hc9.10 wp50", a /. Last[FindMaximum[f10[a, WorkingPrecision -> 60], {a, 142/100}, WorkingPrecision -> 60]], 48, 1, 60];
