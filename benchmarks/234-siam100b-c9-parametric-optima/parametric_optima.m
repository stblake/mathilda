(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 234 -- SIAM 100-Digit Challenge, HOLDOUT B for Problem 9
   (optimal parameter of an oscillatory integral).  Ten new problems that
   need the same capability: the global argmax (or argmin) over an interval
   of a parameter-dependent integral evaluated by NIntegrate inside the
   optimiser (a ?NumericQ black box).  Most integrands oscillate infinitely
   often at an endpoint or on an infinite range.

     bc9.1  F(a) = Int_0^1 sin(a/x) dx;                          argmax, 0 <= a <= 4
     bc9.2  F(a) = (2 + sin 10a) Int_0^1 x^a sin(a/(1-x)) dx;   argmax, 0 <= a <= 5
     bc9.3  F(a) = (2 + cos 7a) Int_0^oo x^-a sin x dx;         argmax, 1/10 <= a <= 3/2
     bc9.4  F(a) = (1/pi) Int_0^pi cos(a t - 10 sin t) dt (Anger J_a(10)); argmax, 0 <= a <= 5
     bc9.5  F(a) = Int_0^oo exp(-x^2) sin(a x) dx;              argmax, 0 <= a <= 10
     bc9.6  F(a) = Int_0^oo cos(a x)/(1 + x^6) dx;              ARGMIN, 0 <= a <= 6
     bc9.7  F(x) = Int_0^1 cos(x t)/sqrt(1 - t^2) dt;           ARGMIN, 0 <= x <= 10
     bc9.8  F(a) = Int_0^1 sin(a x^2) dx;                       argmax, 0 <= a <= 20
     bc9.9  F(x) = (1/pi) Int_0^oo cos(t^3/3 + x t) dt (Airy Ai); argmax, -4 <= x <= 0
            (oscillation with growing frequency on an infinite range: hard)
     bc9.10 F(b) = Int_0^oo sin(x + b/x) dx;                    argmax, 1/20 <= b <= 3
            (infinitely oscillating at BOTH ends: hard)

   Tiers:
     "<id> machine"  NMaximize / NMinimize over the interval constraint -- the
                     question is the GLOBAL optimum on the interval (several of
                     these have other local optima inside it);
     "<id> wp50"     FindMaximum / FindMinimum from a nearby start at
                     WorkingPrecision 60 (a global search at 60 digits is out
                     of reach; this tier asks for the local refinement only).
   The special-function closed forms are NOT used as cases (they are one of
   the two reference methods: siam100_refs/bc9_ref.py).  True values:
   truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate", "NMaximize", "NMinimize", "FindMaximum", "FindMinimum"}];

f1[a_?NumericQ, p_] := NIntegrate[Sin[a/x], {x, 0, 1}, WorkingPrecision -> p];
f2[a_?NumericQ, p_] := (2 + Sin[10 a]) NIntegrate[x^a Sin[a/(1 - x)], {x, 0, 1}, WorkingPrecision -> p];
f3[a_?NumericQ, p_] := (2 + Cos[7 a]) NIntegrate[x^-a Sin[x], {x, 0, Infinity}, WorkingPrecision -> p];
f4[a_?NumericQ, p_] := NIntegrate[Cos[a t - 10 Sin[t]], {t, 0, Pi}, WorkingPrecision -> p]/Pi;
f5[a_?NumericQ, p_] := NIntegrate[Exp[-x^2] Sin[a x], {x, 0, Infinity}, WorkingPrecision -> p];
f6[a_?NumericQ, p_] := NIntegrate[Cos[a x]/(1 + x^6), {x, 0, Infinity}, WorkingPrecision -> p];
f7[x_?NumericQ, p_] := NIntegrate[Cos[x t]/Sqrt[1 - t^2], {t, 0, 1}, WorkingPrecision -> p];
f8[a_?NumericQ, p_] := NIntegrate[Sin[a x^2], {x, 0, 1}, WorkingPrecision -> p];
f9[x_?NumericQ, p_] := NIntegrate[Cos[t^3/3 + x t], {t, 0, Infinity}, WorkingPrecision -> p]/Pi;
f10[b_?NumericQ, p_] := NIntegrate[Sin[x + b/x], {x, 0, Infinity}, WorkingPrecision -> p];

mp = MachinePrecision;
vrun["bc9.1 machine", a /. Last[NMaximize[{f1[a, mp], 0 <= a <= 4}, a]], 10, 1, 85];
vrun["bc9.1 wp50", a /. Last[FindMaximum[f1[a, 60], {a, 6/10}, WorkingPrecision -> 60]], 48, 1, 85];
vrun["bc9.2 machine", a /. Last[NMaximize[{f2[a, mp], 0 <= a <= 5}, a]], 10, 1, 85];
vrun["bc9.2 wp50", a /. Last[FindMaximum[f2[a, 60], {a, 2/10}, WorkingPrecision -> 60]], 48, 1, 85];
vrun["bc9.3 machine", a /. Last[NMaximize[{f3[a, mp], 1/10 <= a <= 3/2}, a]], 10, 1, 85];
vrun["bc9.3 wp50", a /. Last[FindMaximum[f3[a, 60], {a, 93/100}, WorkingPrecision -> 60]], 48, 1, 85];
vrun["bc9.4 machine", a /. Last[NMaximize[{f4[a, mp], 0 <= a <= 5}, a]], 10, 1, 85];
vrun["bc9.4 wp50", a /. Last[FindMaximum[f4[a, 60], {a, 214/100}, WorkingPrecision -> 60]], 48, 1, 85];
vrun["bc9.5 machine", a /. Last[NMaximize[{f5[a, mp], 0 <= a <= 10}, a]], 10, 1, 85];
vrun["bc9.5 wp50", a /. Last[FindMaximum[f5[a, 60], {a, 185/100}, WorkingPrecision -> 60]], 48, 1, 85];
vrun["bc9.6 machine", a /. Last[NMinimize[{f6[a, mp], 0 <= a <= 6}, a]], 10, 1, 85];
vrun["bc9.6 wp50", a /. Last[FindMinimum[f6[a, 60], {a, 43/10}, WorkingPrecision -> 60]], 48, 1, 85];
vrun["bc9.7 machine", x /. Last[NMinimize[{f7[x, mp], 0 <= x <= 10}, x]], 10, 1, 85];
vrun["bc9.7 wp50", x /. Last[FindMinimum[f7[x, 60], {x, 383/100}, WorkingPrecision -> 60]], 48, 1, 85];
vrun["bc9.8 machine", a /. Last[NMaximize[{f8[a, mp], 0 <= a <= 20}, a]], 10, 1, 85];
vrun["bc9.8 wp50", a /. Last[FindMaximum[f8[a, 60], {a, 26/10}, WorkingPrecision -> 60]], 48, 1, 85];
vrun["bc9.9 machine", x /. Last[NMaximize[{f9[x, mp], -4 <= x <= 0}, x]], 10, 1, 85];
vrun["bc9.9 wp50", x /. Last[FindMaximum[f9[x, 60], {x, -1}, WorkingPrecision -> 60]], 48, 1, 85];
vrun["bc9.10 machine", b /. Last[NMaximize[{f10[b, mp], 1/20 <= b <= 3}, b]], 10, 1, 85];
vrun["bc9.10 wp50", b /. Last[FindMaximum[f10[b, 60], {b, 2/10}, WorkingPrecision -> 60]], 48, 1, 85];
