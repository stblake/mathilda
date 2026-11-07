(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 201 -- SIAM 100-Digit Challenge, HOLDOUT A for Problem 6
   ("Biased Random Walk"): ten new lattice-walk return problems.  As in
   experiment 171, the return probability is p = 1 - 1/G with G the expected
   number of visits to the origin, stated as the Brillouin-zone integral of
   Re 1/(1 - phi(k)) (phi = characteristic function of one step),
       Re 1/(1 - phi) = (1 - c)/((1 - c)^2 + s^2),  c = Re phi,  s = Im phi,
   by NIntegrate, and -- for "which eps" -- FindRoot on G(eps) == 1/(1 - p).
   The integrand is (nearly) singular at k = 0.  Folding: when phi is even in
   each k_i separately the domain is [0,pi]^d / pi^d; when the bias couples
   the axes only (k -> -k) is a symmetry and it is [0,pi] x [-pi,pi] / (2 pi^2).

     hc6.1  Z^2, E/W = 1/4 +- 1/10, N = S = 1/4: p.                (no root)
     hc6.2  Z^2, E/W = 1/4 +- eps: eps with p = 2/5 (G = 5/3).
     hc6.3  Z^2, E/W = 1/4 +- eps: eps with p = 3/5 (G = 5/2): smaller eps,
            a sharper corner.
     hc6.4  Z^2, E and N = 1/4 + eps, W and S = 1/4 - eps: p = 1/2.
     hc6.5  Z^2 anisotropic, E/W = 1/3 +- eps, N = S = 1/6: p = 1/2.
     hc6.6  Z^2, E/W = 1/4 +- eps: eps with p = 1/5 (G = 5/4): a strongly
            biased walk, a milder corner (the easy one).  (Not a 1-D walk: in
            Z^1 the integral of Re 1/(1 - phi) misses a delta mass at k = 0
            and is not G.)
     hc6.7  Z^3 simple cubic walk: Polya's return probability (3-D integral of
            1/(1 - c), singular at 0; no root).
     hc6.8  Z^3, E/W = 1/6 +- eps, others 1/6: eps with p = 3/10 (a 3-D
            integral inside FindRoot: the harder one).
     hc6.9  Z^2, E/W = 1/4 +- 1/20, N/S = 1/4 +- 1/10: p.          (no root)
     hc6.10 triangular lattice (steps +-e1, +-e2, +-(e1+e2)), +-e1 = 1/6 +- eps,
            the rest 1/6: eps with p = 1/2.
   Rows: "<id> machine" = NIntegrate (+ FindRoot) at machine precision;
   "<id> wp50" = both with WorkingPrecision -> 60.  Every case is capped at
   60 s.  The book's AGM closed form / the Bessel-Laplace representation are
   NOT used as cases (they are the references).  True values: truth.tsv
   (siam100_refs/hc6_ref.py).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindRoot", "NIntegrate"}];

(* Same request in every system: a machine row (no options) asks for 12
   digits, as the .py (dblquad/tplquad epsrel 1e-12) and .mac (quad_qags
   epsrel 1e-12) halves do; a wp50 row passes its own options unchanged. *)
goal6[] := PrecisionGoal -> 12;
goal6[o__] := Sequence[o];

(* square lattice, phi even in each k_i: [0,pi]^2 *)
gsq[cc_, ss_, o___] := NIntegrate[With[{c = cc[k1, k2], s = ss[k1, k2]}, (1 - c)/((1 - c)^2 + s^2)],
   {k1, 0, Pi}, {k2, 0, Pi}, Evaluate[goal6[o]]]/Pi^2;
(* bias coupling the axes: [0,pi] x [-pi,pi] *)
gsq2[cc_, ss_, o___] := NIntegrate[With[{c = cc[k1, k2], s = ss[k1, k2]}, (1 - c)/((1 - c)^2 + s^2)],
   {k1, 0, Pi}, {k2, -Pi, Pi}, Evaluate[goal6[o]]]/(2 Pi^2);

cSq[k1_, k2_] := (Cos[k1] + Cos[k2])/2;
gEW[e_?NumericQ, o___]    := gsq[cSq, Function[{k1, k2}, 2 e Sin[k1]], Evaluate[goal6[o]]];
gDiag[e_?NumericQ, o___]  := gsq2[cSq, Function[{k1, k2}, 2 e (Sin[k1] + Sin[k2])], Evaluate[goal6[o]]];
gAni[e_?NumericQ, o___]   := gsq[Function[{k1, k2}, 2 Cos[k1]/3 + Cos[k2]/3], Function[{k1, k2}, 2 e Sin[k1]], Evaluate[goal6[o]]];
gEW3[e_?NumericQ, o___]   := NIntegrate[With[{c = (Cos[k1] + Cos[k2] + Cos[k3])/3},
     (1 - c)/((1 - c)^2 + (2 e Sin[k1])^2)], {k1, 0, Pi}, {k2, 0, Pi}, {k3, 0, Pi}, Evaluate[goal6[o]]]/Pi^3;
gTri[e_?NumericQ, o___]   := gsq2[Function[{k1, k2}, (Cos[k1] + Cos[k2] + Cos[k1 + k2])/3],
     Function[{k1, k2}, 2 e Sin[k1]], Evaluate[goal6[o]]];
s9[k1_, k2_] := Sin[k1]/10 + Sin[k2]/5;
u3[o___] := NIntegrate[1/(1 - (Cos[k1] + Cos[k2] + Cos[k3])/3), {k1, 0, Pi}, {k2, 0, Pi}, {k3, 0, Pi}, Evaluate[goal6[o]]]/Pi^3;

vrun["hc6.1 machine", 1 - 1/gEW[1/10], 10, 1, 60];
vrun["hc6.1 wp50", 1 - 1/gEW[1/10, WorkingPrecision -> 60], 48, 1, 60];
vrun["hc6.2 machine", e /. FindRoot[gEW[e] == 5/3, {e, 0.10, 0.11}], 10, 1, 60];
vrun["hc6.2 wp50", e /. FindRoot[gEW[e, WorkingPrecision -> 60] == 5/3, {e, 10/100, 11/100},
                                 WorkingPrecision -> 60], 48, 1, 60];
vrun["hc6.3 machine", e /. FindRoot[gEW[e] == 5/2, {e, 0.027, 0.029}], 10, 1, 60];
vrun["hc6.3 wp50", e /. FindRoot[gEW[e, WorkingPrecision -> 60] == 5/2, {e, 27/1000, 29/1000},
                                 WorkingPrecision -> 60], 48, 1, 60];
vrun["hc6.4 machine", e /. FindRoot[gDiag[e] == 2, {e, 0.04, 0.05}], 10, 1, 60];
vrun["hc6.4 wp50", e /. FindRoot[gDiag[e, WorkingPrecision -> 60] == 2, {e, 4/100, 5/100},
                                 WorkingPrecision -> 60], 48, 1, 60];
vrun["hc6.5 machine", e /. FindRoot[gAni[e] == 2, {e, 0.075, 0.085}], 10, 1, 60];
vrun["hc6.5 wp50", e /. FindRoot[gAni[e, WorkingPrecision -> 60] == 2, {e, 75/1000, 85/1000},
                                 WorkingPrecision -> 60], 48, 1, 60];
vrun["hc6.6 machine", e /. FindRoot[gEW[e] == 5/4, {e, 0.2, 0.22}], 10, 1, 60];
vrun["hc6.6 wp50", e /. FindRoot[gEW[e, WorkingPrecision -> 60] == 5/4, {e, 20/100, 22/100},
                                 WorkingPrecision -> 60], 48, 1, 60];
vrun["hc6.7 machine", 1 - 1/u3[], 10, 1, 60];
vrun["hc6.7 wp50", 1 - 1/u3[WorkingPrecision -> 60], 48, 1, 60];
vrun["hc6.8 machine", e /. FindRoot[gEW3[e] == 10/7, {e, 0.03, 0.035}], 10, 1, 60];
vrun["hc6.8 wp50", e /. FindRoot[gEW3[e, WorkingPrecision -> 60] == 10/7, {e, 30/1000, 35/1000},
                                 WorkingPrecision -> 60], 48, 1, 60];
vrun["hc6.9 machine", 1 - 1/gsq2[cSq, s9], 10, 1, 60];
vrun["hc6.9 wp50", 1 - 1/gsq2[cSq, s9, WorkingPrecision -> 60], 48, 1, 60];
vrun["hc6.10 machine", e /. FindRoot[gTri[e] == 2, {e, 0.045, 0.05}], 10, 1, 60];
vrun["hc6.10 wp50", e /. FindRoot[gTri[e, WorkingPrecision -> 60] == 2, {e, 45/1000, 5/100},
                                  WorkingPrecision -> 60], 48, 1, 60];
