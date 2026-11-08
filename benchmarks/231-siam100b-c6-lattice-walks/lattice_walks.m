(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 231 -- SIAM 100-Digit Challenge, HOLDOUT B for Problem 6
   ("Biased Random Walk").  Ten new problems that need the same capability:
   a 2-D (or 3-D) lattice Green's function integral with an integrable
   singularity at k = 0, alone or inside a root solve.

   For a walk with one-step characteristic function phi(k), the expected
   number of visits to a site x is
       G(x) = (2 pi)^-d Int_{[-pi,pi]^d} Re[ e^{-i k.x} / (1 - phi(k)) ] dk,
   the return probability is 1 - 1/G(0) and the probability of ever visiting
   x != 0 is G(x)/G(0).  With phi = c + i s (c, s real),
       Re 1/(1 - phi) = (1 - c)/((1 - c)^2 + s^2)        (re2 below).
   Each integrand is folded by its own symmetries onto [0,pi]^d (where a half
   is folded by k2 -> -k2 the two halves appear as a sum of two terms).
   (In 1-D the formula would need an extra point mass at k = 0; every walk
   here is 2-D or 3-D, where it does not.)

     bc6.1  Z^2, E 1/4+e, W 1/4-e, N = S = 1/4: return probability 0.6 (G = 5/2); e?
     bc6.2  Z^2, E = N = 1/4+e, W = S = 1/4-e: return probability 1/2; e?
     bc6.3  Z^2, N = S = 1/6, E 1/3+e, W 1/3-e: return probability 1/2; e?
     bc6.4  triangular lattice (steps +-e1, +-e2, +-(e1+e2), 1/6 each; +e1 gets
            +e, -e1 gets -e): return probability 1/2; e?
     bc6.5  Polya's constant: return probability of the unbiased simple cubic
            walk (Watson's triple integral, no root).
     bc6.6  king walk (8 neighbours, 1/8 each; E 1/8+e, W 1/8-e): return
            probability 1/2; e?
     bc6.7  the original walk at e = 1/10: return probability (no root).
     bc6.8  Z^3, E 1/6+e, W 1/6-e, N,S,U,D 1/6: return probability 1/4; e?
            (a root solve around a singular TRIPLE integral -- the hard one)
     bc6.9  the original walk at e = 1/20: probability of ever visiting the
            east neighbour (1,0) = G((1,0))/G(0) (ratio of two integrals).
     bc6.10 diagonal-step walk (NE 1/4+e, SW 1/4-e, NW = SE = 1/4):
            phi = cos k1 cos k2 + 2 i e sin(k1+k2), singular at BOTH (0,0)
            and (pi,pi); return probability 0.4; e?

   Tiers: "<id> machine" = FindRoot / NIntegrate at machine precision;
   "<id> wp50" = the same with WorkingPrecision -> 60.  The book's AGM closed
   form and every other specialised reduction are NOT used as cases (they
   are the references: siam100_refs/bc6_ref.py).  True values: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindRoot", "NIntegrate"}];

re2[c_, s_] := (1 - c)/((1 - c)^2 + s^2);

g1[e_?NumericQ, wp_] := NIntegrate[re2[(Cos[k1] + Cos[k2])/2, 2 e Sin[k1]],
   {k1, 0, Pi}, {k2, 0, Pi}, WorkingPrecision -> wp]/Pi^2;
g2[e_?NumericQ, wp_] := NIntegrate[With[{c = (Cos[k1] + Cos[k2])/2},
     re2[c, 2 e (Sin[k1] + Sin[k2])] + re2[c, 2 e (Sin[k1] - Sin[k2])]],
   {k1, 0, Pi}, {k2, 0, Pi}, WorkingPrecision -> wp]/(2 Pi^2);
g3[e_?NumericQ, wp_] := NIntegrate[re2[(2 Cos[k1] + Cos[k2])/3, 2 e Sin[k1]],
   {k1, 0, Pi}, {k2, 0, Pi}, WorkingPrecision -> wp]/Pi^2;
g4[e_?NumericQ, wp_] := NIntegrate[
     re2[(Cos[k1] + Cos[k2] + Cos[k1 + k2])/3, 2 e Sin[k1]] +
     re2[(Cos[k1] + Cos[k2] + Cos[k1 - k2])/3, 2 e Sin[k1]],
   {k1, 0, Pi}, {k2, 0, Pi}, WorkingPrecision -> wp]/(2 Pi^2);
g5[wp_] := NIntegrate[1/(1 - (Cos[k1] + Cos[k2] + Cos[k3])/3),
   {k1, 0, Pi}, {k2, 0, Pi}, {k3, 0, Pi}, WorkingPrecision -> wp]/Pi^3;
g6[e_?NumericQ, wp_] := NIntegrate[re2[(Cos[k1] + Cos[k2] + 2 Cos[k1] Cos[k2])/4, 2 e Sin[k1]],
   {k1, 0, Pi}, {k2, 0, Pi}, WorkingPrecision -> wp]/Pi^2;
g8[e_?NumericQ, wp_] := NIntegrate[re2[(Cos[k1] + Cos[k2] + Cos[k3])/3, 2 e Sin[k1]],
   {k1, 0, Pi}, {k2, 0, Pi}, {k3, 0, Pi}, WorkingPrecision -> wp]/Pi^3;
h9[e_, wp_] := NIntegrate[With[{c = (Cos[k1] + Cos[k2])/2, s = 2 e Sin[k1]},
     (Cos[k1] (1 - c) + s Sin[k1])/((1 - c)^2 + s^2)],
   {k1, 0, Pi}, {k2, 0, Pi}, WorkingPrecision -> wp]/Pi^2;
g10[e_?NumericQ, wp_] := NIntegrate[With[{c = Cos[k1] Cos[k2]},
     re2[c, 2 e Sin[k1 + k2]] + re2[c, 2 e Sin[k1 - k2]]],
   {k1, 0, Pi}, {k2, 0, Pi}, WorkingPrecision -> wp]/(2 Pi^2);

mp = MachinePrecision;
vrun["bc6.1 machine", e /. FindRoot[g1[e, mp] == 5/2, {e, 0.02, 0.04}], 10, 1, 90];
vrun["bc6.1 wp50", e /. FindRoot[g1[e, 60] == 5/2, {e, 2/100, 4/100}, WorkingPrecision -> 60], 48, 1, 75];
vrun["bc6.2 machine", e /. FindRoot[g2[e, mp] == 2, {e, 0.03, 0.06}], 10, 1, 90];
vrun["bc6.2 wp50", e /. FindRoot[g2[e, 60] == 2, {e, 3/100, 6/100}, WorkingPrecision -> 60], 48, 1, 75];
vrun["bc6.3 machine", e /. FindRoot[g3[e, mp] == 2, {e, 0.06, 0.1}], 10, 1, 90];
vrun["bc6.3 wp50", e /. FindRoot[g3[e, 60] == 2, {e, 6/100, 1/10}, WorkingPrecision -> 60], 48, 1, 75];
vrun["bc6.4 machine", e /. FindRoot[g4[e, mp] == 2, {e, 0.03, 0.06}], 10, 1, 90];
vrun["bc6.4 wp50", e /. FindRoot[g4[e, 60] == 2, {e, 3/100, 6/100}, WorkingPrecision -> 60], 48, 1, 75];
vrun["bc6.5 machine", 1 - 1/g5[mp], 10, 1, 90];
vrun["bc6.5 wp50", 1 - 1/g5[60], 48, 1, 75];
vrun["bc6.6 machine", e /. FindRoot[g6[e, mp] == 2, {e, 0.02, 0.05}], 10, 1, 90];
vrun["bc6.6 wp50", e /. FindRoot[g6[e, 60] == 2, {e, 2/100, 5/100}, WorkingPrecision -> 60], 48, 1, 75];
vrun["bc6.7 machine", 1 - 1/g1[0.1, mp], 10, 1, 90];
vrun["bc6.7 wp50", 1 - 1/g1[1/10, 60], 48, 1, 75];
vrun["bc6.8 machine", e /. FindRoot[g8[e, mp] == 4/3, {e, 0.05, 0.1}], 10, 1, 90];
vrun["bc6.8 wp50", e /. FindRoot[g8[e, 60] == 4/3, {e, 5/100, 1/10}, WorkingPrecision -> 60], 48, 1, 75];
vrun["bc6.9 machine", h9[0.05, mp]/g1[0.05, mp], 10, 1, 90];
vrun["bc6.9 wp50", h9[1/20, 60]/g1[1/20, 60], 48, 1, 75];
vrun["bc6.10 machine", e /. FindRoot[g10[e, mp] == 5/3, {e, 0.08, 0.13}], 10, 1, 90];
vrun["bc6.10 wp50", e /. FindRoot[g10[e, 60] == 5/3, {e, 8/100, 13/100}, WorkingPrecision -> 60], 48, 1, 75];
