(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 212 -- SIAM 100-Digit Challenge HOLDOUT A, family D7
   (generalises experiment 182, "area of a region").

   The capability: the area of a planar region given by inequalities, stated
   (as in 182) in its 1-D form -- the length of each vertical / horizontal
   chord integrated across the region, with the ends of the region located by
   FindRoot where they are transcendental.  Ten instances, varied in kind:

     hd7.1   x^4 + y^4 < 1                       Int_-1^1 2 (1-x^4)^(1/4) dx (easy;
                                                 = 4 Gamma(5/4)^2/Gamma(3/2))
     hd7.2   x^2 < y < cos x                     ends +-a, cos a = a^2 (FindRoot)
     hd7.3   e^{-(y+x^3)^2} > y^2/16 + e^{sin y} one island [y1, 0], cube-root kink
     hd7.4   e^{-(y+sinh x)^2} > y^2/32 + e^{sin y}    chord = ArcSinh difference
     hd7.5   e^{-(y+x^5)^2}  > y^2/32 + e^{sin y}      fifth roots
     hd7.6   e^{-(y+x^3)^2/2} > y^2/32 + e^{sin y}     wider band, kink moves
     hd7.7   sin x + sin y > 1 in [0,pi]^2      Int_0^pi (pi - 2 ArcSin[1 - sin x]) dx
     hd7.8   e^x + e^y < 4, x + y > 0           Int (Log[4 - e^x] + x) dx between
                                                 the exact ends +-Log[2 + Sqrt[3]]
     hd7.9   e^x - 2 < y < 1 - x^2              two transcendental ends (FindRoot)
     hd7.10  e^{-(y+x^3)^2} > y^2/100 + e^{sin y}  THREE islands (harder)

   For the f > g rows, f = e^{-(y + p(x))^2/sc} > g(y) <=> |y + p(x)| < s(y),
   s = Sqrt[-sc Log[g(y)]] (needs g < 1), so each horizontal line meets the
   region in ONE x-interval of length
       w(y) = p^-1(-y + s) - p^-1(-y - s),
   integrated over the islands where g < 1 (FindRoot for their ends from
   starting values read off a plot; y = 0 is an exact end, g(0) = 1), with a
   breakpoint at the kink where -y - s = 0 (the root p^-1 of zero has a
   vertical tangent) -- exactly the 182 recipe.
   machine: default goals; wp50: WorkingPrecision -> 60 in every FindRoot and
   NIntegrate, PrecisionGoal -> 50.  Caps: 30 s per machine run (warm-up + 3
   timed), 60 s per wp50 run.  True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate", "FindRoot", "Surd"}];

(* machine (wp = MachinePrecision) or 60-digit spellings of FindRoot / NIntegrate *)
fr[eq_, v_, v0_, MachinePrecision] := v /. FindRoot[eq, {v, N[v0]}];
fr[eq_, v_, v0_, wp_] := v /. FindRoot[eq, {v, v0}, WorkingPrecision -> wp];
ni[f_, rng_, MachinePrecision] := NIntegrate[f, rng];
ni[f_, rng_, wp_] := NIntegrate[f, rng, WorkingPrecision -> wp, PrecisionGoal -> 50];

g16[y_] := y^2/16 + Exp[Sin[y]];
g32[y_] := y^2/32 + Exp[Sin[y]];
g100[y_] := y^2/100 + Exp[Sin[y]];
cbrt[t_] := Surd[t, 3]; root5[t_] := Surd[t, 5];
(* area of {Exp[-(y + p(x))^2/sc] > g(y)}, inv = p^-1;
   islands = {{start of lower end (or 0), start of upper end (or 0), {kink starts}}, ...} *)
fgArea[g_, inv_, sc_, islands_, wp_] := Total[Function[isl, Module[{a, b, ks},
    a = If[isl[[1]] === 0, 0, fr[g[y] == 1, y, isl[[1]], wp]];
    b = If[isl[[2]] === 0, 0, fr[g[y] == 1, y, isl[[2]], wp]];
    ks = fr[-sc Log[g[y]] == y^2, y, #, wp] & /@ isl[[3]];
    ni[inv[-y + Sqrt[-sc Log[g[y]]]] - inv[-y - Sqrt[-sc Log[g[y]]]], Join[{y, a}, ks, {b}], wp]]] /@ islands];

area[1, wp_] := ni[2 (1 - x^4)^(1/4), {x, -1, 1}, wp];
area[2, wp_] := With[{a = fr[Cos[x] == x^2, x, 8/10, wp]}, ni[Cos[x] - x^2, {x, -a, a}, wp]];
area[3, wp_] := fgArea[g16, cbrt, 1, {{-26/10, 0, {-8/10}}}, wp];
area[4, wp_] := fgArea[g32, ArcSinh, 1, {{-28/10, 0, {}}, {38/10, 445/100, {}}}, wp];
area[5, wp_] := fgArea[g32, root5, 1, {{-28/10, 0, {-8/10}}, {38/10, 445/100, {}}}, wp];
area[6, wp_] := fgArea[g32, cbrt, 2, {{-28/10, 0, {-13/10}}, {38/10, 445/100, {}}}, wp];
area[7, wp_] := ni[Pi - 2 ArcSin[1 - Sin[x]], {x, 0, Pi}, wp];
area[8, wp_] := ni[Log[4 - Exp[x]] + x, {x, -Log[2 + Sqrt[3]], Log[2 + Sqrt[3]]}, wp];
area[9, wp_] := With[{a = fr[3 - x^2 - Exp[x] == 0, x, -17/10, wp], b = fr[3 - x^2 - Exp[x] == 0, x, 8/10, wp]},
  ni[3 - x^2 - Exp[x], {x, a, b}, wp]];
area[10, wp_] := fgArea[g100, cbrt, 1, {{-794/100, -704/100, {}}, {-304/100, 0, {-86/100}},
                                        {325/100, 585/100, {}}}, wp];

vrun["hd7.1 machine", area[1, MachinePrecision], 10, 3, 30];
vrun["hd7.1 wp50", area[1, 60], 48, 1, 60];
vrun["hd7.2 machine", area[2, MachinePrecision], 10, 3, 30];
vrun["hd7.2 wp50", area[2, 60], 48, 1, 60];
vrun["hd7.3 machine", area[3, MachinePrecision], 10, 3, 30];
vrun["hd7.3 wp50", area[3, 60], 48, 1, 60];
vrun["hd7.4 machine", area[4, MachinePrecision], 10, 3, 30];
vrun["hd7.4 wp50", area[4, 60], 48, 1, 60];
vrun["hd7.5 machine", area[5, MachinePrecision], 10, 3, 30];
vrun["hd7.5 wp50", area[5, 60], 48, 1, 60];
vrun["hd7.6 machine", area[6, MachinePrecision], 10, 3, 30];
vrun["hd7.6 wp50", area[6, 60], 48, 1, 60];
vrun["hd7.7 machine", area[7, MachinePrecision], 10, 3, 30];
vrun["hd7.7 wp50", area[7, 60], 48, 1, 60];
vrun["hd7.8 machine", area[8, MachinePrecision], 10, 3, 30];
vrun["hd7.8 wp50", area[8, 60], 48, 1, 60];
vrun["hd7.9 machine", area[9, MachinePrecision], 10, 3, 30];
vrun["hd7.9 wp50", area[9, 60], 48, 1, 60];
vrun["hd7.10 machine", area[10, MachinePrecision], 10, 3, 30];
vrun["hd7.10 wp50", area[10, 60], 48, 1, 60];
