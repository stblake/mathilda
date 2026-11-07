(* bench-timeout: 400 *)
(* ==========================================================================
   Experiment 174 -- SIAM 100-Digit Challenge, Problem 9 (an optimal
   parameter of an oscillatory integral).

       I(a) = Int_0^2 (2 + sin 10a) x^a sin(a/(2 - x)) dx.
       For which a in [0, 5] is I(a) maximal?             = 0.7859336743...

   Idiomatic statement only: the integral as written, by NIntegrate, inside
   the optimiser (a ?NumericQ helper, so the optimiser sees a black box).
     c9.machine  NMaximize over the constraint 0 <= a <= 5 -- the question is a
                 GLOBAL maximum (2 + sin 10a makes I multimodal: there are
                 local maxima at a = 1.385, 2.0, ...).
     c9.wp50     FindMaximum from a = 4/5 at WorkingPrecision 60 (a global
                 search at 60 digits is out of reach for every system; the
                 50-digit tier asks for the local refinement only).
   The integrand oscillates infinitely often as x -> 2.  The book's Meijer-G
   closed form is NOT used as a case.  True value and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate", "NMaximize", "FindMaximum"}];

f[a_?NumericQ] := (2 + Sin[10 a]) NIntegrate[x^a Sin[a/(2 - x)], {x, 0, 2}];
f60[a_?NumericQ] := (2 + Sin[10 a]) NIntegrate[x^a Sin[a/(2 - x)], {x, 0, 2},
                                                WorkingPrecision -> 60];

vrun["c9.machine", a /. Last[NMaximize[{f[a], 0 <= a <= 5}, a]], 10];
vrun["c9.wp50", a /. Last[FindMaximum[f60[a], {a, 4/5}, WorkingPrecision -> 60]], 48];
