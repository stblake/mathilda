(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 200 -- SIAM 100-Digit Challenge, holdout A for Problem 5
   (family hc5: best uniform polynomial approximation, posed as C5 poses it
   in experiment 170 -- the DISCRETE minimax over sample points, as a
   NON-SMOOTH minimisation over the coefficients):

       E = min over c of  max_k | f(x_k) - (c_0 + c_1 x_k + ... + c_n x_k^n) |

     k   f            samples x_k                         n   coefficients
     1   sqrt(x)      k/128, k = 0..128                   1   real   (E = 1/8 exactly; easy)
     2   e^x          cos(k pi/64), k = 0..64             3   real
     3   |x|          -1 + k/64, k = 0..128               4   real   (kink; rational E)
     4   1/Gamma(x)   2k/100, k = 0..100                  3   real
     5   arctan(5x)   cos(k pi/100), k = 0..100           5   real   (6 unknowns; harder)
     6   log(1+x)     k/64, k = 0..64                     2   real
     7   1/Gamma(z)   e^(2 pi i k/128), k = 0..127        2   complex (6 real unknowns)
     8   Gamma(z+2)   e^(2 pi i k/128)                    2   complex
     9   log(2+z)     e^(2 pi i k/64), k = 0..63          2   complex
     10  1/Gamma(z)   e^(2 pi i k/128)                    4   complex (10 real unknowns; harder)

   Every row is NMinimize over the real (and, for k >= 7, imaginary) parts
   of the coefficients, with no symmetry reduction, at machine precision and
   with WorkingPrecision -> 60; scored against the discrete E itself.
   True values: truth.tsv; references siam100_refs/hc5_ref.{py,wl}.
   Python: differential_evolution on a box (fixed seed); mpmath has no
   minimiser (wp50 SKIP).  Maxima: fmin_cobyla from the zero polynomial
   (machine), no bigfloat minimiser (wp50 SKIP).  FriCAS: SKIP.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NMinimize", "Gamma"}];

prob[1] = {Sqrt[#] &, Table[k/128, {k, 0, 128}], 1, False};
prob[2] = {Exp[#] &, Table[Cos[k Pi/64], {k, 0, 64}], 3, False};
prob[3] = {Abs[#] &, Table[-1 + k/64, {k, 0, 128}], 4, False};
prob[4] = {1/Gamma[#] &, Table[2 k/100, {k, 0, 100}], 3, False};
prob[5] = {ArcTan[5 #] &, Table[Cos[k Pi/100], {k, 0, 100}], 5, False};
prob[6] = {Log[1 + #] &, Table[k/64, {k, 0, 64}], 2, False};
prob[7] = {1/Gamma[#] &, Table[Exp[2 Pi I k/128], {k, 0, 127}], 2, True};
prob[8] = {Gamma[# + 2] &, Table[Exp[2 Pi I k/128], {k, 0, 127}], 2, True};
prob[9] = {Log[2 + #] &, Table[Exp[2 Pi I k/64], {k, 0, 63}], 2, True};
prob[10] = {1/Gamma[#] &, Table[Exp[2 Pi I k/128], {k, 0, 127}], 4, True};

av = Table[Symbol["a" <> ToString[j]], {j, 0, 5}];
bv = Table[Symbol["b" <> ToString[j]], {j, 0, 5}];

(* the objective, its variables *)
obj[k_, prec_] := Module[{f, xs, n, cplx, xk, fk},
  {f, xs, n, cplx} = prob[k];
  xk = N[xs, prec]; fk = N[f /@ xs, prec];
  If[cplx,
    {Max[Abs[fk - (av[[1]] + I bv[[1]]) - Sum[(av[[j + 1]] + I bv[[j + 1]]) xk^j, {j, 1, n}]]],
     Join[av[[1 ;; n + 1]], bv[[1 ;; n + 1]]]},
    {Max[Abs[fk - av[[1]] - Sum[av[[j + 1]] xk^j, {j, 1, n}]]], av[[1 ;; n + 1]]}]];

Do[
  With[{om = obj[k, MachinePrecision], o60 = obj[k, 60]},
    vrun["hc5." <> ToString[k] <> " machine", First[NMinimize[om[[1]], om[[2]]]], 10, 1, 90];
    vrun["hc5." <> ToString[k] <> " wp50", First[NMinimize[o60[[1]], o60[[2]],
                                                         WorkingPrecision -> 60]], 48, 1, 90]],
  {k, 10}];
