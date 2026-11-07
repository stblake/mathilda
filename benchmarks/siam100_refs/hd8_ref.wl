(* Holdout-A D8 family (experiment 213), method B (independent of hd8_ref.py's
   residue radial integrals):
   1-D rows:  exact Integrate over R, then N;
   2-D Sobolev rows: polar, NESTED 1-D NIntegrates (r inner, theta outer, both
     numerical) at WorkingPrecision 90, as d8_small_ref.wl;
   3-D Sobolev row c_2(R^3): the radial integral Int_0^inf t^(1/2)/(1+t+q t^2) dt
     by exact symbolic Integrate (as a function of q, 2/3 <= q <= 1), the two
     angles by a 2-D NIntegrate at WorkingPrecision 80;
   Bessel-potential and anisotropic rows: exact Integrate. *)
psum[vars_, m_] := Total[Times @@@ (vars^(2 #) & /@ Select[Tuples[Range[0, m], Length[vars]], Total[#] <= m &])];
c1[m_] := Sqrt[Integrate[1/psum[{x}, m], {x, -Infinity, Infinity}]/(2 Pi)];
Print["hd8.1 ", N[c1[1], 70]]; Print["hd8.2 ", N[c1[2], 70]];
Print["hd8.3 ", N[c1[5], 70]]; Print["hd8.4 ", N[c1[10], 70]];
c2d[m_] := Module[{rad, P = psum[{x, y}, m]},
  rad[th_?NumericQ] := NIntegrate[r/(P /. {x -> r Cos[th], y -> r Sin[th]}), {r, 0, 1, Infinity},
     WorkingPrecision -> 90, PrecisionGoal -> 75, MaxRecursion -> 40];
  Sqrt[4 NIntegrate[rad[th], {th, 0, Pi/4, Pi/2}, WorkingPrecision -> 90, PrecisionGoal -> 72]/(2 Pi)^2]];
Print["hd8.5 ", N[c2d[3], 70]];
Print["hd8.6 ", N[c2d[4], 70]];
J = Integrate[Sqrt[t]/(1 + t + q t^2), {t, 0, Infinity}, Assumptions -> 2/3 <= q <= 1];
Print["   radial J(q) = ", InputForm[J]];
Jf[q0_?NumericQ] := J /. q -> q0;
qq[th_, ph_] := (1 + (Sin[th] Cos[ph])^4 + (Sin[th] Sin[ph])^4 + Cos[th]^4)/2;
(* nested 1-D NIntegrates (phi inner, theta outer): a 2-D NIntegrate stalls (eincr) at ~18 digits *)
inner7[th_?NumericQ] := NIntegrate[Jf[qq[th, ph]], {ph, 0, Pi/4, Pi/2}, WorkingPrecision -> 80, PrecisionGoal -> 65, MaxRecursion -> 40];
c7 = Sqrt[8 NIntegrate[inner7[th]/2 Sin[th], {th, 0, Pi/4, Pi/2}, WorkingPrecision -> 80, PrecisionGoal -> 62, MaxRecursion -> 40]/(2 Pi)^3];
Print["hd8.7 ", N[c7, 65]];
Print["hd8.8 ", N[Sqrt[Integrate[2 Pi r (1 + r^2)^(-3/2), {r, 0, Infinity}]/(2 Pi)^2], 70]];
Print["hd8.9 ", N[Sqrt[Integrate[4 Pi r^2 (1 + r^2)^(-5/2), {r, 0, Infinity}]/(2 Pi)^3], 70]];
Print["hd8.10 ", N[Sqrt[Integrate[Integrate[1/(1 + x^2 + y^4), {x, -Infinity, Infinity}, Assumptions -> y \[Element] Reals], {y, -Infinity, Infinity}]/(2 Pi)^2], 70]];
