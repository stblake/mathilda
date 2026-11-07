(* Sanity/reduced rows of d8, method 2 (independent of d8_small_ref.py's
   closed-form radial integral): c_3(R^1) by exact Integrate; c_2(R^2) in polar
   coordinates as NESTED 1-D NIntegrates (r inner, theta outer, both numerical)
   at WorkingPrecision 90.  (A Cartesian or polar 2-D NIntegrate stalls at
   15-18 digits: NIntegrate::eincr.) *)
Print["c3R1 ", Sqrt[Integrate[1/(1 + x^2 + x^4 + x^6), {x, -Infinity, Infinity}]/(2 Pi)]];
w[x_, y_] := 1 + x^2 + y^2 + x^4 + x^2 y^2 + y^4;
rad[th_?NumericQ] := NIntegrate[r/w[r Cos[th], r Sin[th]], {r, 0, 1, Infinity},
   WorkingPrecision -> 90, PrecisionGoal -> 75];
c2 = Sqrt[4 NIntegrate[rad[th], {th, 0, Pi/2}, WorkingPrecision -> 90, PrecisionGoal -> 72]/(2 Pi)^2];
Print["c2R2 ", N[c2, 75]];
