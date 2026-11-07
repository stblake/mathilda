(* Holdout A, family d15 (experiment 220), WL cross-checks of hd15_ref.py:
   hd15.2 by BesselJZero, hd15.4 by NIntegrate (GlobalAdaptive, WP 80). *)
Print["hd15.2 ", N[Sqrt[2 BesselJZero[-1/4, 1]], 70]];
Print["hd15.4 ", NIntegrate[1/(Exp[y] - y), {y, 0, 1, 10, Infinity}, WorkingPrecision -> 80, PrecisionGoal -> 70]];
