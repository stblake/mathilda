(* Experiment 20 -- Curve fitting and regression.
   MEASURES src/fit.c and src/linalg/lstsq.c: the normal-equations / QR least
   squares path.  Separate from experiment 15 (optimisation) because linear least
   squares is a DIRECT solve, not an iteration -- different code, different
   failure mode.  Checks are the fitted VALUE at a point. *)

Get["../harness.m"];
Get["../data.m"];

require[{"Fit", "LeastSquares", "PseudoInverse", "Design", "LinearSolve"}];

(* Deterministic data, so all three systems fit the same points. *)
(* MACHINE FLOATS, NOT EXACT RATIONALS.  `i/100` in the Wolfram language is an
   exact Rational, so `Fit` over such points does EXACT rational linear algebra
   while the numpy column does float64 -- two different problems, and the first
   run of this experiment reported 239000x because of it, not because of any
   defect.  N[] pins both columns to machine doubles. *)
pts1 = N[Table[{i/100, 3 (i/100) + 2}, {i, 1, 5000}]];
bench["Fit linear, 5000 points", Fit[pts1, {1, t}, t];];
check["Fit linear, 5000 points", Round[10^4 (Fit[pts1, {1, t}, t] /. t -> 7)]];

pts2 = N[Table[{i/500, (i/500)^2 - (i/500) + 1}, {i, 1, 5000}]];
bench["Fit quadratic, 5000 points", Fit[pts2, {1, t, t^2}, t];];
check["Fit quadratic, 5000 points",
  Round[10^4 (Fit[pts2, {1, t, t^2}, t] /. t -> 3)]];

bench["Fit degree-8 polynomial, 5000 points",
  Fit[pts2, Table[t^k, {k, 0, 8}], t], 1];
check["Fit degree-8 polynomial, 5000 points",
  Round[10^2 (Fit[pts2, Table[t^k, {k, 0, 8}], t] /. t -> 3)]];

(* LeastSquares on an overdetermined dense system: the QR path directly. *)
m = rand01[{2000, 20}];
rhs = rand01[{2000}];
bench["LeastSquares 2000x20", LeastSquares[m, rhs];];
check["LeastSquares 2000x20",
  Round[10^4 Total[LeastSquares[N[spdIntMatrix[8]], N[Table[i, {i, 8}]]]]]];

(* Fitting a trig basis rather than a polynomial one. *)
bench["Fit trig basis, 5000 points",
  Fit[pts1, {1, Sin[t], Cos[t]}, t];];
check["Fit trig basis, 5000 points",
  Round[10^2 (Fit[pts1, {1, Sin[t], Cos[t]}, t] /. t -> 1)]];
