(* ==========================================================================
   Experiment 61 -- Least squares & regularized fitting (group E)
   ==========================================================================
   WHAT IT MEASURES.  LeastSquares (src/linalg/lstsq.c) on large overdetermined
   systems, and Fit (src/fit.c) as a polynomial least-squares and with Tikhonov
   (ridge) regularization.  Experiment 20 times a plain polynomial Fit; this adds
   the overdetermined-solve scaling and the regularized path.

   Baseline is compiled: numpy.linalg.lstsq / scipy.linalg.lstsq (LAPACK gelsd),
   numpy.polyfit, and a numpy closed-form ridge (X^T X + lambda I)^-1 X^T y --
   which Mathilda's FitRegularization -> {"Tikhonov", lambda} matches exactly
   (fit.c reduces ridge to ordinary least squares on the augmented system).

   FindFit (nonlinear model fitting) is ABSENT in Mathilda; the benchIf row
   records it as a feature gap while scipy.optimize.curve_fit measures it.

   TIMING on large random systems; CHECKS on the small deterministic
   Vandermonde system, so the three systems pin the same coefficients.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"LeastSquares", "Fit", "FitRegularization", "FindFit"}];

seed[];
a1 = rand01[{2000, 50}];   b1 = rand01[{2000}];
a2 = rand01[{8000, 200}];  b2 = rand01[{8000}];
fitBig = N[Table[{i/100, 2 + 3 (i/100) + (i/100)^2 + Sin[1. i]}, {i, 1, 2000}]];

(* deterministic check system *)
des = N[Table[x^j, {x, 1, 12}, {j, 0, 4}]];
bc = N[Table[2 + 3 x + x^2 + Mod[x, 3], {x, 1, 12}]];
dataC = N[Table[{x, 2 + 3 x + x^2 + Mod[x, 3]}, {x, 1, 12}]];

bench["LeastSquares overdetermined 2000x50", LeastSquares[a1, b1];];
check["LeastSquares overdetermined 2000x50", Round[10^6 Total[LeastSquares[des, bc]]]];

bench["LeastSquares overdetermined 8000x200", LeastSquares[a2, b2];, 1];
check["LeastSquares overdetermined 8000x200", Round[10^6 Total[LeastSquares[des, bc]]]];

bench["Fit polynomial degree 4", Fit[fitBig, {1, x, x^2, x^3, x^4}, x];];
check["Fit polynomial degree 4",
      Round[10^6 (Fit[dataC, {1, x, x^2, x^3, x^4}, x] /. x -> 7)]];

bench["Fit ridge (Tikhonov) degree 4",
      Fit[fitBig, {1, x, x^2, x^3, x^4}, x, FitRegularization -> {"Tikhonov", 1/2}];];
check["Fit ridge (Tikhonov) degree 4",
      Round[10^6 (Fit[dataC, {1, x, x^2, x^3, x^4}, x,
                      FitRegularization -> {"Tikhonov", 1/2}] /. x -> 7)]];

(* FindFit is ABSENT -> SKIP; scipy.optimize.curve_fit measures the Python side. *)
benchIf["FindFit exponential model", "FindFit",
        FindFit[dataC, a Exp[b x], {a, b}, x];];
checkIf["FindFit exponential model", "FindFit",
        Round[10^3 (a /. FindFit[dataC, a Exp[b x], {a, b}, x])]];
