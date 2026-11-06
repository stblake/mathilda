(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 134 -- ISSAC '97 Systems Challenge, Problem 9
   ==========================================================================

       cd benchmarks/134-issac97-p9-integral-operator
       ../../Mathilda -file integral_operator.m

   PROBLEM.  The largest eigenvalue of the integral operator
       (A f)(x) = Integrate[K(x, y) f(y), {y, 0, 1}],
       K(x, y) = E^(x + y + x^2 + y^2 + x y + x^2 y^2).
   K is smooth and positive, so the spectrum decays geometrically and a
   Nystrom discretisation converges exponentially in the node count.

   METHODS (the notebook's two):
     m1  Nystrom: composite 32-point Gauss-Legendre on c = 2^i subintervals,
         eigenvalues of K(x_i, x_j) w_j, i increased until the top eigenvalue
         is stationary to 1e-14.  The notebook took its nodes from the
         obsolete NumericalMath`GaussianQuadrature` package; here they come
         from a Golub-Welsch helper (Eigensystem of the Jacobi matrix) that
         runs identically in both systems, and GaussianQuadratureWeights is
         `require`d so its absence is still registered.  Machine, 20 and 40
         digits.
     m2  power iteration of the operator through FunctionInterpolation of
         an inner NIntegrate (`require`d; Mathilda has no
         FunctionInterpolation, so the case is SKIPped there).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Eigensystem", "Eigenvalues", "GaussianQuadratureWeights", "FunctionInterpolation"}];
(* Captured BEFORE any definition below mentions the symbol: parsing a
   definition creates the symbol, after which Names[] would report it. *)
haveFI = haveQ["FunctionInterpolation"];
(* In a system where GaussianQuadratureWeights lives in a loadable package
   (the 1997 NumericalMath`GaussianQuadrature`, today
   NumericalDifferentialEquationAnalysis`), load it first; Mathilda has it built
   in, and its $VersionNumber (< 1) skips the load. *)
If[$VersionNumber >= 1, Quiet[Needs["NumericalDifferentialEquationAnalysis`"]]];
haveGQW = haveQ["GaussianQuadratureWeights"];

K[x_, y_] := E^(x + y + x^2 + y^2 + x y + x^2 y^2);

(* n-point Gauss-Legendre on (0, 1) at precision prec, by Golub-Welsch. *)
gaussLegendre01[n_, prec_] := Module[{J, ev, evec, b},
  b[k_] := k/Sqrt[4 k^2 - 1];
  J = N[Table[Which[j == i + 1, b[i], i == j + 1, b[j], True, 0], {i, n}, {j, n}], prec];
  {ev, evec} = Eigensystem[J];
  Transpose[Sort[Transpose[{(1 + ev)/2, evec[[All, 1]]^2}]]]];

nystrom[prec_] := Module[{xs, ws, eOld = 1, eNew, i = 0, c, pts, wts},
  {xs, ws} = gaussLegendre01[32, prec];
  While[c = 2^i;
    pts = Flatten[Table[xs + s, {s, 0, c - 1}]/c];
    wts = Flatten[Table[ws, {c}]/c];
    eNew = Max[Re[Eigenvalues[Outer[K, pts, pts] . DiagonalMatrix[wts]]]];
    Abs[(eNew - eOld)/eNew] > 10^(-14) && i < 4,
    i = i + 1; eOld = eNew];
  eNew];

(* Check the largest node, not the weight sum: the sum is 1 to working
   precision, so its leading digits sit on the 0.999.../1.000... rounding edge. *)
run["p9 m1 Gauss nodes Golub-Welsch wp20", gaussLegendre01[32, 20], sig[Max[#[[1]]], 15] &, 3];
(* The notebook's own call, from the GaussianQuadrature package. *)
If[haveGQW,
  run["p9 m1 Gauss nodes GaussianQuadratureWeights wp20",
    Transpose[GaussianQuadratureWeights[32, 0, 1, 20]], sig[Max[#[[1]]], 15] &, 3],
  Print["SKIP\tp9 m1 Gauss nodes GaussianQuadratureWeights wp20\tGaussianQuadratureWeights"]];
run["p9 m1 Nystrom machine", nystrom[MachinePrecision], sig[#, 12] &, 3];
run["p9 m1 Nystrom wp20", nystrom[20], sig[#, 14] &, 3];
run["p9 m1 Nystrom wp40", nystrom[40], sig[#, 14] &];

(* m2: the notebook's operator iteration, verbatim apart from symbol names. *)
setIQ[p_, io_, ag_, ip_] := SetOptions[FunctionInterpolation, InterpolationPoints -> p,
  InterpolationOrder -> io, AccuracyGoal -> ag, InterpolationPrecision -> ip];
m2[] := Module[{i = 0, a, b, fpl},
  Quiet[fpl = FixedPointList[Function[fn,
     setIQ[i = i + 1; 20 + 10 i, i + 3, i, 16];
     a = FunctionInterpolation[NIntegrate[Evaluate[SetPrecision[K[x, y] fn[[1]][y], 25]],
           {y, 0, 1}, PrecisionGoal -> i, WorkingPrecision -> 25], {x, 0, 1}];
     b = FunctionInterpolation[a[x]/a[1/2], {x, 0, 1}];
     {b, Plus @@ Table[a[x]/b[x], {x, 1/100, 1, 1/100}]/100}],
    {1 &, 1}, 12, SameTest -> (Abs[(#1[[2]] - #2[[2]])/#2[[2]]] < 10^(-14) &)]];
  fpl[[-1, -1]]];
If[haveFI,
  run["p9 m2 operator iteration FunctionInterpolation", m2[], sig[#, 8] &, 1, 600],
  Print["SKIP\tp9 m2 operator iteration FunctionInterpolation\tFunctionInterpolation"]];
