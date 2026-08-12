(* ==========================================================================
   Experiment 12 -- Numerical quadrature (NIntegrate)
   ==========================================================================
   WHAT IT MEASURES.  src/numerical_calculus/: Gauss-Kronrod adaptive
   (gkadapt.c), double-exponential (dequad.c), oscillatory (oscint.c) and
   multidimensional cubature (cubature.c).

   FOUR INTEGRAND CLASSES, because one adaptive rule cannot be good at all four
   and a single smooth row would hide that: smooth, endpoint-singular,
   oscillatory, and 2-D.

   Checks are the integral VALUE rounded -- these have known answers, so the
   check tests accuracy as well as agreement.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"NIntegrate", "Integrate", "AccuracyGoal", "PrecisionGoal",
         "WorkingPrecision", "Method"}];

(* 1. Smooth: the baseline any rule handles. Exact answer 2. *)
bench["NIntegrate Sin[x] on [0,Pi]", NIntegrate[Sin[x], {x, 0, Pi}];];
check["NIntegrate Sin[x] on [0,Pi]", Round[10^6 NIntegrate[Sin[x], {x, 0, Pi}]]];

(* 2. Endpoint singularity: 1/Sqrt[x] on [0,1], exact answer 2.  Adaptive
   bisection alone converges slowly here; a double-exponential rule does not. *)
bench["NIntegrate 1/Sqrt[x] on [0,1]", NIntegrate[1/Sqrt[x], {x, 0, 1}];];
check["NIntegrate 1/Sqrt[x] on [0,1]",
  Round[10^4 NIntegrate[1/Sqrt[x], {x, 0, 1}]]];

(* 3. Oscillatory: 40 half-periods.  Cost scales with the number of zeros. *)
bench["NIntegrate Sin[40 x] on [0,Pi]",
  NIntegrate[Sin[40 x], {x, 0, Pi}];];
check["NIntegrate Sin[40 x] on [0,Pi]",
  Round[10^4 NIntegrate[Sin[40 x], {x, 0, Pi}]]];

(* 4. A peaked integrand: a narrow Gaussian the sampler must find. *)
bench["NIntegrate narrow Gaussian",
  NIntegrate[Exp[-1000 (x - 1/2)^2], {x, 0, 1}];];
check["NIntegrate narrow Gaussian",
  Round[10^6 NIntegrate[Exp[-1000 (x - 1/2)^2], {x, 0, 1}]]];

(* 5. Two dimensions: cubature rather than a product rule. *)
bench["NIntegrate 2-D smooth",
  NIntegrate[Sin[x] Cos[y], {x, 0, Pi}, {y, 0, Pi/2}];];
check["NIntegrate 2-D smooth",
  Round[10^4 NIntegrate[Sin[x] Cos[y], {x, 0, Pi}, {y, 0, Pi/2}]]];

(* 6. An infinite interval: exact answer Sqrt[Pi]/2. *)
bench["NIntegrate Exp[-x^2] on [0,Infinity]",
  NIntegrate[Exp[-x^2], {x, 0, Infinity}];];
check["NIntegrate Exp[-x^2] on [0,Infinity]",
  Round[10^4 NIntegrate[Exp[-x^2], {x, 0, Infinity}]]];
