(* Experiment 31 -- NIntegrate, a difficulty ladder.
   ==========================================================================
   WHY A LADDER.  Experiment 12 measures four integrand classes at one
   difficulty each, which answers "is the quadrature there" but not "where does
   it break".  An adaptive rule does not degrade smoothly: it is flat until the
   integrand out-runs the rule's assumption, then falls off a cliff.  Every
   family below is therefore run at INCREASING difficulty along its own axis,
   so a cliff shows up as a jump between two adjacent rows of the same family
   rather than as one number with nothing to compare it to.

   THE AXES
     oscillation   Sin[k x], k = 40 -> 200 -> 1000.  Cost should grow like the
                   number of half-periods unless a dedicated oscillatory rule
                   is reached; a rule that merely bisects pays much more.
     peak width    Exp[-w (x-1/2)^2], w = 1e3 -> 1e5.  The peak occupies
                   1/Sqrt[w] of the interval, so a bisecting rule must find a
                   feature 10x narrower in the second row than the first.
     endpoint      x^-1/2 -> x^-9/10.  Both integrable, the second far closer
                   to divergence; this is what double-exponential is for.
     dimension     1-D -> 2-D -> 3-D.  Cubature cost is exponential in d, so
                   these rows say what the base is.
     precision     machine -> 30 digits -> 50 digits.  This is where a fixed
                   double-precision rule has to hand over to MPFR.

   CHECKS ARE ACCURACY, NOT AGREEMENT.  Every integrand here has a known
   closed-form value, and each check rounds the computed result against it. A
   system that is fast because it stopped early fails the check rather than
   winning the row -- the distinction the ranked report cannot otherwise make.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"NIntegrate", "WorkingPrecision", "AccuracyGoal", "PrecisionGoal"}];

(* ---- oscillation: Integrate[Sin[k x], {x, 0, Pi}] = (1 - Cos[k Pi])/k ----
   k even => exactly 0, so the check is the SCALED value: rounding a true zero
   to four places would pass for any answer under 5e-5 and measure nothing. *)
bench["NI oscillatory k=40", NIntegrate[Sin[40 x], {x, 0, Pi}];];
check["NI oscillatory k=40", Round[10^6 NIntegrate[Sin[40 x], {x, 0, Pi}]]];

bench["NI oscillatory k=200", NIntegrate[Sin[200 x], {x, 0, Pi}];];
check["NI oscillatory k=200", Round[10^6 NIntegrate[Sin[200 x], {x, 0, Pi}]]];

bench["NI oscillatory k=1000", NIntegrate[Sin[1000 x], {x, 0, Pi}];];
check["NI oscillatory k=1000", Round[10^6 NIntegrate[Sin[1000 x], {x, 0, Pi}]]];

(* Odd k, so the answer is 2/k rather than 0 -- a nonzero oscillatory value
   the check can actually pin down. *)
bench["NI oscillatory k=1001 nonzero", NIntegrate[Sin[1001 x], {x, 0, Pi}];];
check["NI oscillatory k=1001 nonzero",
  Round[10^6 NIntegrate[Sin[1001 x], {x, 0, Pi}]]];

(* ---- peak width: Integrate[Exp[-w (x-1/2)^2], {x,0,1}] -> Sqrt[Pi/w] ----
   The ladder stops at w = 1e5. At w = 1e6 Mathematica returns 0.000886, i.e.
   exactly HALF the correct 0.00177245, having resolved one side of the peak
   and not the other; Mathilda and scipy both get it right. That is a real
   finding, but as a row it would only ever be a check disagreement, and the
   harness withholds timings for those -- so it is recorded here rather than
   left as a permanently red row measuring nothing. *)
bench["NI narrow peak w=1e3", NIntegrate[Exp[-1000 (x - 1/2)^2], {x, 0, 1}];];
check["NI narrow peak w=1e3",
  Round[10^8 NIntegrate[Exp[-1000 (x - 1/2)^2], {x, 0, 1}]]];

bench["NI narrow peak w=1e5",
  NIntegrate[Exp[-100000 (x - 1/2)^2], {x, 0, 1}];];
check["NI narrow peak w=1e5",
  Round[10^8 NIntegrate[Exp[-100000 (x - 1/2)^2], {x, 0, 1}]]];

(* ---- endpoint singularity: 1/Sqrt[x] -> 2, x^(-9/10) -> 10 ---- *)
bench["NI endpoint x^-1/2", NIntegrate[1/Sqrt[x], {x, 0, 1}];];
check["NI endpoint x^-1/2", Round[10^6 NIntegrate[1/Sqrt[x], {x, 0, 1}]]];

bench["NI endpoint x^-9/10", NIntegrate[x^(-9/10), {x, 0, 1}];];
check["NI endpoint x^-9/10", Round[10^4 NIntegrate[x^(-9/10), {x, 0, 1}]]];

(* ---- interior kink: Integrate[Abs[x - 1/3], {x,0,1}] = 5/18 ----
   Not a singularity, but the derivative jumps, so an error estimate built on a
   smoothness assumption is wrong exactly at the panel containing 1/3. *)
bench["NI interior kink", NIntegrate[Abs[x - 1/3], {x, 0, 1}];];
check["NI interior kink", Round[10^6 NIntegrate[Abs[x - 1/3], {x, 0, 1}]]];

(* ---- semi-infinite range: Integrate[Exp[-x^2], {x,0,Infinity}] = Sqrt[Pi]/2 *)
bench["NI semi-infinite Gaussian", NIntegrate[Exp[-x^2], {x, 0, Infinity}];];
check["NI semi-infinite Gaussian",
  Round[10^6 NIntegrate[Exp[-x^2], {x, 0, Infinity}]]];

(* ---- dimension: 2-D and 3-D cubature, both with closed forms ----
   2-D: Integrate[x y, {x,0,1},{y,0,1}] = 1/4;  3-D: x y z -> 1/8. *)
bench["NI 2-D product", NIntegrate[x y, {x, 0, 1}, {y, 0, 1}];];
check["NI 2-D product", Round[10^6 NIntegrate[x y, {x, 0, 1}, {y, 0, 1}]]];

bench["NI 3-D product", NIntegrate[x y z, {x, 0, 1}, {y, 0, 1}, {z, 0, 1}];];
check["NI 3-D product",
  Round[10^6 NIntegrate[x y z, {x, 0, 1}, {y, 0, 1}, {z, 0, 1}]]];

(* 2-D with a ridge rather than a product: no separability to exploit. *)
bench["NI 2-D ridge", NIntegrate[Exp[-100 (x - y)^2], {x, 0, 1}, {y, 0, 1}];];
check["NI 2-D ridge",
  Round[10^6 NIntegrate[Exp[-100 (x - y)^2], {x, 0, 1}, {y, 0, 1}]]];

(* ---- precision: machine -> 30 -> 50 digits ----
   Integrate[Sin[x], {x,0,Pi}] = 2 exactly, so the check is the same value at
   every precision and any digit lost is the row's own fault. *)
bench["NI 30-digit Sin", NIntegrate[Sin[x], {x, 0, Pi}, WorkingPrecision -> 30];];
check["NI 30-digit Sin",
  Round[10^6 NIntegrate[Sin[x], {x, 0, Pi}, WorkingPrecision -> 30]]];

bench["NI 50-digit Gaussian",
  NIntegrate[Exp[-x^2], {x, 0, 3}, WorkingPrecision -> 50];];
check["NI 50-digit Gaussian",
  Round[10^6 NIntegrate[Exp[-x^2], {x, 0, 3}, WorkingPrecision -> 50]]];
