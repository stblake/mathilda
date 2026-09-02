(* Experiment 16 -- Interpolation.
   MEASURES src/interp.c and interp_mpfr.c: building an InterpolatingFunction and
   evaluating it, for both `Interpolation` (points -> function) and its value-only
   companion `ListInterpolation` (values on a regular grid -> function).  BUILD
   AND EVALUATE ARE SEPARATE ROWS on purpose -- a system can be fast to construct
   and slow to evaluate, and the two land in different files.

   EVERY CHECK EVALUATES AT A GRID NODE.  A node is a data point, so every
   interpolation method returns the stored value there regardless of the scheme
   (Wolfram's local Hermite vs scipy's global not-a-knot spline give different
   values BETWEEN nodes -- comparing them there would be a spurious CHECK-FAIL).
   Timing still measures the real construction / evaluation over the whole grid. *)

Get["../harness.m"];
Get["../data.m"];

(* InterpolationOrder is an OPTION, not a head -- `Names["InterpolationOrder"]`
   is {} even though `Interpolation[..., InterpolationOrder -> 1]` works, so it
   does not belong in the absence register. *)
require[{"Interpolation", "ListInterpolation", "InterpolatingFunction", "Fit"}];

(* MACHINE FLOATS, NOT EXACT RATIONALS.  `i/100` in the Wolfram language is an
   exact Rational, so `Interpolation` over such points keeps Sin[1/10] SYMBOLIC
   while the numpy column does float64 -- two different problems, and the first
   run of this experiment reported 70000x because of it, not because of any
   defect.  N[] pins both columns to machine doubles. *)
pts = N[Table[{i/10, Sin[i/10]}, {i, 1, 2000}]];

bench["Interpolation build, 2000 knots", Interpolation[pts];];
check["Interpolation build, 2000 knots", Length[pts]];

ifn = Interpolation[pts];
bench["Interpolation evaluate, 20000 points",
  Table[ifn[N[j/1000]], {j, 1000, 20999}];];
check["Interpolation evaluate, 20000 points", Round[10^6 N[ifn[2.5]]]];

(* Linear (order-1) build.  Keyed on nothing -- `Interpolation` is already
   required above, and InterpolationOrder is a supported option; the earlier
   benchIf gate on the option NAME silently SKIP'd this row on every run. *)
bench["Interpolation order 1 (linear) build",
  Interpolation[pts, InterpolationOrder -> 1];];
check["Interpolation order 1 (linear) build",
  Round[10^6 N[Interpolation[pts, InterpolationOrder -> 1][2.5]]]];

(* A 2-D grid: separable interpolation, a different code path. *)
grid = N[Table[Sin[i/10] Cos[j/10], {i, 1, 60}, {j, 1, 60}]];
benchIf["ListInterpolation 2-D 60x60", "ListInterpolation",
  ListInterpolation[grid];];
checkIf["ListInterpolation 2-D 60x60", "ListInterpolation",
  Round[10^6 N[ListInterpolation[grid][5., 5.]]]];

(* Evaluating on a packed array rather than element by element: this is the row
   that shows whether the interpolant has a vectorised path at all. *)
xs = rand01[{100000}] * 100 + 1;
bench["Interpolation over 10^5 array", ifn[xs];];
check["Interpolation over 10^5 array", Round[10^6 N[ifn[7.5]]]];

(* ==========================================================================
   ListInterpolation -- the value-only companion (src/interp.c).  A raw array of
   values on a regular grid, abscissae synthesised, delegating to the same
   Interpolation engine (so the same vectorised InterpolatingFunction object and
   the buffer-direct construction path off a packed value tensor).  The scipy
   analog is a CubicSpline over an arange grid (1-D), or RectBivariateSpline for
   the tensor case -- both node-exact interpolating cubics, so a grid-node check
   agrees to machine precision.
   ========================================================================== *)

(* 1-D build from 10^5 values on the integer grid 1, 2, ... .  The value tensor
   is auto-packed, so this exercises the buffer-direct construction path. *)
lvals = N[Table[Sin[i/10], {i, 1, 100000}]];
benchIf["ListInterpolation 1-D build, 10^5 values", "ListInterpolation",
  ListInterpolation[lvals];];
checkIf["ListInterpolation 1-D build, 10^5 values", "ListInterpolation",
  Round[10^6 N[ListInterpolation[lvals][25]]]];       (* node 25 -> Sin[2.5] *)

(* 1-D evaluate: the vectorised application path over a 10^5 packed query array,
   the ListInterpolation analog of the "Interpolation over 10^5 array" row. *)
lifn = ListInterpolation[lvals];
lxs = rand01[{100000}] * 99000 + 1;                    (* in-domain [1, 99001] *)
benchIf["ListInterpolation 1-D evaluate, 10^5 array", "ListInterpolation",
  lifn[lxs];];
checkIf["ListInterpolation 1-D evaluate, 10^5 array", "ListInterpolation",
  Round[10^6 N[lifn[7500]]]];                          (* node 7500 -> Sin[750] *)

(* 2-D build: a 200x200 separable grid, the tensor construction path at scale. *)
lgrid = N[Table[Sin[i/10] Cos[j/10], {i, 1, 200}, {j, 1, 200}]];
benchIf["ListInterpolation 2-D build, 200x200", "ListInterpolation",
  ListInterpolation[lgrid];];
checkIf["ListInterpolation 2-D build, 200x200", "ListInterpolation",
  Round[10^6 N[ListInterpolation[lgrid][50, 80]]]];    (* node -> Sin[5] Cos[8] *)
