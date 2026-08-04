(* Experiment 16 -- Interpolation.
   MEASURES src/interp.c and interp_mpfr.c: building an InterpolatingFunction and
   evaluating it.  BUILD AND EVALUATE ARE SEPARATE ROWS on purpose -- a system can
   be fast to construct and slow to evaluate, and the two land in different files.
   Checks are the interpolant's VALUE at a point between knots. *)

Get["../harness.m"];
Get["../data.m"];

require[{"Interpolation", "InterpolationOrder", "ListInterpolation",
         "InterpolatingFunction", "Fit"}];

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

benchIf["Interpolation order 1 (linear) build", "InterpolationOrder",
  Interpolation[pts, InterpolationOrder -> 1];];
checkIf["Interpolation order 1 (linear) build", "InterpolationOrder",
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
