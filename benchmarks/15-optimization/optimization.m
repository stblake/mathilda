(* Experiment 15 -- Numerical optimisation.
   MEASURES src/findmin.c and src/fit.c: unconstrained minimisation and linear
   least squares.  Checks are the objective VALUE at the minimum, rounded --
   invariant under which descent path was taken. *)

Get["../harness.m"];
Get["../data.m"];

require[{"FindMinimum", "FindMaximum", "Minimize", "NMinimize", "Fit",
         "FindFit", "LeastSquares"}];

(* 1. A quadratic bowl: the easiest possible case, and the floor. *)
bench["FindMinimum quadratic", FindMinimum[(x - 3)^2 + 2, {x, 0}];];
check["FindMinimum quadratic",
  Round[10^6 First[FindMinimum[(x - 3)^2 + 2, {x, 0}]]]];

(* 2. Rosenbrock in 2-D: the classic curved-valley test. *)
bench["FindMinimum Rosenbrock 2-D",
  FindMinimum[(1 - x)^2 + 100 (y - x^2)^2, {{x, -1}, {y, 1}}];];
check["FindMinimum Rosenbrock 2-D",
  Round[10^4 First[FindMinimum[(1 - x)^2 + 100 (y - x^2)^2,
                               {{x, -1}, {y, 1}}]]]];

(* 3. A 1-D transcendental objective. *)
bench["FindMinimum x Sin[x] on [2,6]",
  FindMinimum[x Sin[x], {x, 4}];];
check["FindMinimum x Sin[x] on [2,6]",
  Round[10^4 First[FindMinimum[x Sin[x], {x, 4}]]]];

(* 4. Linear least squares through Fit: a normal-equations solve. *)
(* MACHINE FLOATS, NOT EXACT RATIONALS.  `i/100` in the Wolfram language is an
   exact Rational, so `Fit` over such points does EXACT rational linear algebra
   while the numpy column does float64 -- two different problems, and the first
   run of this experiment reported 239000x because of it, not because of any
   defect.  N[] pins both columns to machine doubles. *)
pts = N[Table[{i, 2 i + 1}, {i, 1, 200}]];
bench["Fit linear, 200 points", Fit[pts, {1, t}, t];];
check["Fit linear, 200 points", Round[10^6 (Fit[pts, {1, t}, t] /. t -> 10)]];

(* 5. Polynomial least squares, degree 5 over 500 points. *)
pts2 = N[Table[{i/100, (i/100)^3 - 2 (i/100) + 1}, {i, 1, 500}]];
bench["Fit degree-5 polynomial, 500 points",
  Fit[pts2, {1, t, t^2, t^3, t^4, t^5}, t];];
check["Fit degree-5 polynomial, 500 points",
  Round[10^4 (Fit[pts2, {1, t, t^2, t^3, t^4, t^5}, t] /. t -> 2)]];
