(* Experiment 83 -- DIRECT (DIviding RECTangles) global optimization
   (NMinimize, Method -> {"DIRECT", ...}).

   MEASURES src/numerical_calculus/nm_direct.c -- a deterministic Lipschitzian
   search that normalizes the box to the unit hypercube and repeatedly subdivides
   the "potentially optimal" cells (Jones, Perttunen & Stuckman 1993; locally-
   biased DIRECT-L of Gablonsky & Kelley 2001) -- against scipy.optimize.direct.
   Same problems, same labels, same order, scipy's default parameters on both
   sides (locally_biased True, eps 1e-4, maxiter 1000, maxfun 1000 n). DIRECT is
   deterministic: there is no seed. The race is raw-vs-raw -- scipy's direct does
   no local polish, so "PostProcess" -> False disables Mathilda's for an apples-
   to-apples comparison of the same algorithm. The check is the objective at the
   reported point rounded to 10^3: two raw DIRECT runs agree at the basin, not to
   machine precision, so the rounding is coarser than the polished-engine
   benchmarks. Objectives over several variables are built with Sum/Table over
   fresh global symbols and passed through Evaluate[] past NMinimize's argument
   handling. *)

Get["../harness.m"];

require[{"NMinimize"}];

(* n-D Ackley over fresh global symbols, precomputed out of the timed region. *)
ackley[n_] := With[{v = Table[Symbol["av" <> ToString[i]], {i, 1, n}]},
                   -20 Exp[-0.2 Sqrt[Sum[v[[i]]^2, {i, 1, n}]/n]]
                   - Exp[Sum[Cos[2 Pi v[[i]]], {i, 1, n}]/n] + 20 + E];
ackleyV[n_] := Table[Symbol["av" <> ToString[i]], {i, 1, n}];
ackleyB[n_] := And @@ Table[-5 <= Symbol["av" <> ToString[i]] <= 5, {i, 1, n}];

ac5 = ackley[5]; ac5v = ackleyV[5]; ac5b = ackleyB[5];

(* 2-D benchmarks (bare x, y so the objective hits the compiled machine-precision
   path -- the fair analogue of scipy's numpy objective). *)
bench["01 Himmelblau 2D", NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DIRECT","PostProcess"->False}];];
check["01 Himmelblau 2D", Round[10^3 First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DIRECT","PostProcess"->False}]]]];

bench["02 Booth 2D", NMinimize[{(x+2y-7)^2+(2x+y-5)^2, -10<=x<=10 && -10<=y<=10}, {x,y}, Method->{"DIRECT","PostProcess"->False}];];
check["02 Booth 2D", Round[10^3 First[NMinimize[{(x+2y-7)^2+(2x+y-5)^2, -10<=x<=10 && -10<=y<=10}, {x,y}, Method->{"DIRECT","PostProcess"->False}]]]];

bench["03 Beale 2D", NMinimize[{(1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2, -4.5<=x<=4.5 && -4.5<=y<=4.5}, {x,y}, Method->{"DIRECT","PostProcess"->False}];];
check["03 Beale 2D", Round[10^3 First[NMinimize[{(1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2, -4.5<=x<=4.5 && -4.5<=y<=4.5}, {x,y}, Method->{"DIRECT","PostProcess"->False}]]]];

bench["04 Rosenbrock 2D", NMinimize[{100(y-x^2)^2+(1-x)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DIRECT","PostProcess"->False}];];
check["04 Rosenbrock 2D", Round[10^3 First[NMinimize[{100(y-x^2)^2+(1-x)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DIRECT","PostProcess"->False}]]]];

bench["05 Six-hump camel 2D", NMinimize[{(4-2.1 x^2+x^4/3)x^2+x y+(-4+4 y^2)y^2, -3<=x<=3 && -2<=y<=2}, {x,y}, Method->{"DIRECT","PostProcess"->False}];];
check["05 Six-hump camel 2D", Round[10^3 First[NMinimize[{(4-2.1 x^2+x^4/3)x^2+x y+(-4+4 y^2)y^2, -3<=x<=3 && -2<=y<=2}, {x,y}, Method->{"DIRECT","PostProcess"->False}]]]];

bench["06 Ackley 2D", NMinimize[{-20 Exp[-0.2 Sqrt[(x^2+y^2)/2]] - Exp[(Cos[2 Pi x]+Cos[2 Pi y])/2] + 20 + E, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DIRECT","PostProcess"->False}];];
check["06 Ackley 2D", Round[10^3 First[NMinimize[{-20 Exp[-0.2 Sqrt[(x^2+y^2)/2]] - Exp[(Cos[2 Pi x]+Cos[2 Pi y])/2] + 20 + E, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DIRECT","PostProcess"->False}]]]];

(* n-D benchmark (fresh global symbols, injected with Evaluate[]). *)
bench["07 Ackley 5D", NMinimize[{Evaluate[ac5], Evaluate[ac5b]}, Evaluate[ac5v], Method->{"DIRECT","PostProcess"->False}];];
check["07 Ackley 5D", Round[10^3 First[NMinimize[{Evaluate[ac5], Evaluate[ac5b]}, Evaluate[ac5v], Method->{"DIRECT","PostProcess"->False}]]]];

(* Bare x, y (the 2-D form written out) for the compiled machine-precision path. *)
bench["08 Styblinski-Tang 2D", NMinimize[{(x^4-16 x^2+5 x + y^4-16 y^2+5 y)/2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DIRECT","PostProcess"->False}];];
check["08 Styblinski-Tang 2D", Round[10^3 First[NMinimize[{(x^4-16 x^2+5 x + y^4-16 y^2+5 y)/2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DIRECT","PostProcess"->False}]]]];

bench["09 Rastrigin 2D", NMinimize[{20 + x^2-10 Cos[2 Pi x] + y^2-10 Cos[2 Pi y], -5.12<=x<=5.12 && -5.12<=y<=5.12}, {x,y}, Method->{"DIRECT","PostProcess"->False}];];
check["09 Rastrigin 2D", Round[10^3 First[NMinimize[{20 + x^2-10 Cos[2 Pi x] + y^2-10 Cos[2 Pi y], -5.12<=x<=5.12 && -5.12<=y<=5.12}, {x,y}, Method->{"DIRECT","PostProcess"->False}]]]];
