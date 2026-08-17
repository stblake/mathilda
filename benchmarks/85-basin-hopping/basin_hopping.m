(* Experiment 85 -- Basin Hopping (Monte-Carlo minimization) global optimization
   (NMinimize, Method -> {"BasinHopping", ...}).

   MEASURES src/numerical_calculus/nm_basin_hopping.c -- each hop is a uniform
   random displacement + a local minimization ("quench") + a Metropolis accept on
   the two locally-minimized energies, with an adaptive step size targeting a fixed
   acceptance rate -- against scipy.optimize.basinhopping. Same problems, same
   labels, same order, scipy's default parameters on both sides (T 1, stepsize 0.5,
   interval 50, target_accept_rate 0.5, stepwise_factor 0.9, niter 100), with a
   fixed RandomSeed/seed so both runs are reproducible. These are eight standard
   bounded multimodal benchmarks whose global both a single Mathilda run and a
   single scipy run reach and polish to agreement; the check is the objective at the
   global optimum rounded to 10^6. Objectives over several variables are built with
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

(* 2-D benchmarks (bare x, y -> compiled machine-precision objective, the fair
   analogue of scipy's numpy objective). *)
bench["01 Himmelblau 2D", NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}];];
check["01 Himmelblau 2D", Round[10^6 First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}]]]];

bench["02 Booth 2D", NMinimize[{(x+2y-7)^2+(2x+y-5)^2, -10<=x<=10 && -10<=y<=10}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}];];
check["02 Booth 2D", Round[10^6 First[NMinimize[{(x+2y-7)^2+(2x+y-5)^2, -10<=x<=10 && -10<=y<=10}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}]]]];

bench["03 Beale 2D", NMinimize[{(1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2, -4.5<=x<=4.5 && -4.5<=y<=4.5}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}];];
check["03 Beale 2D", Round[10^6 First[NMinimize[{(1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2, -4.5<=x<=4.5 && -4.5<=y<=4.5}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}]]]];

bench["04 Rosenbrock 2D", NMinimize[{100(y-x^2)^2+(1-x)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}];];
check["04 Rosenbrock 2D", Round[10^6 First[NMinimize[{100(y-x^2)^2+(1-x)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}]]]];

bench["05 Six-hump camel 2D", NMinimize[{(4-2.1 x^2+x^4/3)x^2+x y+(-4+4 y^2)y^2, -3<=x<=3 && -2<=y<=2}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}];];
check["05 Six-hump camel 2D", Round[10^6 First[NMinimize[{(4-2.1 x^2+x^4/3)x^2+x y+(-4+4 y^2)y^2, -3<=x<=3 && -2<=y<=2}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}]]]];

bench["06 Ackley 2D", NMinimize[{-20 Exp[-0.2 Sqrt[(x^2+y^2)/2]] - Exp[(Cos[2 Pi x]+Cos[2 Pi y])/2] + 20 + E, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}];];
check["06 Ackley 2D", Round[10^6 First[NMinimize[{-20 Exp[-0.2 Sqrt[(x^2+y^2)/2]] - Exp[(Cos[2 Pi x]+Cos[2 Pi y])/2] + 20 + E, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}]]]];

bench["07 Sphere 2D", NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}];];
check["07 Sphere 2D", Round[10^6 First[NMinimize[{x^2+y^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"BasinHopping","RandomSeed"->1}]]]];

(* n-D benchmark: fresh symbols + Evaluate[] so the objective compiles. *)
bench["08 Ackley 5D", NMinimize[{Evaluate[ac5], Evaluate[ac5b]}, Evaluate[ac5v], Method->{"BasinHopping","RandomSeed"->1}];];
check["08 Ackley 5D", Round[10^6 First[NMinimize[{Evaluate[ac5], Evaluate[ac5b]}, Evaluate[ac5v], Method->{"BasinHopping","RandomSeed"->1}]]]];
