(* Experiment 81 -- Dual Annealing (Generalized Simulated Annealing) global
   optimization (NMinimize, Method -> {"DualAnnealing", ...}).

   MEASURES src/numerical_calculus/findmin.c's nm_dual_annealing -- a heavy-tailed
   Tsallis visiting distribution + generalized Metropolis acceptance with
   reannealing, plus a local search after each Markov chain -- against
   scipy.optimize.dual_annealing. Same problems, same labels, same order, scipy's
   default parameters on both sides (visit 2.62, accept -5, initial_temp 5230,
   maxiter 1000), with a fixed RandomSeed/seed so both runs are reproducible.
   These are standard bounded multimodal benchmarks whose global both systems
   reach and polish to agreement; the check is the objective at the global optimum
   rounded to 10^6. Objectives over several variables are built with Sum/Table
   over fresh global symbols and passed through Evaluate[] past NMinimize's
   argument handling. *)

Get["../harness.m"];

require[{"NMinimize"}];

(* n-D Ackley over fresh global symbols, precomputed out of the timed region. *)
ackley[n_] := With[{v = Table[Symbol["av" <> ToString[i]], {i, 1, n}]},
                   -20 Exp[-0.2 Sqrt[Sum[v[[i]]^2, {i, 1, n}]/n]]
                   - Exp[Sum[Cos[2 Pi v[[i]]], {i, 1, n}]/n] + 20 + E];
ackleyV[n_] := Table[Symbol["av" <> ToString[i]], {i, 1, n}];
ackleyB[n_] := And @@ Table[-5 <= Symbol["av" <> ToString[i]] <= 5, {i, 1, n}];

ac5 = ackley[5]; ac5v = ackleyV[5]; ac5b = ackleyB[5];

(* 2-D benchmarks. *)
bench["01 Himmelblau 2D", NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}];];
check["01 Himmelblau 2D", Round[10^6 First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}]]]];

bench["02 Booth 2D", NMinimize[{(x+2y-7)^2+(2x+y-5)^2, -10<=x<=10 && -10<=y<=10}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}];];
check["02 Booth 2D", Round[10^6 First[NMinimize[{(x+2y-7)^2+(2x+y-5)^2, -10<=x<=10 && -10<=y<=10}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}]]]];

bench["03 Beale 2D", NMinimize[{(1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2, -4.5<=x<=4.5 && -4.5<=y<=4.5}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}];];
check["03 Beale 2D", Round[10^6 First[NMinimize[{(1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2, -4.5<=x<=4.5 && -4.5<=y<=4.5}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}]]]];

bench["04 Rosenbrock 2D", NMinimize[{100(y-x^2)^2+(1-x)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}];];
check["04 Rosenbrock 2D", Round[10^6 First[NMinimize[{100(y-x^2)^2+(1-x)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}]]]];

bench["05 Six-hump camel 2D", NMinimize[{(4-2.1 x^2+x^4/3)x^2+x y+(-4+4 y^2)y^2, -3<=x<=3 && -2<=y<=2}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}];];
check["05 Six-hump camel 2D", Round[10^6 First[NMinimize[{(4-2.1 x^2+x^4/3)x^2+x y+(-4+4 y^2)y^2, -3<=x<=3 && -2<=y<=2}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}]]]];

bench["06 Ackley 2D", NMinimize[{-20 Exp[-0.2 Sqrt[(x^2+y^2)/2]] - Exp[(Cos[2 Pi x]+Cos[2 Pi y])/2] + 20 + E, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}];];
check["06 Ackley 2D", Round[10^6 First[NMinimize[{-20 Exp[-0.2 Sqrt[(x^2+y^2)/2]] - Exp[(Cos[2 Pi x]+Cos[2 Pi y])/2] + 20 + E, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}]]]];

(* n-D benchmark. *)
bench["07 Ackley 5D", NMinimize[{Evaluate[ac5], Evaluate[ac5b]}, Evaluate[ac5v], Method->{"DualAnnealing","RandomSeed"->1}];];
check["07 Ackley 5D", Round[10^6 First[NMinimize[{Evaluate[ac5], Evaluate[ac5b]}, Evaluate[ac5v], Method->{"DualAnnealing","RandomSeed"->1}]]]];

(* Bare x, y (the 2-D forms written out) so the objective hits the compiled
   machine-precision path -- the fair analogue of scipy's numpy objective. An
   inline Sum[...x[i]...] would run on the interpreter and is ~1000x slower here,
   a Mathilda indexed-variable limitation unrelated to the Dual Annealing engine. *)
bench["08 Styblinski-Tang 2D", NMinimize[{(x^4-16 x^2+5 x + y^4-16 y^2+5 y)/2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}];];
check["08 Styblinski-Tang 2D", Round[10^6 First[NMinimize[{(x^4-16 x^2+5 x + y^4-16 y^2+5 y)/2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}]]]];

bench["09 Rastrigin 2D", NMinimize[{20 + x^2-10 Cos[2 Pi x] + y^2-10 Cos[2 Pi y], -5.12<=x<=5.12 && -5.12<=y<=5.12}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}];];
check["09 Rastrigin 2D", Round[10^6 First[NMinimize[{20 + x^2-10 Cos[2 Pi x] + y^2-10 Cos[2 Pi y], -5.12<=x<=5.12 && -5.12<=y<=5.12}, {x,y}, Method->{"DualAnnealing","RandomSeed"->1}]]]];
