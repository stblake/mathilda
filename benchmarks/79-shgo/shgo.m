(* Experiment 79 -- SHGO (Simplicial Homology Global Optimization) global
   optimization (NMinimize, Method -> {"SHGO", ...}).

   MEASURES src/numerical_calculus/findmin.c's nm_shgo -- the method samples the
   bounded box (the exact Kuhn-triangulation 1-skeleton for "Simplicial", a
   k-nearest-neighbour graph over a Sobol'/Halton low-discrepancy set otherwise),
   identifies the minimizer pool (vertices better than every graph neighbour),
   and starts one local polish from each pool vertex -- against
   scipy.optimize.shgo. Same problems, same labels, same order, same
   sampling_method and sample count on both sides, so the race is
   SHGO-vs-shgo. These are the standard bounded multimodal benchmarks where the
   pool + homology construction shines: Himmelblau, Branin, six-hump camel,
   Cross-in-tray, Shubert, Rastrigin, Ackley, Styblinski-Tang. The check is the
   objective at the global optimum rounded to 10^6, which both systems reach and
   polish to agreement. Objectives over several variables are built with
   Sum/Table over fresh global symbols and passed through Evaluate[] past
   NMinimize's argument handling. *)

Get["../harness.m"];

require[{"NMinimize"}];

(* n-D objectives over fresh global symbols. *)
rastrigin[n_] := With[{v = Table[Symbol["rv" <> ToString[i]], {i, 1, n}]},
                      10 n + Sum[v[[i]]^2 - 10 Cos[2 Pi v[[i]]], {i, 1, n}]];
rastriginV[n_] := Table[Symbol["rv" <> ToString[i]], {i, 1, n}];
rastriginB[n_] := And @@ Table[-5.12 <= Symbol["rv" <> ToString[i]] <= 5.12, {i, 1, n}];

ackley[n_] := With[{v = Table[Symbol["av" <> ToString[i]], {i, 1, n}]},
                   -20 Exp[-0.2 Sqrt[Sum[v[[i]]^2, {i, 1, n}]/n]]
                   - Exp[Sum[Cos[2 Pi v[[i]]], {i, 1, n}]/n] + 20 + E];
ackleyV[n_] := Table[Symbol["av" <> ToString[i]], {i, 1, n}];
ackleyB[n_] := And @@ Table[-5 <= Symbol["av" <> ToString[i]] <= 5, {i, 1, n}];

styb[n_] := With[{v = Table[Symbol["sv" <> ToString[i]], {i, 1, n}]},
                 Sum[v[[i]]^4 - 16 v[[i]]^2 + 5 v[[i]], {i, 1, n}]/2];
stybV[n_] := Table[Symbol["sv" <> ToString[i]], {i, 1, n}];
stybB[n_] := And @@ Table[-5 <= Symbol["sv" <> ToString[i]] <= 5, {i, 1, n}];

(* Precompute out of the timed region. *)
ra3 = rastrigin[3]; ra3v = rastriginV[3]; ra3b = rastriginB[3];
ac4 = ackley[4];    ac4v = ackleyV[4];    ac4b = ackleyB[4];
st4 = styb[4];      st4v = stybV[4];      st4b = stybB[4];

(* 2-D benchmarks, simplicial sampling. *)
bench["01 Himmelblau simplicial", NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->100}];];
check["01 Himmelblau simplicial", Round[10^6 First[NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->100}]]]];

bench["02 Branin simplicial", NMinimize[{(y-5.1/(4 Pi^2) x^2+5/Pi x-6)^2+10(1-1/(8 Pi))Cos[x]+10, -5<=x<=10 && 0<=y<=15}, {x,y}, Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->100}];];
check["02 Branin simplicial", Round[10^6 First[NMinimize[{(y-5.1/(4 Pi^2) x^2+5/Pi x-6)^2+10(1-1/(8 Pi))Cos[x]+10, -5<=x<=10 && 0<=y<=15}, {x,y}, Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->100}]]]];

bench["03 Six-hump camel simplicial", NMinimize[{(4-2.1 x^2+x^4/3)x^2+x y+(-4+4 y^2)y^2, -3<=x<=3 && -2<=y<=2}, {x,y}, Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->200}];];
check["03 Six-hump camel simplicial", Round[10^6 First[NMinimize[{(4-2.1 x^2+x^4/3)x^2+x y+(-4+4 y^2)y^2, -3<=x<=3 && -2<=y<=2}, {x,y}, Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->200}]]]];

bench["04 Cross-in-tray simplicial", NMinimize[{-0.0001(Abs[Sin[x]Sin[y]Exp[Abs[100-Sqrt[x^2+y^2]/Pi]]]+1)^0.1, -10<=x<=10 && -10<=y<=10}, {x,y}, Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->500}];];
check["04 Cross-in-tray simplicial", Round[10^6 First[NMinimize[{-0.0001(Abs[Sin[x]Sin[y]Exp[Abs[100-Sqrt[x^2+y^2]/Pi]]]+1)^0.1, -10<=x<=10 && -10<=y<=10}, {x,y}, Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->500}]]]];

bench["05 Shubert simplicial", NMinimize[{(Sum[j Cos[(j+1)x+j],{j,1,5}])(Sum[j Cos[(j+1)y+j],{j,1,5}]), -10<=x<=10 && -10<=y<=10}, {x,y}, Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->600}];];
check["05 Shubert simplicial", Round[10^6 First[NMinimize[{(Sum[j Cos[(j+1)x+j],{j,1,5}])(Sum[j Cos[(j+1)y+j],{j,1,5}]), -10<=x<=10 && -10<=y<=10}, {x,y}, Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->600}]]]];

(* n-D benchmarks, Sobol'/Simplicial. *)
bench["06 Rastrigin 3D simplicial", NMinimize[{Evaluate[ra3], Evaluate[ra3b]}, Evaluate[ra3v], Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->200}];];
check["06 Rastrigin 3D simplicial", Round[10^6 First[NMinimize[{Evaluate[ra3], Evaluate[ra3b]}, Evaluate[ra3v], Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->200}]]]];

bench["07 Ackley 4D simplicial", NMinimize[{Evaluate[ac4], Evaluate[ac4b]}, Evaluate[ac4v], Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->200}];];
check["07 Ackley 4D simplicial", Round[10^6 First[NMinimize[{Evaluate[ac4], Evaluate[ac4b]}, Evaluate[ac4v], Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->200}]]]];

bench["08 Styblinski-Tang 4D sobol", NMinimize[{Evaluate[st4], Evaluate[st4b]}, Evaluate[st4v], Method->{"SHGO","SamplingMethod"->"Sobol","SearchPoints"->256}];];
check["08 Styblinski-Tang 4D sobol", Round[10^6 First[NMinimize[{Evaluate[st4], Evaluate[st4b]}, Evaluate[st4v], Method->{"SHGO","SamplingMethod"->"Sobol","SearchPoints"->256}]]]];
