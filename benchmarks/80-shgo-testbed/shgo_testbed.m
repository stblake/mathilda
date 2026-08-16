(* Experiment 80 -- SHGO across the NMinimize test corpus (NMinimize,
   Method -> {"SHGO", ...}).

   This runs the *prior NMinimize stress/unit test problems* (from
   tests/test_nminimize.c) through the new SHGO engine and races each against
   scipy.optimize.shgo on the identical problem (same box, sampling method, and
   sample count). It is the cleanly-comparable set -- the multimodal stress
   functions and the inequality-constrained problems within SHGO's applicability
   envelope, where both implementations reach the same global and the timing
   comparison is meaningful. The rest of the corpus is analysed in README.md:
   the constrained problems that scipy.optimize.shgo mishandles (Mathilda solves
   them via its augmented-Lagrangian polish), the high-dimensional / sharp-needle
   functions outside SHGO's envelope (both implementations miss the global), and
   the combinatorial-integer problems (scipy.optimize.shgo is continuous-only).

   The check is the objective at the global optimum, rounded per problem so that
   both implementations, reaching the same basin, agree. Indexed objectives use
   x[i] exactly as the original tests do. *)

Get["../harness.m"];

require[{"NMinimize"}];

(* ================= Multimodal battery (low-dimensional) ================= *)

bench["01 Eggholder 2D", NMinimize[{-(y+47)Sin[Sqrt[Abs[y+x/2+47]]]-x Sin[Sqrt[Abs[x-(y+47)]]], -512<=x<=512 && -512<=y<=512}, {x,y}, Method->{"SHGO","SearchPoints"->500}];];
check["01 Eggholder 2D", Round[First[NMinimize[{-(y+47)Sin[Sqrt[Abs[y+x/2+47]]]-x Sin[Sqrt[Abs[x-(y+47)]]], -512<=x<=512 && -512<=y<=512}, {x,y}, Method->{"SHGO","SearchPoints"->500}]]]];

bench["02 Schwefel 2D", NMinimize[{837.9658-x Sin[Sqrt[Abs[x]]]-y Sin[Sqrt[Abs[y]]], -500<=x<=500 && -500<=y<=500}, {x,y}, Method->{"SHGO","SearchPoints"->300}];];
check["02 Schwefel 2D", Round[10^3 First[NMinimize[{837.9658-x Sin[Sqrt[Abs[x]]]-y Sin[Sqrt[Abs[y]]], -500<=x<=500 && -500<=y<=500}, {x,y}, Method->{"SHGO","SearchPoints"->300}]]]];

bench["03 Schaffer N2 2D", NMinimize[{0.5+(Sin[x^2-y^2]^2-0.5)/(1+0.001(x^2+y^2))^2, -100<=x<=100 && -100<=y<=100}, {x,y}, Method->{"SHGO","SearchPoints"->400}];];
check["03 Schaffer N2 2D", Round[10^4 First[NMinimize[{0.5+(Sin[x^2-y^2]^2-0.5)/(1+0.001(x^2+y^2))^2, -100<=x<=100 && -100<=y<=100}, {x,y}, Method->{"SHGO","SearchPoints"->400}]]]];

bench["04 Mishra Bird disk", NMinimize[{Sin[y]Exp[(1-Cos[x])^2]+Cos[x]Exp[(1-Sin[y])^2]+(x-y)^2, (x+5)^2+(y+5)^2<25}, {x,y}, Method->{"SHGO","SearchPoints"->300}];];
check["04 Mishra Bird disk", Round[10^3 First[NMinimize[{Sin[y]Exp[(1-Cos[x])^2]+Cos[x]Exp[(1-Sin[y])^2]+(x-y)^2, (x+5)^2+(y+5)^2<25}, {x,y}, Method->{"SHGO","SearchPoints"->300}]]]];

bench["05 Rugged sine 2D", NMinimize[{x^2+y^2+10 Sin[3x]^2+10 Sin[3y]^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"SHGO","SearchPoints"->200}];];
check["05 Rugged sine 2D", Round[10^4 First[NMinimize[{x^2+y^2+10 Sin[3x]^2+10 Sin[3y]^2, -5<=x<=5 && -5<=y<=5}, {x,y}, Method->{"SHGO","SearchPoints"->200}]]]];

bench["06 Ackley 3D", NMinimize[{-20 Exp[-0.2 Sqrt[1/3 Sum[x[i]^2,{i,1,3}]]]-Exp[1/3 Sum[Cos[2 Pi x[i]],{i,1,3}]]+20+E, Table[-32<=x[i]<=32,{i,1,3}]}, Table[x[i],{i,1,3}], Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->200}];];
check["06 Ackley 3D", Round[10^4 First[NMinimize[{-20 Exp[-0.2 Sqrt[1/3 Sum[x[i]^2,{i,1,3}]]]-Exp[1/3 Sum[Cos[2 Pi x[i]],{i,1,3}]]+20+E, Table[-32<=x[i]<=32,{i,1,3}]}, Table[x[i],{i,1,3}], Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->200}]]]];

bench["07 Rastrigin 5D", NMinimize[{50+Sum[x[i]^2-10 Cos[2 Pi x[i]],{i,1,5}], Table[-5.12<=x[i]<=5.12,{i,1,5}]}, Table[x[i],{i,1,5}], Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->300}];];
check["07 Rastrigin 5D", Round[10^4 First[NMinimize[{50+Sum[x[i]^2-10 Cos[2 Pi x[i]],{i,1,5}], Table[-5.12<=x[i]<=5.12,{i,1,5}]}, Table[x[i],{i,1,5}], Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->300}]]]];

bench["08 Griewank 5D", NMinimize[{Sum[x[i]^2/4000,{i,1,5}]-Product[Cos[x[i]/Sqrt[i]],{i,1,5}]+1, Table[-15<=x[i]<=15,{i,1,5}]}, Table[x[i],{i,1,5}], Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->300}];];
check["08 Griewank 5D", Round[10^4 First[NMinimize[{Sum[x[i]^2/4000,{i,1,5}]-Product[Cos[x[i]/Sqrt[i]],{i,1,5}]+1, Table[-15<=x[i]<=15,{i,1,5}]}, Table[x[i],{i,1,5}], Method->{"SHGO","SamplingMethod"->"Simplicial","SearchPoints"->300}]]]];

(* ================= Inequality-constrained ================= *)

bench["09 Disk linear", NMinimize[{x+y, x^2+y^2<=9}, {x,y}, Method->{"SHGO","SearchPoints"->100}];];
check["09 Disk linear", Round[10^4 First[NMinimize[{x+y, x^2+y^2<=9}, {x,y}, Method->{"SHGO","SearchPoints"->100}]]]];

bench["10 Chained trig 1D", NMinimize[{Sin[2x]+Cos[x], -2<=x<=3}, {x}, Method->{"SHGO","SearchPoints"->100}];];
check["10 Chained trig 1D", Round[10^4 First[NMinimize[{Sin[2x]+Cos[x], -2<=x<=3}, {x}, Method->{"SHGO","SearchPoints"->100}]]]];

bench["11 Quadratic halfplane", NMinimize[{x^2+y^2, x+y>=1}, {x,y}, Method->{"SHGO","SearchPoints"->100}];];
check["11 Quadratic halfplane", Round[10^4 First[NMinimize[{x^2+y^2, x+y>=1}, {x,y}, Method->{"SHGO","SearchPoints"->100}]]]];
