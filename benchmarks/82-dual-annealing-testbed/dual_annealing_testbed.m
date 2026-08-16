(* Experiment 82 -- Dual Annealing across the NMinimize test corpus.

   Runs the prior NMinimize stress / unit-test problems (from
   tests/test_nminimize.c) through the new Dual Annealing engine
   (Method -> {"DualAnnealing", ...}) and compares each to
   scipy.optimize.dual_annealing on the identical bounded box, with scipy's
   default parameters and a matched seed on both sides. This is the
   cleanly-comparable subset -- the continuous, box-bounded multimodal stress
   functions where both implementations reach the same global. The rest of the
   corpus (constrained, equality-constrained, disjunctive, combinatorial /
   mixed-integer, and the sharp-needle functions both miss) is analysed in
   README.md; scipy.optimize.dual_annealing is BOX-ONLY (no constraints
   argument at all), so those are not Dual-Annealing-vs-dual_annealing races.

   The check is the objective at the global optimum, rounded per problem so that
   both implementations, reaching the same basin, agree. Multi-variable
   objectives are built over fresh global symbols and passed through Evaluate[]
   past NMinimize's argument handling (and so they compile to the machine path).

   NB: the precomputed-objective accumulators (pAck3, pRos10, ...) are named so
   they can never collide with a problem variable Symbol[prefix <> ToString[i]] --
   assigning e.g. `ro10 = rosenbrock[10]` would make the symbol ro10 (the 10th
   variable) self-referential and hang the file. *)

Get["../harness.m"];

require[{"NMinimize"}];

M[sd_] := {"DualAnnealing", "RandomSeed" -> sd};

(* Fresh-symbol objective / variable / box triples, precomputed out of the
   timed region. Separate obj/var/box functions (the 79-shgo idiom). *)
ackleyObj[n_] := With[{v = Table[Symbol["ak" <> ToString[i]], {i, 1, n}]},
    -20 Exp[-0.2 Sqrt[Total[Table[v[[i]]^2, {i, 1, n}]]/n]]
    - Exp[Total[Table[Cos[2 Pi v[[i]]], {i, 1, n}]]/n] + 20 + E];
ackleyVar[n_] := Table[Symbol["ak" <> ToString[i]], {i, 1, n}];
ackleyBox[n_, lo_, hi_] := And @@ Table[lo <= Symbol["ak" <> ToString[i]] <= hi, {i, 1, n}];

rastObj[n_] := With[{v = Table[Symbol["rn" <> ToString[i]], {i, 1, n}]},
    10 n + Total[Table[v[[i]]^2 - 10 Cos[2 Pi v[[i]]], {i, 1, n}]]];
rastVar[n_] := Table[Symbol["rn" <> ToString[i]], {i, 1, n}];
rastBox[n_] := And @@ Table[-5.12 <= Symbol["rn" <> ToString[i]] <= 5.12, {i, 1, n}];

rosObj[n_] := With[{v = Table[Symbol["ro" <> ToString[i]], {i, 1, n}]},
    Total[Table[100 (v[[i + 1]] - v[[i]]^2)^2 + (1 - v[[i]])^2, {i, 1, n - 1}]]];
rosVar[n_] := Table[Symbol["ro" <> ToString[i]], {i, 1, n}];
rosBox[n_] := And @@ Table[-5 <= Symbol["ro" <> ToString[i]] <= 5, {i, 1, n}];

schwObj[n_] := With[{v = Table[Symbol["sw" <> ToString[i]], {i, 1, n}]},
    418.9829 n - Total[Table[v[[i]] Sin[Sqrt[Abs[v[[i]]]]], {i, 1, n}]]];
schwVar[n_] := Table[Symbol["sw" <> ToString[i]], {i, 1, n}];
schwBox[n_] := And @@ Table[-500 <= Symbol["sw" <> ToString[i]] <= 500, {i, 1, n}];

pAck3 = ackleyObj[3];  pAck3v = ackleyVar[3];  pAck3b = ackleyBox[3, -32, 32];
pRas5 = rastObj[5];    pRas5v = rastVar[5];    pRas5b = rastBox[5];
pRas8 = rastObj[8];    pRas8v = rastVar[8];    pRas8b = rastBox[8];
pRos10 = rosObj[10];   pRos10v = rosVar[10];   pRos10b = rosBox[10];
pSch10 = schwObj[10];  pSch10v = schwVar[10];  pSch10b = schwBox[10];

(* --- 1-D / 2-D, bare variables (compiled machine path) --- *)
bench["01 Chained trig 1D", NMinimize[{Sin[2 x] + Cos[x], -2 <= x <= 3}, {x}, Method -> M[1]];];
check["01 Chained trig 1D", Round[10^4 First[NMinimize[{Sin[2 x] + Cos[x], -2 <= x <= 3}, {x}, Method -> M[1]]]]];

bench["02 Schwefel 2D", NMinimize[{837.9658 - x Sin[Sqrt[Abs[x]]] - y Sin[Sqrt[Abs[y]]], -500 <= x <= 500 && -500 <= y <= 500}, {x, y}, Method -> M[1]];];
check["02 Schwefel 2D", Round[10^3 First[NMinimize[{837.9658 - x Sin[Sqrt[Abs[x]]] - y Sin[Sqrt[Abs[y]]], -500 <= x <= 500 && -500 <= y <= 500}, {x, y}, Method -> M[1]]]]];

bench["03 Schaffer N2 2D", NMinimize[{0.5 + (Sin[x^2 - y^2]^2 - 0.5)/(1 + 0.001 (x^2 + y^2))^2, -100 <= x <= 100 && -100 <= y <= 100}, {x, y}, Method -> M[1]];];
check["03 Schaffer N2 2D", Round[10^4 First[NMinimize[{0.5 + (Sin[x^2 - y^2]^2 - 0.5)/(1 + 0.001 (x^2 + y^2))^2, -100 <= x <= 100 && -100 <= y <= 100}, {x, y}, Method -> M[1]]]]];

bench["04 Rugged sine 2D", NMinimize[{x^2 + y^2 + 10 Sin[3 x]^2 + 10 Sin[3 y]^2, -5 <= x <= 5 && -5 <= y <= 5}, {x, y}, Method -> M[1]];];
check["04 Rugged sine 2D", Round[10^4 First[NMinimize[{x^2 + y^2 + 10 Sin[3 x]^2 + 10 Sin[3 y]^2, -5 <= x <= 5 && -5 <= y <= 5}, {x, y}, Method -> M[1]]]]];

(* --- 3-D..10-D, fresh-symbol objectives --- *)
bench["05 Ackley 3D", NMinimize[{Evaluate[pAck3], Evaluate[pAck3b]}, Evaluate[pAck3v], Method -> M[1]];];
check["05 Ackley 3D", Round[10^4 First[NMinimize[{Evaluate[pAck3], Evaluate[pAck3b]}, Evaluate[pAck3v], Method -> M[1]]]]];

bench["06 Rastrigin 5D", NMinimize[{Evaluate[pRas5], Evaluate[pRas5b]}, Evaluate[pRas5v], Method -> M[1]];];
check["06 Rastrigin 5D", Round[10^4 First[NMinimize[{Evaluate[pRas5], Evaluate[pRas5b]}, Evaluate[pRas5v], Method -> M[1]]]]];

bench["07 Rastrigin 8D", NMinimize[{Evaluate[pRas8], Evaluate[pRas8b]}, Evaluate[pRas8v], Method -> M[1]];];
check["07 Rastrigin 8D", Round[10^4 First[NMinimize[{Evaluate[pRas8], Evaluate[pRas8b]}, Evaluate[pRas8v], Method -> M[1]]]]];

bench["08 Rosenbrock 10D", NMinimize[{Evaluate[pRos10], Evaluate[pRos10b]}, Evaluate[pRos10v], Method -> M[1]];];
check["08 Rosenbrock 10D", Round[10^4 First[NMinimize[{Evaluate[pRos10], Evaluate[pRos10b]}, Evaluate[pRos10v], Method -> M[1]]]]];

bench["09 Schwefel 10D", NMinimize[{Evaluate[pSch10], Evaluate[pSch10b]}, Evaluate[pSch10v], Method -> M[1]];];
check["09 Schwefel 10D", Round[10^4 First[NMinimize[{Evaluate[pSch10], Evaluate[pSch10b]}, Evaluate[pSch10v], Method -> M[1]]]]];
