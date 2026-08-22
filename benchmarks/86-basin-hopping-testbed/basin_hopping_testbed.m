(* Experiment 86 -- Basin Hopping across the NMinimize test corpus.

   Runs the prior NMinimize stress / unit-test problems (from
   tests/test_nminimize.c) through the Basin Hopping engine
   (Method -> {"BasinHopping", ...}) and compares each to
   scipy.optimize.basinhopping on the identical bounded box, with a MATCHED
   parameterisation on both sides: K = 8 independent multi-start runs (Mathilda's
   "SearchPoints" -> 8, and 8 seeded scipy basinhopping runs kept at the min), and
   a box-scaled "StepSize" for the wide-box deceptive functions (the uniform 0.5
   displacement cannot traverse a [-500,500] box in 100 hops). "Appropriate
   parameterisation" = the same multi-start budget and step size on both sides, so
   the race stays like-for-like on the algorithm.

   This file is the cleanly-comparable subset -- continuous, box-bounded problems
   where both implementations, with the matched parameterisation, reach the same
   global. The rest of the corpus (the funnel cases where Mathilda's single-family
   walk reaches the global and scipy's stalls even at K=8; the deceptive high-D
   Schwefel both miss; the constrained / equality / disjunctive cases scipy's
   box-only basinhopping cannot express; and the combinatorial / mixed-integer
   cases) is analysed in README.md.

   The check is the objective at the global optimum, rounded per problem so both
   implementations, reaching the same basin, agree. Multi-variable objectives are
   built over fresh global symbols and passed through Evaluate[] past NMinimize's
   argument handling (so they compile to the machine path). The precomputed-
   objective accumulators (pAck3, ...) are named so they can never collide with a
   problem variable Symbol[prefix <> ToString[i]]. *)

Get["../harness.m"];

require[{"NMinimize"}];

(* K = 8 multi-start; box-scaled StepSize for the wide-box deceptive functions. *)
M[sd_]           := {"BasinHopping", "RandomSeed" -> sd, "SearchPoints" -> 8};
Ms[sd_, step_]   := {"BasinHopping", "RandomSeed" -> sd, "SearchPoints" -> 8, "StepSize" -> step};

(* Fresh-symbol objective / variable / box triples (the 79-shgo / 82 idiom). *)
ackleyObj[n_] := With[{v = Table[Symbol["ak" <> ToString[i]], {i, 1, n}]},
    -20 Exp[-0.2 Sqrt[Total[Table[v[[i]]^2, {i, 1, n}]]/n]]
    - Exp[Total[Table[Cos[2 Pi v[[i]]], {i, 1, n}]]/n] + 20 + E];
ackleyVar[n_] := Table[Symbol["ak" <> ToString[i]], {i, 1, n}];
ackleyBox[n_, lo_, hi_] := And @@ Table[lo <= Symbol["ak" <> ToString[i]] <= hi, {i, 1, n}];

rosObj[n_] := With[{v = Table[Symbol["ro" <> ToString[i]], {i, 1, n}]},
    Total[Table[100 (v[[i + 1]] - v[[i]]^2)^2 + (1 - v[[i]])^2, {i, 1, n - 1}]]];
rosVar[n_] := Table[Symbol["ro" <> ToString[i]], {i, 1, n}];
rosBox[n_] := And @@ Table[-5 <= Symbol["ro" <> ToString[i]] <= 5, {i, 1, n}];

pAck3 = ackleyObj[3];  pAck3v = ackleyVar[3];  pAck3b = ackleyBox[3, -32, 32];
pRos10 = rosObj[10];   pRos10v = rosVar[10];   pRos10b = rosBox[10];

(* --- 1-D / 2-D, bare variables (compiled machine path) --- *)
bench["01 Chained trig 1D", NMinimize[{Sin[2 x] + Cos[x], -2 <= x <= 3}, {x}, Method -> M[1]];];
check["01 Chained trig 1D", Round[10^4 First[NMinimize[{Sin[2 x] + Cos[x], -2 <= x <= 3}, {x}, Method -> M[1]]]]];

bench["02 Rugged sine 2D", NMinimize[{x^2 + y^2 + 10 Sin[3 x]^2 + 10 Sin[3 y]^2, -5 <= x <= 5 && -5 <= y <= 5}, {x, y}, Method -> M[1]];];
check["02 Rugged sine 2D", Round[10^4 First[NMinimize[{x^2 + y^2 + 10 Sin[3 x]^2 + 10 Sin[3 y]^2, -5 <= x <= 5 && -5 <= y <= 5}, {x, y}, Method -> M[1]]]]];

bench["03 Schwefel 2D", NMinimize[{837.9658 - x Sin[Sqrt[Abs[x]]] - y Sin[Sqrt[Abs[y]]], -500 <= x <= 500 && -500 <= y <= 500}, {x, y}, Method -> Ms[1, 150]];];
check["03 Schwefel 2D", Round[10^3 First[NMinimize[{837.9658 - x Sin[Sqrt[Abs[x]]] - y Sin[Sqrt[Abs[y]]], -500 <= x <= 500 && -500 <= y <= 500}, {x, y}, Method -> Ms[1, 150]]]]];

bench["04 Schaffer N2 2D", NMinimize[{0.5 + (Sin[x^2 - y^2]^2 - 0.5)/(1 + 0.001 (x^2 + y^2))^2, -100 <= x <= 100 && -100 <= y <= 100}, {x, y}, Method -> Ms[1, 40]];];
check["04 Schaffer N2 2D", Round[10^4 First[NMinimize[{0.5 + (Sin[x^2 - y^2]^2 - 0.5)/(1 + 0.001 (x^2 + y^2))^2, -100 <= x <= 100 && -100 <= y <= 100}, {x, y}, Method -> Ms[1, 40]]]]];

(* --- 3-D / 10-D, fresh-symbol objectives --- *)
bench["05 Ackley 3D", NMinimize[{Evaluate[pAck3], Evaluate[pAck3b]}, Evaluate[pAck3v], Method -> M[1]];];
check["05 Ackley 3D", Round[10^4 First[NMinimize[{Evaluate[pAck3], Evaluate[pAck3b]}, Evaluate[pAck3v], Method -> M[1]]]]];

bench["06 Rosenbrock 10D", NMinimize[{Evaluate[pRos10], Evaluate[pRos10b]}, Evaluate[pRos10v], Method -> M[1]];];
check["06 Rosenbrock 10D", Round[10^4 First[NMinimize[{Evaluate[pRos10], Evaluate[pRos10b]}, Evaluate[pRos10v], Method -> M[1]]]]];
