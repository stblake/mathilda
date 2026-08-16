(* Experiment 84 -- DIRECT across the NMinimize test corpus.

   Runs the RACEABLE (continuous, box-bounded) NMinimize stress / unit-test
   problems from tests/test_nminimize.c through the DIRECT engine
   (Method -> {"DIRECT", ...}) and compares each to scipy.optimize.direct on the
   identical box. DIRECT is deterministic -- no seed. Both sides use scipy's
   default parameters except where a case is annotated LF/N in its label
   ("LocallyBiased" -> False, "MaxFunctionEvaluations" -> N) -- the same
   appropriate parameterisation on both sides. The race is raw-vs-raw --
   scipy's direct does no local polish, so "PostProcess" -> False disables
   Mathilda's.

   Every case reaches the SAME point in both systems: the global for most, a
   shared non-global basin for Rosenbrock-10D (~8.7916) and Gaussian-well-10D
   (both raw DIRECTs stall at 1.0 -- the well at 1.2345 is too narrow for
   center-based subdivision). Schwefel-10D is excluded (the two deterministic
   subdivisions diverge into different missed basins).

   The rest of the corpus is analysed in README.md and is NOT a race:
   scipy.optimize.direct is BOX-ONLY, so the constrained / equality /
   disjunctive / integer problems have no scipy column (Mathilda's DIRECT solves
   the low-dimensional ones via its penalty + augmented-Lagrangian polish, which
   scipy's direct cannot do); the high-dimensional equality-constrained A-series
   and the combinatorial B-series are beyond DIRECT in either system.

   The check is the objective at the reported point, rounded per problem so both
   implementations agree at the basin. Multi-variable objectives are built over
   fresh global symbols and passed through Evaluate[] past NMinimize's argument
   handling (so they compile to the machine path). The accumulator globals are
   named so they can never collide with a problem variable Symbol[prefix <> i]. *)

Get["../harness.m"];

require[{"NMinimize"}];

R  = {"DIRECT", "PostProcess" -> False};                                  (* default params *)
Regg = {"DIRECT", "PostProcess" -> False, "LocallyBiased" -> False, "MaxFunctionEvaluations" -> 20000};
Rgw  = {"DIRECT", "PostProcess" -> False, "LocallyBiased" -> False, "MaxFunctionEvaluations" -> 40000};

(* Fresh-symbol objective / variable / box triples, precomputed out of the
   timed region (the 79-shgo idiom). *)
griObj[n_] := With[{v = Table[Symbol["gr" <> ToString[i]], {i, 1, n}]},
    Total[Table[v[[i]]^2/4000, {i, 1, n}]] - Product[Cos[v[[i]]/Sqrt[i]], {i, 1, n}] + 1];
griVar[n_] := Table[Symbol["gr" <> ToString[i]], {i, 1, n}];
griBox[n_, w_] := And @@ Table[-w <= Symbol["gr" <> ToString[i]] <= w, {i, 1, n}];

sqObj[n_] := With[{v = Table[Symbol["sq" <> ToString[i]], {i, 1, n}]},
    Total[Table[10^(6 (i - 1)/9) v[[i]]^2, {i, 1, n}]]];
sqVar[n_] := Table[Symbol["sq" <> ToString[i]], {i, 1, n}];
sqBox[n_] := And @@ Table[-10 <= Symbol["sq" <> ToString[i]] <= 10, {i, 1, n}];

qbObj[n_] := With[{v = Table[Symbol["qb" <> ToString[i]], {i, 1, n}]},
    Total[Table[v[[i]]^2, {i, 1, n}]]];
qbVar[n_] := Table[Symbol["qb" <> ToString[i]], {i, 1, n}];
qbBox[n_] := And @@ Table[-5 <= Symbol["qb" <> ToString[i]] <= 5, {i, 1, n}];

ackObj[n_] := With[{v = Table[Symbol["ak" <> ToString[i]], {i, 1, n}]},
    -20 Exp[-0.2 Sqrt[Total[Table[v[[i]]^2, {i, 1, n}]]/n]]
    - Exp[Total[Table[Cos[2 Pi v[[i]]], {i, 1, n}]]/n] + 20 + E];
ackVar[n_] := Table[Symbol["ak" <> ToString[i]], {i, 1, n}];
ackBox[n_, w_] := And @@ Table[-w <= Symbol["ak" <> ToString[i]] <= w, {i, 1, n}];

katObj[n_] := With[{v = Table[Symbol["kt" <> ToString[i]], {i, 1, n}]},
    (10/n^2) Product[1 + i Sum[Abs[2^k v[[i]] - Round[2^k v[[i]]]]/2^k, {k, 1, 25}], {i, 1, n}]^(10/n^1.2) - (10/n^2)];
katVar[n_] := Table[Symbol["kt" <> ToString[i]], {i, 1, n}];
katBox[n_] := And @@ Table[-100 <= Symbol["kt" <> ToString[i]] <= 100, {i, 1, n}];

rasObj[n_] := With[{v = Table[Symbol["rs" <> ToString[i]], {i, 1, n}]},
    10 n + Total[Table[v[[i]]^2 - 10 Cos[2 Pi v[[i]]], {i, 1, n}]]];
rasVar[n_] := Table[Symbol["rs" <> ToString[i]], {i, 1, n}];
rasBox[n_] := And @@ Table[-5.12 <= Symbol["rs" <> ToString[i]] <= 5.12, {i, 1, n}];

rosObj[n_] := With[{v = Table[Symbol["ro" <> ToString[i]], {i, 1, n}]},
    Total[Table[100 (v[[i + 1]] - v[[i]]^2)^2 + (1 - v[[i]])^2, {i, 1, n - 1}]]];
rosVar[n_] := Table[Symbol["ro" <> ToString[i]], {i, 1, n}];
rosBox[n_] := And @@ Table[-5 <= Symbol["ro" <> ToString[i]] <= 5, {i, 1, n}];

maObj[n_] := With[{v = Table[Symbol["ma" <> ToString[i]], {i, 1, n}]},
    -Exp[-0.2 Sqrt[Total[Table[v[[i]]^2, {i, 1, n}]]/n]] Product[Cos[20 v[[i]]], {i, 1, n}]
    + 0.05 Total[Table[v[[i]]^2, {i, 1, n}]]];
maVar[n_] := Table[Symbol["ma" <> ToString[i]], {i, 1, n}];
maBox[n_] := And @@ Table[-5 <= Symbol["ma" <> ToString[i]] <= 5, {i, 1, n}];

gwObj[n_] := With[{v = Table[Symbol["gw" <> ToString[i]], {i, 1, n}]},
    1 - Exp[-2 Total[Table[(v[[i]] - 1.2345)^2, {i, 1, n}]]] + 10^-5 Total[Table[v[[i]]^2, {i, 1, n}]]];
gwVar[n_] := Table[Symbol["gw" <> ToString[i]], {i, 1, n}];
gwBox[n_] := And @@ Table[-5 <= Symbol["gw" <> ToString[i]] <= 5, {i, 1, n}];

pGri5 = griObj[5];   pGri5v = griVar[5];   pGri5b = griBox[5, 15];
pGri10 = griObj[10]; pGri10v = griVar[10]; pGri10b = griBox[10, 600];
pSq10 = sqObj[10];   pSq10v = sqVar[10];   pSq10b = sqBox[10];
pQb3 = qbObj[3];     pQb3v = qbVar[3];     pQb3b = qbBox[3];
pAck3 = ackObj[3];   pAck3v = ackVar[3];   pAck3b = ackBox[3, 32];
pKat8 = katObj[8];   pKat8v = katVar[8];   pKat8b = katBox[8];
pRas5 = rasObj[5];   pRas5v = rasVar[5];   pRas5b = rasBox[5];
pRas8 = rasObj[8];   pRas8v = rasVar[8];   pRas8b = rasBox[8];
pRos10 = rosObj[10]; pRos10v = rosVar[10]; pRos10b = rosBox[10];
pMa10 = maObj[10];   pMa10v = maVar[10];   pMa10b = maBox[10];
pGw10 = gwObj[10];   pGw10v = gwVar[10];   pGw10b = gwBox[10];

(* --- 1-D / 2-D, bare variables (compiled machine path) --- *)
bench["01 Chained trig 1D", NMinimize[{Sin[2 x] + Cos[x], -2 <= x <= 3}, {x}, Method -> R];];
check["01 Chained trig 1D", Round[10^4 First[NMinimize[{Sin[2 x] + Cos[x], -2 <= x <= 3}, {x}, Method -> R]]]];

bench["02 Quartic 1D", NMinimize[{x^4 - 3 x^2 - x, -5 <= x <= 5}, {x}, Method -> R];];
check["02 Quartic 1D", Round[10^4 First[NMinimize[{x^4 - 3 x^2 - x, -5 <= x <= 5}, {x}, Method -> R]]]];

bench["03 Gamma 1D", NMinimize[{Gamma[x], 1 <= x <= 2}, {x}, Method -> R];];
check["03 Gamma 1D", Round[10^5 First[NMinimize[{Gamma[x], 1 <= x <= 2}, {x}, Method -> R]]]];

bench["04 Schwefel 2D", NMinimize[{837.9658 - x Sin[Sqrt[Abs[x]]] - y Sin[Sqrt[Abs[y]]], -500 <= x <= 500 && -500 <= y <= 500}, {x, y}, Method -> R];];
check["04 Schwefel 2D", Round[10^3 First[NMinimize[{837.9658 - x Sin[Sqrt[Abs[x]]] - y Sin[Sqrt[Abs[y]]], -500 <= x <= 500 && -500 <= y <= 500}, {x, y}, Method -> R]]]];

bench["05 Schaffer N2 2D", NMinimize[{0.5 + (Sin[x^2 - y^2]^2 - 0.5)/(1 + 0.001 (x^2 + y^2))^2, -100 <= x <= 100 && -100 <= y <= 100}, {x, y}, Method -> R];];
check["05 Schaffer N2 2D", Round[10^4 First[NMinimize[{0.5 + (Sin[x^2 - y^2]^2 - 0.5)/(1 + 0.001 (x^2 + y^2))^2, -100 <= x <= 100 && -100 <= y <= 100}, {x, y}, Method -> R]]]];

bench["06 Rugged sine 2D", NMinimize[{x^2 + y^2 + 10 Sin[3 x]^2 + 10 Sin[3 y]^2, -5 <= x <= 5 && -5 <= y <= 5}, {x, y}, Method -> R];];
check["06 Rugged sine 2D", Round[10^4 First[NMinimize[{x^2 + y^2 + 10 Sin[3 x]^2 + 10 Sin[3 y]^2, -5 <= x <= 5 && -5 <= y <= 5}, {x, y}, Method -> R]]]];

bench["07 Eggholder 2D LF/20000", NMinimize[{-(y + 47) Sin[Sqrt[Abs[y + x/2 + 47]]] - x Sin[Sqrt[Abs[x - (y + 47)]]], -512 <= x <= 512 && -512 <= y <= 512}, {x, y}, Method -> Regg];];
check["07 Eggholder 2D LF/20000", Round[10^3 First[NMinimize[{-(y + 47) Sin[Sqrt[Abs[y + x/2 + 47]]] - x Sin[Sqrt[Abs[x - (y + 47)]]], -512 <= x <= 512 && -512 <= y <= 512}, {x, y}, Method -> Regg]]]];

(* --- 3-D..10-D, fresh-symbol objectives --- *)
bench["08 Griewank 5D", NMinimize[{Evaluate[pGri5], Evaluate[pGri5b]}, Evaluate[pGri5v], Method -> R];];
check["08 Griewank 5D", Round[10^4 First[NMinimize[{Evaluate[pGri5], Evaluate[pGri5b]}, Evaluate[pGri5v], Method -> R]]]];

bench["09 Griewank 10D", NMinimize[{Evaluate[pGri10], Evaluate[pGri10b]}, Evaluate[pGri10v], Method -> R];];
check["09 Griewank 10D", Round[10^4 First[NMinimize[{Evaluate[pGri10], Evaluate[pGri10b]}, Evaluate[pGri10v], Method -> R]]]];

bench["10 Scaled-quadratic 10D", NMinimize[{Evaluate[pSq10], Evaluate[pSq10b]}, Evaluate[pSq10v], Method -> R];];
check["10 Scaled-quadratic 10D", Round[10^4 First[NMinimize[{Evaluate[pSq10], Evaluate[pSq10b]}, Evaluate[pSq10v], Method -> R]]]];

bench["11 Quadratic bowl 3D", NMinimize[{Evaluate[pQb3], Evaluate[pQb3b]}, Evaluate[pQb3v], Method -> R];];
check["11 Quadratic bowl 3D", Round[10^4 First[NMinimize[{Evaluate[pQb3], Evaluate[pQb3b]}, Evaluate[pQb3v], Method -> R]]]];

bench["12 Ackley 3D", NMinimize[{Evaluate[pAck3], Evaluate[pAck3b]}, Evaluate[pAck3v], Method -> R];];
check["12 Ackley 3D", Round[10^4 First[NMinimize[{Evaluate[pAck3], Evaluate[pAck3b]}, Evaluate[pAck3v], Method -> R]]]];

bench["13 Katsuura 8D", NMinimize[{Evaluate[pKat8], Evaluate[pKat8b]}, Evaluate[pKat8v], Method -> R];];
check["13 Katsuura 8D", Round[10^4 First[NMinimize[{Evaluate[pKat8], Evaluate[pKat8b]}, Evaluate[pKat8v], Method -> R]]]];

bench["14 Rastrigin 5D", NMinimize[{Evaluate[pRas5], Evaluate[pRas5b]}, Evaluate[pRas5v], Method -> R];];
check["14 Rastrigin 5D", Round[10^4 First[NMinimize[{Evaluate[pRas5], Evaluate[pRas5b]}, Evaluate[pRas5v], Method -> R]]]];

bench["15 Rastrigin 8D", NMinimize[{Evaluate[pRas8], Evaluate[pRas8b]}, Evaluate[pRas8v], Method -> R];];
check["15 Rastrigin 8D", Round[10^4 First[NMinimize[{Evaluate[pRas8], Evaluate[pRas8b]}, Evaluate[pRas8v], Method -> R]]]];

bench["16 Rosenbrock 10D", NMinimize[{Evaluate[pRos10], Evaluate[pRos10b]}, Evaluate[pRos10v], Method -> R];];
check["16 Rosenbrock 10D", Round[10^3 First[NMinimize[{Evaluate[pRos10], Evaluate[pRos10b]}, Evaluate[pRos10v], Method -> R]]]];

bench["17 Modified Ackley 10D", NMinimize[{Evaluate[pMa10], Evaluate[pMa10b]}, Evaluate[pMa10v], Method -> R];];
check["17 Modified Ackley 10D", Round[10^4 First[NMinimize[{Evaluate[pMa10], Evaluate[pMa10b]}, Evaluate[pMa10v], Method -> R]]]];

bench["18 Gaussian well 10D LF/40000", NMinimize[{Evaluate[pGw10], Evaluate[pGw10b]}, Evaluate[pGw10v], Method -> Rgw];];
check["18 Gaussian well 10D LF/40000", Round[10^3 First[NMinimize[{Evaluate[pGw10], Evaluate[pGw10b]}, Evaluate[pGw10v], Method -> Rgw]]]];
