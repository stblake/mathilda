(* Holdout A, family C7 (experiment 202), method B: the same SparseArrays as
   202-siam100h-c7-sparse-solve/sparse_solve.m, LinearSolve at WorkingPrecision 80
   (N[A, 80]); the functional of the solution printed to 75 digits.
   Usage: wolframscript -file hc7_ref.wl *)
mat[n_, offs_, d_] := SparseArray[Join[{Band[{1, 1}] -> d},
   Flatten[Table[{Band[{1, 1 + o}] -> 1, Band[{1 + o, 1}] -> 1}, {o, offs}]]], {n, n}];
pow2[n_] := 2^Range[0, Floor[Log2[n - 1]]];
lap2[m_] := With[{n = m^2, r = Table[If[Mod[i, m] == 0, 0, -1], {i, m^2 - 1}]},
   SparseArray[{Band[{1, 1}] -> 4, Band[{1, 2}] -> r, Band[{2, 1}] -> r,
                Band[{1, 1 + m}] -> -1, Band[{1 + m, 1}] -> -1}, {n, n}]];
lap3[m_] := With[{n = m^3,
    r1 = Table[If[Mod[i, m] == 0, 0, -1], {i, m^3 - 1}],
    rm = Table[If[Mod[i, m^2] > m^2 - m || Mod[i, m^2] == 0, 0, -1], {i, m^3 - m}]},
   SparseArray[{Band[{1, 1}] -> 6, Band[{1, 2}] -> r1, Band[{2, 1}] -> r1,
                Band[{1, 1 + m}] -> rm, Band[{1 + m, 1}] -> rm,
                Band[{1, 1 + m^2}] -> -1, Band[{1 + m^2, 1}] -> -1}, {n, n}]];
fib[n_] := Select[Union[Fibonacci[Range[2, 30]]], # < n &];
mats = {
  mat[20000, {1}, Prime[Range[20000]]],
  mat[10000, {1, 2}, Prime[Range[10000]]],
  mat[5000, 3^Range[0, 7], Prime[Range[5000]]],
  mat[3000, fib[3000], Prime[Range[3000]]],
  mat[1000, pow2[1000], Range[1000]],
  mat[3000, pow2[3000], Prime[Range[3000]]],
  mat[1000, pow2[1000], Prime[Range[1000]]],
  lap2[60],
  mat[1000, pow2[1000], Prime[Range[1000]]],
  lap3[20]};
fn = {First, First, First, First, First, Last, #[[2]] &, First, Total, First};
rhs = {1, 1, 1, 1, 1, "n", 1, 1, "ones", 1};
Do[
  With[{a = N[mats[[k]], 80], n = Length[mats[[k]]]},
    b = Which[rhs[[k]] === "n", UnitVector[n, n], rhs[[k]] === "ones", ConstantArray[1, n], True, UnitVector[n, 1]];
    {t, x} = AbsoluteTiming[LinearSolve[a, b]];
    Print["hc7.", k, " ", ToString[N[fn[[k]][x], 75], InputForm], "  (", Round[t, 0.1], " s)"]],
  {k, 1, 10}];
