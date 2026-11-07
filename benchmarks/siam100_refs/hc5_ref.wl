(* hc5 (holdout A, family C5) reference, method B.
   Real problems k = 1..6: the discrete minimax LP  min E  s.t.
   -E <= f_k - p(x_k) <= E  solved EXACTLY by LinearProgramming (simplex over
   the rationals) on the sample points and data rounded to 100 digits
   (Rationalize at 10^-100), so the optimum is E to ~95 digits.
   Complex problems k = 7..10 (conjugation-symmetric data, so the optimum has
   real coefficients): the equal-modulus square system |e_k| = E on the active
   set modulo conjugation (n+2 points; read off hc5_ref.py's SLSQP solve) in
   the n+1 real coefficients and E, by FindRoot at WorkingPrecision 120; then
   every other sample is checked to lie below E.  Method A: hc5_ref.py (scipy
   LP / SLSQP active set + mpmath alternation / KKT Newton at dps 90). *)
rat[x_] := Rationalize[N[x, 130], 10^-100];
realP = <|
  1 -> {Sqrt[#] &, Table[k/128, {k, 0, 128}], 1},
  2 -> {Exp[#] &, Table[Cos[k Pi/64], {k, 0, 64}], 3},
  3 -> {Abs[#] &, Table[-1 + k/64, {k, 0, 128}], 4},
  4 -> {1/Gamma[#] &, Table[2 k/100, {k, 0, 100}], 3},
  5 -> {ArcTan[5 #] &, Table[Cos[k Pi/100], {k, 0, 100}], 5},
  6 -> {Log[1 + #] &, Table[k/64, {k, 0, 64}], 2}|>;
Do[
  {f, xs, n} = realP[k];
  xr = rat /@ xs; fr = rat /@ (f /@ xs);
  V = Table[If[j == 0, 1, xr[[i]]^j], {i, Length[xr]}, {j, 0, n}];
  m = Join[MapThread[Append, {V, ConstantArray[1, Length[xr]]}],
           MapThread[Append, {-V, ConstantArray[1, Length[xr]]}]];
  bb = Join[fr, -fr];
  sol = LinearProgramming[Append[ConstantArray[0, n + 1], 1], m, bb,
                          ConstantArray[{-Infinity, Infinity}, n + 2]];
  Print[k, " ", ToString[NumberForm[N[Last[sol], 75], 70]], "   exact LP (", Head[Last[sol]], ")"],
  {k, 6}];
cplxP = <|
  7 -> {1/Gamma[#] &, 128, 2, {19, 20, 43, 44}, {-0.0230522184895, 1.12741168532919, 0.598999433929804}},
  8 -> {Gamma[# + 2] &, 128, 2, {10, 11, 37, 38}, {1.00181182525426, 0.414294456049011, 0.455481135261200}},
  9 -> {Log[2 + #] &, 64, 2, {0, 19, 20, 32}, {0.692986587777837, 0.500498136335710, -0.143680443443782}},
  10 -> {1/Gamma[#] &, 128, 4, {13, 14, 35, 36, 53, 54},
         {0.000256618973940404, 0.999321680655094, 0.577599648038791, -0.676342057525015, -0.0778734933174255}}|>;
Do[
  {f, M, n, act, c0} = cplxP[k];
  zs = Table[Exp[2 Pi I j/M], {j, 0, M - 1}];
  zn = N[zs, 140]; fn = N[f /@ zs, 140];
  as = Table[Symbol["a" <> ToString[j]], {j, 0, n}];
  e2[j_] := Abs[fn[[j + 1]] - Sum[as[[i + 1]] zn[[j + 1]]^i, {i, 0, n}]]^2;
  r = FindRoot[Table[ComplexExpand[e2[j]] == EE^2, {j, act}],
     Append[Transpose[{as, SetPrecision[c0, 120]}], {EE, SetPrecision[Sqrt[e2[First[act]] /. Thread[as -> c0]], 120]}],
     WorkingPrecision -> 120, PrecisionGoal -> 100, AccuracyGoal -> 100, MaxIterations -> 200];
  others = Max[Table[Sqrt[e2[j]] /. r, {j, Complement[Range[0, M - 1], act, Mod[-act, M]]}]];
  Print[k, " ", ToString[NumberForm[N[EE /. r, 75], 70]], "   max other/E = ", N[others/(EE /. r), 10]],
  {k, 7, 10}];
