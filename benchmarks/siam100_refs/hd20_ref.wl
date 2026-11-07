(* Reference B for holdout family hd20 (exp 225): regular half-iterates via the
   truncated SCHROEDER system (the statement asked of the systems), in exact
   (or 150-digit) arithmetic, at several truncation orders n, to show the
   truncation error and to confirm Reference A (siam100_refs/hd20_ref.py:
   closed forms / Koenigs limits, no power series).
     sigma_n(z) = Sum_{k=1}^n c_k z^k,  c_1 = 1,
     [z^j] (sigma_n(f(z)) - lam sigma_n(z)) = 0,  j = 2..n,
     h_n(z0) = the root x of sigma_n(x) == Sqrt[lam] sigma_n(z0).
   Run: wolframscript -file hd20_ref.wl *)
schroederH[f_, lam_, z0_, n_, prec_] := Module[{p1, pw, A, c, s},
  (* powers of the truncated Taylor polynomial of f (exact, or at prec digits) *)
  p1 = Normal[Series[f, {z, 0, n}]]; If[prec =!= Infinity, p1 = N[p1, prec]];
  pw = NestList[Expand[Normal[Series[# p1, {z, 0, n}]]] &, p1, n - 1];
  pw = CoefficientList[#, z, n + 1] & /@ pw;
  A = Table[If[j == 1, KroneckerDelta[k, 1], pw[[k, j + 1]] - lam KroneckerDelta[j, k]], {j, n}, {k, n}];
  c = LinearSolve[If[prec === Infinity, A, N[A, prec]], UnitVector[n, 1]];
  s[x_] := c . x^Range[n];
  x /. FindRoot[s[x] == Sqrt[lam] s[z0], {x, z0}, WorkingPrecision -> If[prec === Infinity, 150, prec - 10]]];
probs = {
  {"hd20.1", 16 z/(4 + 3 z), 4, 1/2, 0},
  {"hd20.2", Tan[4 ArcTan[z]], 4, 1/10, 0},
  {"hd20.3", 4 z (1 - z), 4, 1/10, 0},
  {"hd20.4", 4 z (1 + z), 4, 1/10, 0},
  {"hd20.5", 1 - (4 (1 - z)^3 - 3 (1 - z)), 9, 1/10, 1},   (* u = 1 - x; answer 1 - u *)
  {"hd20.6", 2 z/(1 + z^2), 2, 1/5, 0},
  {"hd20.7", 2 Sin[z], 2, 1/5, 0},
  {"hd20.8", 2 z Exp[z], 2, 1/10, 0},
  {"hd20.9", Sqrt[2]^(z + 2) - 2, Log[2], -1/2, 2}};          (* w = x - 2; answer 2 + w *)
closed = {8/9, 20/99, Sin[Sqrt[2] ArcSin[Sqrt[1/10]]]^2, Sinh[Sqrt[2] ArcSinh[Sqrt[1/10]]]^2,
  Cos[Sqrt[3] ArcCos[9/10]], Tanh[Sqrt[2] ArcTanh[1/5]], None, None, None};
Do[With[{p = probs[[i]]},
   Module[{vals, post},
    post[v_] := Which[p[[5]] === 1, 1 - v, p[[5]] === 2, 2 + v, True, v];
    vals = Table[{n, N[post[schroederH[p[[2]], p[[3]], p[[4]], n, If[i == 9, 150, Infinity]]], 75]},
       {n, If[i == 9, {25, 30, 90, 110, 140}, {20, 25, 30, 60, 90, 110, 130}]}];
    Print[p[[1]]];
    Do[Print["  n=", v[[1]], "  ", v[[2]]], {v, vals}];
    If[closed[[i]] =!= None, Print["  closed  ", N[closed[[i]], 75],
       "   |h_maxn - closed| = ", N[Abs[vals[[-1, 2]] - closed[[i]]], 3]]]]], {i, If[Length[$ScriptCommandLine] > 1, ToExpression /@ Rest[$ScriptCommandLine], Range[Length[probs]]]}];
robbins[n_, b_] := Module[{A, v},
  A = Table[(k Log[b])^j/j! - KroneckerDelta[j, k], {j, 0, n - 1}, {k, 1, n}];
  v = LinearSolve[N[A, 200], UnitVector[n, 1]];
  x /. FindRoot[-1 + v . x^Range[n] == -1/2, {x, 1/2}, WorkingPrecision -> 150]];
Print["hd20.10 robbins base 2 n=16 (LinearSolve at 200 digits): ", N[robbins[16, 2], 75]];
