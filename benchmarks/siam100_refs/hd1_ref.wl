(* Holdout A, family hd1: the same digit-automaton power-sum recursion as
   hd1_ref.py, independently coded in WL: exact rational moments of the explicit
   part (numbers with <= k0 digits enumerated by Select over Range), the moment
   matrix in 160-digit arithmetic, LinearSolve for the closed-form tail.
   Run: wolframscript -file hd1_ref.wl *)
kmp[pat_List] := Module[{m = Length[pat]},
  {Range[0, m - 1], 0, Function[{q, d}, Module[{s = Append[Take[pat, q], d], L},
     L = SelectFirst[Range[Min[Length[s], m], 0, -1],
        # == 0 || Take[s, -#] === Take[pat, #] &];
     If[L == m, None, L]]]}];
nod[set_] := {{0}, 0, Function[{q, d}, If[MemberQ[set, d], None, 0]]};
nondec = {Range[0, 9], 0, Function[{q, d}, If[d >= Max[q, 1], d, None]]};
probs = <|1 -> {10, nod[{9}], 1}, 2 -> {10, nod[{0}], 1}, 3 -> {10, kmp[{0, 0}], 1},
   4 -> {10, kmp[{3, 1, 4}], 1}, 5 -> {2, kmp[{1, 1}], 1}, 6 -> {3, nod[{2}], 1},
   7 -> {10, nod[{9}], 2}, 8 -> {10, nod[{0, 2, 4, 6, 8}], 1}, 9 -> {10, nondec, 1},
   10 -> {10, kmp[{2, 0, 2, 6}], 1}|>;
run[{b_, {states_, q0_, delta_}, s_}, k0_, jj_, prec_] :=
 Module[{st, nums, cur, tot = 0, psi, ns = Length[states], idx, M, v, y, nxt},
  idx = AssociationThread[states -> Range[ns]];
  cur = Select[Table[{d, delta[q0, d]}, {d, 1, b - 1}], #[[2]] =!= None &];
  Do[tot += Total[N[1/cur[[All, 1]]^s, prec]];
     If[k < k0, cur = Flatten[Table[With[{q2 = delta[c[[2]], d]},
          If[q2 === None, Nothing, {b c[[1]] + d, q2}]], {c, cur}, {d, 0, b - 1}], 1]],
   {k, 1, k0}];
  psi = Table[0, {ns}, {jj + 1}];
  Do[psi[[idx[c[[2]]], j + 1]] += N[1/c[[1]]^(s + j), prec], {c, cur}, {j, 0, jj}];
  M = Table[0, {ns (jj + 1)}, {ns (jj + 1)}];
  Do[With[{q2 = delta[q, d]}, If[q2 =!= None,
      Do[M[[(idx[q2] - 1) (jj + 1) + j + 1, (idx[q] - 1) (jj + 1) + j + w + 1]] +=
          N[Binomial[-(s + j), w] If[w == 0, 1, d^w] / b^(s + j + w), prec], {j, 0, jj}, {w, 0, jj - j}]]],
   {q, states}, {d, 0, b - 1}];
  v = N[Flatten[psi], prec];
  y = LinearSolve[IdentityMatrix[ns (jj + 1)] - M, M . v];
  N[tot, prec] + Sum[y[[(i - 1) (jj + 1) + 1]], {i, ns}]];
Do[With[{p = probs[k]}, Print[k, " ",
   NumberForm[run[p, If[p[[1]] == 10, 3, If[p[[1]] == 2, 14, 9]], 60, 160], 75]]], {k, 1, 10}];
