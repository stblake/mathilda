(* Holdout-A D10 family (experiment 215), method B (generalises d10_ref.wl):
   lambda is the constant in the Poisson equation of the ratio chain
       h(r) - Sum_i p_i h(c_i + s_i/r) = log r - lambda    on J = [lo, hi],
   solved by Chebyshev collocation at N first-kind nodes,
   h = Sum_{k>=1} a_k T_k((r - mid)/half), at WorkingPrecision wp.
   Usage: wolframscript -file hd10_ref.wl  (prints N, exp(lambda) per case) *)
cases = <|
  2 -> {{{3, -1, 1/2}, {3, 1, 1/2}}, {2, 4}},
  3 -> {{{4, -1, 1/2}, {4, 1, 1/2}}, {3, 5}},
  4 -> {{{5/2, -1, 1/2}, {5/2, 1, 1/2}}, {3/2, 7/2}},
  5 -> {{{3, -2, 1/2}, {3, 2, 1/2}}, {2, 4}},
  6 -> {{{1, 1, 1/2}, {2, 1, 1/2}}, {1, 3}},
  7 -> {{{3, -1, 1/3}, {3, 0, 1/3}, {3, 1, 1/3}}, {2, 4}},
  8 -> {{{3, -1, 1/4}, {3, 1, 3/4}}, {2, 4}},
  9 -> {{{1, 1, 1/3}, {2, 1, 1/3}, {3, 1, 1/3}}, {1, 4}},
  10 -> {{{2, -1, 1/4}, {2, 1, 3/4}}, {1, 3}}|>;
lam[{maps_, {lo_, hi_}}, nn_, wp_] := Module[{mid = (lo + hi)/2, half = (hi - lo)/2, xs, row, A, b},
  xs = Table[N[Cos[Pi (j - 1/2)/nn], wp], {j, nn}];
  row[z_] := Table[ChebyshevT[k, z], {k, 1, nn - 1}];
  A = Table[With[{r = mid + half x},
      Append[row[x] - Sum[m[[3]] row[(m[[1]] + m[[2]]/r - mid)/half], {m, maps}], 1]], {x, xs}];
  b = Log[mid + half xs];
  Exp[Last[LinearSolve[A, b]]]];
ns = ToExpression[If[Length[$ScriptCommandLine] > 1, $ScriptCommandLine[[2]], "{60, 90}"]];
ks = ToExpression[If[Length[$ScriptCommandLine] > 2, $ScriptCommandLine[[3]], "Range[2, 9]"]];
Do[Do[Print[k, "  N=", nn, "  ", ToString[N[lam[cases[k], nn, 110], 66], InputForm]], {nn, ns}], {k, ks}];
