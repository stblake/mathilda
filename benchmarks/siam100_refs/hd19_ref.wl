(* Reference B for holdout family hd19 (exp 224): an independent DUALITY
   CERTIFICATE (as siam100_refs/d19_ref.wl, generalised).  With element
   energies e_e(v) = v_e^T K_e v_e (unit conductivity),
       J(c) = max_v [2 F.v - Sum_e c_e e_e(v)],
   so for every v:  min_c J >= LB(v) = 2 F.v - max_{c feasible} Sum c_e e_e(v)
   (a fractional knapsack: hi on the largest energies, lo on the smallest).
   UB = J(c_star) by exact-assembly LinearSolve at 100 digits for the design c_star
   from reference A (hd19_ref.py -> hd19_cstar_<k>.txt); LB at v = U(c_star).
   LB == UB to the printed digits proves the minimum, independently of how c_star
   was found.  Run: wolframscript -file hd19_ref.wl k *)
k = ToExpression[Last[$ScriptCommandLine]];
probs = <|1 -> {2, 4, 4, 1, 10, 11/2, 1}, 2 -> {2, 4, 4, 1, 100, 103/4, 1}, 3 -> {2, 6, 6, 1, 100, 101/2, 1},
  4 -> {2, 5, 5, 1, 100, 101/2, 1}, 5 -> {1, 8, 0, 1, 100, 101/2, 1}, 6 -> {1, 16, 0, 1, 4, 5/2, 1},
  7 -> {2, 4, 4, 1, 100, 101/2, xx}, 8 -> {2, 8, 4, 1, 100, 101/2, 1}, 9 -> {2, 4, 4, 1, 2, 3/2, 1},
  10 -> {2, 6, 6, 1, 1000, 1003/4, 1}|>;
{dim, nx, ny, lo, hi, mean, load} = probs[k];
prec = 100;
cstar = SetPrecision[ToExpression[StringReplace[#, "e" -> "*^"]] & /@ StringSplit[Import[
     FileNameJoin[{DirectoryName[$InputFileName], "hd19_cstar_" <> ToString[k] <> ".txt"}], "Text"]], prec];
If[dim == 1,
  h = 1/nx; nIdx[i_] := If[1 <= i <= nx - 1, i, 0];
  elems = Table[{nIdx[a], nIdx[a + 1]}, {a, 0, nx - 1}]; nInt = nx - 1;
  F = Table[h (load /. xx -> h i), {i, nInt}]; loc = {{1, -1}, {-1, 1}}/h,
  h = 1/ny; nIdx[i_, j_] := If[1 <= i <= nx - 1 && 1 <= j <= ny - 1, (i - 1) (ny - 1) + j, 0];
  elems = Flatten[Table[{nIdx[a, b], nIdx[a + 1, b], nIdx[a + 1, b + 1], nIdx[a, b + 1]}, {a, 0, nx - 1}, {b, 0, ny - 1}], 1];
  nInt = (nx - 1) (ny - 1);
  F = Flatten[Table[h^2 (load /. xx -> h i), {i, nx - 1}, {j, ny - 1}]];
  loc = {{4, -1, -2, -1}, {-1, 4, -1, -2}, {-2, -1, 4, -1}, {-1, -2, -1, 4}}/6];
ne = Length[elems];
stiff[c_] := Module[{K = ConstantArray[0, {nInt, nInt}]},
  Do[With[{idx = elems[[e]]}, Do[If[idx[[r]] > 0 && idx[[s]] > 0,
       K[[idx[[r]], idx[[s]]]] += c[[e]] loc[[r, s]]], {r, Length[idx]}, {s, Length[idx]}]], {e, ne}]; K];
U = LinearSolve[stiff[cstar], F];
UB = F . U;
en = Table[With[{ve = Table[If[i > 0, U[[i]], 0], {i, elems[[e]]}]}, ve . loc . ve], {e, ne}];
ord = Ordering[en, All, Greater];
budget = mean ne - lo ne; cmax = ConstantArray[lo, ne];
Do[With[{add = Min[hi - lo, budget]}, cmax[[e]] += add; budget -= add], {e, ord}];
LB = 2 F . U - cmax . en;
Print["hd19.", k, " feasible: ", {Min[cstar] >= lo, Max[cstar] <= hi, N[Total[cstar] - mean ne, 3]}];
Print["hd19.", k, " UB = ", N[UB, 72]];
Print["hd19.", k, " LB = ", N[LB, 72]];
Print["hd19.", k, " UB - LB = ", N[UB - LB, 5]];
