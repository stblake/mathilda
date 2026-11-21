(* Scoping semantics probes: printed results must be identical with and
   without any fast path. *)
p[x_] := Print[InputForm[x]];
g[v_] := Module[{e = 1}, v + e]; p[g[e + 1]];
sh[e_] := Module[{e = 1}, e + 1]; p[sh[99]];
gb[v_] := Block[{e = 1}, v + e]; p[gb[e + 1]];
p[With[{q = 12}, With[{k = q}, k]]];
p[With[{k = 2}, Hold[Table[k, {k, 3}]]]];
gs[f_] := Sum[f, {k, 1, 3}]; p[gs[1/k]];
p[Function[x, (x + #) &][a]];
p[Module[{c = 0}, Scan[(c = c + 1) &, {1, 2, 3}]; c]];
p[Module[{x}, Hold[x]] /. s_Symbol :> StringMatchQ[SymbolName[s], "x$*"]];
p[With[{x = 5}, Function[x, x + 1][2]]];
p[With[{x = 5}, Function[{y}, x + y][2]]];
p[With[{x = 5}, Module[{y = x}, y + x]]];
p[Module[{a = 1}, With[{b = a}, Module[{a = 2}, a + b]]]];
big = Hold @@ {Table[h[k, u, v], {k, 200}]};
f1[x_] := Module[{t = x}, If[t > 0, t, big; ReleaseHold[big]; u + t]];
p[f1[1]]; p[Length[f1[-1]]];
f2[x_] := With[{y = x}, Hold[{y, z, big}]]; p[f2[3][[1, 1 ;; 2]]];
(* a local captured by a value *)
cap[v_] := Module[{t}, t = v; {t, Hold[v]}]; p[cap[t]];
(* nested function shadowing *)
p[Function[{x}, Function[{x}, x][2] + x][10]];
(* With into held + nested scoping that rebinds *)
p[With[{a = 1}, Hold[With[{a = 2}, a], a]]];
p[With[{a = 1}, Hold[Module[{a}, a], a]]];
(* repeated application keeps working after cached summaries *)
k1[x_] := Module[{r = x}, Do[r = r + 1, {3}]; If[r > 100, big, r]];
p[Table[k1[i], {i, 5}]];
(* redefinition: summaries on the old RHS must not affect the new one *)
k2[x_] := Module[{r = x}, {r, w}]; p[k2[1]]; k2[x_] := Module[{w = x}, {w, r}]; p[k2[1] /. s_Symbol :> SymbolName[s]];
(* system-named local *)
p[With[{Plus = Times}, Plus[2, 3]]];
p[Module[{List = 7}, List]];
