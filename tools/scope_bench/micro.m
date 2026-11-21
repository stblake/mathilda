(* Scoping micro-benchmarks. Each builds a function whose body holds a large
   dead branch (If[cond, live, dead]: If is HoldRest, so the dead branch is
   never evaluated but IS walked by the scoping substitution), then times
   n calls. Reported: microseconds per call. *)
dead[m_] := StringJoin[Table["q" <> ToString[k] <> "[u, v, {1, 2, 3}, w + " <> ToString[k] <> "]; ", {k, m}]];
time[expr_, n_] := Module[{best = Infinity, t0}, Do[t0 = TimeUsed[]; Do[expr, {n}]; best = Min[best, TimeUsed[] - t0], {3}]; Round[10^6 best/n, 0.1]];
SetAttributes[time, HoldFirst];
n = 2000;
Do[
  ToExpression["fM[x_] := Module[{a = x, b}, If[x > 0, b = a + 1, " <> dead[m] <> "Null]; b]"];
  ToExpression["fW[x_] := With[{a = x}, If[x > 0, a + 1, " <> dead[m] <> "Null]]"];
  ToExpression["fB[x_] := Block[{a = x, b}, If[x > 0, b = a + 1, " <> dead[m] <> "Null]; b]"];
  ToExpression["fF = Function[{a}, If[a > 0, a + 1, " <> dead[m] <> "Null]]"];
  ToExpression["fP[x_] := If[x > 0, x + 1, " <> dead[m] <> "Null]"];
  Print["deadstmts=", m, " bodyLeaves~", LeafCount[ToExpression["Hold[" <> dead[m] <> "Null]"]],
    "  rule-only ", time[fP[1], n], "  Module ", time[fM[1], n], "  With ", time[fW[1], n],
    "  Block ", time[fB[1], n], "  Function ", time[fF[1], n], "  us/call"],
  {m, {0, 10, 100, 1000}}];
(* nesting depth: d nested Modules around a body with a dead branch of 100 statements *)
Do[
  inner = "If[a1 > 0, a1 + 1, " <> dead[100] <> "Null]";
  s = inner; Do[s = "Module[{a" <> ToString[j] <> " = " <> If[j == 1, "x", "a" <> ToString[j - 1]] <> "}, " <> s <> "]", {j, d, 1, -1}];
  ToExpression["fN[x_] := " <> s];
  Print["nest depth=", d, "  ", time[fN[1], n], " us/call"], {d, {1, 2, 4, 8}}];
(* many locals: L locals, all referenced once, dead branch of 100 statements *)
Do[
  locs = StringRiffle[Table["l" <> ToString[j], {j, L}], ", "];
  ToExpression["fL[x_] := Module[{" <> locs <> "}, l1 = x; If[x > 0, l1 + 1, " <> dead[100] <> locs <> "]]"];
  Print["locals=", L, "  ", time[fL[1], n], " us/call"], {L, {1, 4, 16, 32}}];
