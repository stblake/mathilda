(* prelude_selftest.m -- the PDE corpus verifier judged on hand-made answers
 * with a known verdict.  Run:  ./Mathilda -file DSolve_test_status/pde/prelude_selftest.m
 * (from the repo root); prints one line per case and "selftest: N/N".  Each case
 * pins one verdict path of dsolve_pde_corpus_prelude.m, including the wrong
 * answers it must REFUSE. *)
Get["DSolve_test_status/pde/dsolve_pde_corpus_prelude.m"];

pdCases = {
  (* explicit, arbitrary function symbolic *)
  {"transport general", D[u[x, t], t] + D[u[x, t], x] == 0, u[x, t], {x, t},
     {{u[x, t] -> C[1][t - x]}}, 0},
  {"transport sign flipped", D[u[x, t], t] + D[u[x, t], x] == 0, u[x, t], {x, t},
     {{u[x, t] -> C[1][t + x]}}, 1},
  {"wave two functions", D[u[x, t], {t, 2}] == 4 D[u[x, t], {x, 2}], u, {x, t},
     {{u -> Function[{x, t}, C[1][x - 2 t] + C[2][x + 2 t]]}}, 0},
  {"wave wrong speed", D[u[x, t], {t, 2}] == 4 D[u[x, t], {x, 2}], u, {x, t},
     {{u -> Function[{x, t}, C[1][x - 3 t]]}}, 1},
  {"3-var invariants", D[w[x, y, z], x] + D[w[x, y, z], y] + D[w[x, y, z], z] == 0,
     w[x, y, z], {x, y, z}, {{w[x, y, z] -> C[1][y - x, z - x]}}, 0},
  (* IVPs *)
  {"transport IVP", {D[u[t, x], t] - 3 D[u[t, x], x] == 0, u[0, x] == Exp[-x^2]},
     u[t, x], {t, x}, {{u[t, x] -> Exp[-(x + 3 t)^2]}}, 0},
  {"transport IVP wrong data", {D[u[t, x], t] - 3 D[u[t, x], x] == 0, u[0, x] == Exp[-x^2]},
     u[t, x], {t, x}, {{u[t, x] -> Exp[-(x + 3 t)^2] + t}}, 1},
  {"transport IVP unfitted", {D[u[t, x], t] - 3 D[u[t, x], x] == 0, u[0, x] == Exp[-x^2]},
     u[t, x], {t, x}, {{u[t, x] -> C[1][x + 3 t]}}, 7},
  {"heat Gaussian IVP", {D[u[x, t], t] == D[u[x, t], {x, 2}], u[x, 0] == Exp[-x^2]},
     u[x, t], {x, t}, {{u[x, t] -> Exp[-x^2/(1 + 4 t)]/Sqrt[1 + 4 t]}}, 0},
  (* answers as Mathilda's own methods spell them *)
  {"d'Alembert with an unevaluated integral",
     {D[u[x, t], {t, 2}] == 4 D[u[x, t], {x, 2}], u[x, 0] == f[x], Derivative[0, 1][u][x, 0] == g[x]},
     u, {x, t},
     {{u -> Function[{x, t}, 1/2 (f[-2 t + x] + f[2 t + x]) + 1/4 Integrate[g[K], {K, -2 t + x, 2 t + x}]]}}, 0},
  {"heat kernel: distributional initial condition",
     {D[u[x, t], t] == D[u[x, t], {x, 2}], u[x, 0] == f[x]}, u[x, t], {x, t},
     {{u[x, t] -> Integrate[f[K]/E^((-K + x)^2/(4 t)), {K, -Infinity, Infinity}]/(2 Sqrt[Pi t])}}, 6},
  {"Charpit complete integral (bare constants)", D[u[x, y], x] D[u[x, y], y] == 1, u[x, y], {x, y},
     {{u[x, y] -> y/C[1] + x C[1] + C[2]}}, 0},
  (* implicit *)
  {"inviscid Burgers implicit", u[x, y] D[u[x, y], x] + D[u[x, y], y] == 0, u[x, y], {x, y},
     {{-(x - y u[x, y])/u[x, y] == C[1][u[x, y]]}}, 4},
  (* series *)
  {"heat Dirichlet series",
     {D[u[x, t], t] == D[u[x, t], {x, 2}], u[0, t] == 0, u[Pi, t] == 0, u[x, 0] == x (Pi - x)},
     u[x, t], {x, t},
     {{u[x, t] -> Inactive[Sum][(4 (1 - (-1)^K[1]))/(Pi K[1]^3) E^(-K[1]^2 t) Sin[K[1] x],
                              {K[1], 1, Infinity}]}}, 5},
  {"heat Dirichlet series wrong coefficient",
     {D[u[x, t], t] == D[u[x, t], {x, 2}], u[0, t] == 0, u[Pi, t] == 0, u[x, 0] == x (Pi - x)},
     u[x, t], {x, t},
     {{u[x, t] -> Inactive[Sum][(2 (1 - (-1)^K[1]))/(Pi K[1]^3) E^(-K[1]^2 t) Sin[K[1] x],
                              {K[1], 1, Infinity}]}}, 1},
  {"heat series wrong decay",
     {D[u[x, t], t] == D[u[x, t], {x, 2}], u[0, t] == 0, u[Pi, t] == 0, u[x, 0] == x (Pi - x)},
     u[x, t], {x, t},
     {{u[x, t] -> Inactive[Sum][(4 (1 - (-1)^K[1]))/(Pi K[1]^3) E^(-K[1] t) Sin[K[1] x],
                              {K[1], 1, Infinity}]}}, 1}
};

pdOK = 0;
(* pd-prefixed iterator: the cases spell c, f, g, ... as their own symbols *)
Do[Module[{got = pdeCheckAnswer[pdCase[[2]], pdCase[[3]], pdCase[[4]], pdCase[[5]]]},
    If[got === pdCase[[6]], pdOK++];
    Print[If[got === pdCase[[6]], "ok  ", "BAD "], pdCase[[1]], ": got ", pdeVerdictName[got],
          ", want ", pdeVerdictName[pdCase[[6]]]]],
  {pdCase, pdCases}];
Print["selftest: ", pdOK, "/", Length[pdCases]];
