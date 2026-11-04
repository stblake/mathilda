(* ::Package:: *)
(* ------------------------------------------------------------------ *)
(* corpus.m -- test cases for the DifferentialForms package (issue     *)
(* #125), run by tools/check_differentialforms.py through              *)
(* tools/differentialforms/run_df_cases.m after the package is loaded. *)
(*                                                                    *)
(* Each record is {id, HoldComplete[expr], kind[, expected[, msgs]]}:  *)
(*   "recorded"  answer compared with the expected answer recorded by  *)
(*               `check_differentialforms.py --record --binary B`      *)
(*   "identity"  expr must evaluate to True (a mathematical identity   *)
(*               the package must satisfy)                             *)
(*   "expect"    answer compared with HoldComplete `expected`, derived *)
(*               by hand -- used for the tensor functions, whose       *)
(*               package rules attach to TensorProduct                 *)
(* The optional 5th element lists messages ("Head::tag") the case is   *)
(* allowed to emit.                                                   *)
(*                                                                    *)
(* Comparison is semantic (see dfSameQ in the driver): the difference  *)
(* must vanish once the coefficients of the package's linear objects   *)
(* are simplified.                                                     *)
(* ------------------------------------------------------------------ *)

(* Shorthands used by identity cases; the package's private heads. *)
dfAltform = DifferentialForms`Private`altform;
g2 = t[x, x] + t[y, y];
g3 = t[x, x] + t[y, y] + t[z, z];
gPolar = t[r, r] + r^2 t[theta, theta];
polar = {x -> r Cos[theta], y -> r Sin[theta]};
cylindrical = {x -> r Cos[theta], y -> r Sin[theta], z -> z};
spherical = {x -> r Cos[theta] Sin[phi], y -> r Sin[theta] Sin[phi], z -> r Cos[phi]};
unitSquare = Chain[{x -> s, y -> t}, {s, 0, 1}, {t, 0, 1}];
dfVanishQ[e_] := dfZeroQ[e];
dfEqualQ[a_, b_] := dfZeroQ[a - b];

$DFCases = {

(* ---------------- issue #125: the reporter's session ---------------- *)
{"issue-d-x",          HoldComplete[d[x]], "recorded"},
{"issue-d-xy",         HoldComplete[d[x, y]], "recorded"},
{"issue-power",        HoldComplete[d[x]^d[y]], "recorded"},
{"issue-wedge-infix",  HoldComplete[d[x] \[Wedge] d[y]], "recorded"},
{"issue-wedge-infix3", HoldComplete[d[x] \[Wedge] d[y] \[Wedge] d[z]], "recorded"},

(* ---------------- d: basis forms ---------------- *)
{"d-xyz",        HoldComplete[d[x, y, z]], "recorded"},
{"d-xx",         HoldComplete[d[x, x]], "recorded"},
{"d-yx",         HoldComplete[d[y, x]], "recorded"},
{"d-zyx",        HoldComplete[d[z, y, x]], "recorded"},
{"d-xyx",        HoldComplete[d[x, y, x]], "recorded"},
{"d-sum",        HoldComplete[d[x] + d[x]], "recorded"},
{"d-scaled",     HoldComplete[3 d[x] + x d[x] - d[y]], "recorded"},
{"d-zero-coef",  HoldComplete[0 d[x]], "recorded"},
{"d-cancel",     HoldComplete[d[x, y] + d[y, x]], "recorded"},

(* ---------------- d: exterior derivative of functions ---------------- *)
{"d-x2",         HoldComplete[d[x^2]], "recorded"},
{"d-x2y",        HoldComplete[d[x^2 y]], "recorded"},
{"d-sin",        HoldComplete[d[Sin[x y]]], "recorded"},
{"d-exp",        HoldComplete[d[Exp[x] y z]], "recorded"},
{"d-log",        HoldComplete[d[Log[x^2 + y^2]]], "recorded"},
{"d-f1",         HoldComplete[d[f[x]]], "recorded"},
{"d-f2",         HoldComplete[d[f[x, y]]], "recorded"},
{"d-f3",         HoldComplete[d[f[x, y, z]]], "recorded"},
{"d-const",      HoldComplete[d[5]], "recorded"},
{"d-rational",   HoldComplete[d[x/y]], "recorded"},
{"d-sqrt",       HoldComplete[d[Sqrt[x^2 + y^2 + z^2]]], "recorded"},
{"d-listable",   HoldComplete[d[{x y, x^2, Sin[z]}]], "recorded"},
{"d-matrix",     HoldComplete[d[{{x, y}, {x y, 1}}]], "recorded"},
{"d-coords-opt", HoldComplete[ExteriorDerivative[f[x, y], Coordinates -> {x, y}]], "recorded"},
{"d-coords-opt3",HoldComplete[ExteriorDerivative[x^2 y z, Coordinates -> {x, y, z}]], "recorded"},
{"d-extder",     HoldComplete[ExteriorDerivative[x y]], "recorded"},
{"d-constant-attr", HoldComplete[SetAttributes[kc, Constant]; d[kc x^2]], "recorded"},
{"d-constant-dt",   HoldComplete[SetAttributes[kc2, Constant]; Dt[kc2 x]], "recorded"},
{"d-constant-msg",  HoldComplete[SetAttributes[{kc3, kc4}, Constant]; d[kc3, kc4]], "recorded"},

(* ---------------- d: exterior derivative of forms ---------------- *)
{"d-1form",      HoldComplete[d[f[x, y] d[x]]], "recorded"},
{"d-1form-b",    HoldComplete[d[f[x, y] d[x] + g[x, y] d[y]]], "recorded"},
{"d-1form-poly", HoldComplete[d[x y d[z] - z^2 d[x]]], "recorded"},
{"d-2form",      HoldComplete[d[x d[y, z]]], "recorded"},
{"d-2form-gen",  HoldComplete[d[f[x, y, z] d[y, z] + g[x, y, z] d[z, x] + h[x, y, z] d[x, y]]], "recorded"},
{"d-3form",      HoldComplete[d[x y z d[x, y, z]]], "recorded"},
{"d-exp-form",   HoldComplete[d[Exp[x] d[y]]], "recorded"},
{"d-at",         HoldComplete[d @ (x d[y])], "recorded"},
{"dd-0",         HoldComplete[d[d[x y z]]], "recorded"},
{"dd-f",         HoldComplete[d[d[f[x, y, z]]]], "recorded"},
{"dd-1form",     HoldComplete[d[d[f[x, y, z] d[x] + g[x, y, z] d[y]]]], "recorded"},

(* ---------------- identities for d ---------------- *)
{"id-dd-0form",  HoldComplete[dfVanishQ[d[d[f[x, y, z, w]]]]], "identity"},
{"id-dd-1form",  HoldComplete[dfVanishQ[d[d[f[x, y, z] d[x] + g[x, y, z] d[y] + h[x, y, z] d[z]]]]], "identity"},
{"id-dd-2form",  HoldComplete[dfVanishQ[d[d[f[x, y, z, w] d[x, y] + g[x, y, z, w] d[z, w]]]]], "identity"},
{"id-dd-poly",   HoldComplete[dfVanishQ[d[d[x^3 y d[z] + Sin[x z] d[y, w]]]]], "identity"},
{"id-leibniz-01",
  HoldComplete[dfEqualQ[d[Wedge[f[x, y, z], g[x, y, z] d[y]]],
    Wedge[d[f[x, y, z]], g[x, y, z] d[y]] + f[x, y, z] d[g[x, y, z] d[y]]]], "identity"},
{"id-leibniz-11",
  HoldComplete[Module[{a = f[x, y, z] d[x], b = g[x, y, z] d[y] + h[x, y, z] d[z]},
    dfEqualQ[d[Wedge[a, b]], Wedge[d[a], b] - Wedge[a, d[b]]]]], "identity"},
{"id-leibniz-12",
  HoldComplete[Module[{a = f[x, y, z, w] d[w], b = g[x, y, z, w] d[x, y]},
    dfEqualQ[d[Wedge[a, b]], Wedge[d[a], b] - Wedge[a, d[b]]]]], "identity"},
{"id-d-linear",
  HoldComplete[dfEqualQ[d[2 x y d[z] + Sin[x] d[y]], 2 d[x y d[z]] + d[Sin[x] d[y]]]], "identity"},
{"id-d-coords-agree",
  HoldComplete[dfEqualQ[d[x^2 Sin[y] z], ExteriorDerivative[x^2 Sin[y] z, Coordinates -> {x, y, z}]]], "identity"},

(* ---------------- ExteriorProduct / Wedge ---------------- *)
{"ep-basic",     HoldComplete[ExteriorProduct[d[x], d[y]]], "recorded"},
{"ep-anti",      HoldComplete[ExteriorProduct[d[y], d[x]]], "recorded"},
{"ep-repeat",    HoldComplete[ExteriorProduct[d[x], d[x]]], "recorded"},
{"ep-scalar",    HoldComplete[ExteriorProduct[f[x], d[x], d[y]]], "recorded"},
{"ep-coef",      HoldComplete[ExteriorProduct[x d[y], y d[z]]], "recorded"},
{"ep-sums",      HoldComplete[ExteriorProduct[a d[x] + b d[y], c d[x] + e d[y]]], "recorded"},
{"ep-sums3",     HoldComplete[ExteriorProduct[d[x] + d[y], d[y] + d[z], d[z] + d[x]]], "recorded"},
{"ep-single",    HoldComplete[ExteriorProduct[d[x, y]]], "recorded"},
{"ep-12",        HoldComplete[ExteriorProduct[d[x], d[y, z]]], "recorded"},
{"ep-21",        HoldComplete[ExteriorProduct[d[y, z], d[x]]], "recorded"},
{"ep-times",     HoldComplete[ExteriorProduct[(x + y) (d[x] + d[z]), d[y]]], "recorded"},
{"wedge-3",      HoldComplete[Wedge[d[x], d[y], d[z]]], "recorded"},
{"wedge-2",      HoldComplete[Wedge[d[x], x d[y]]], "recorded"},
{"wedge-lists",  HoldComplete[Wedge[{d[x], d[y]}, {d[y], d[x]}]], "recorded"},
{"wedge-infix-sum", HoldComplete[(d[x] + d[y]) \[Wedge] (d[x] - d[y])], "recorded"},
{"id-wedge-anticomm",
  HoldComplete[Module[{a = f[x, y, z] d[x] + d[y], b = g[x, y, z] d[z] + x d[x]},
    dfVanishQ[Wedge[a, b] + Wedge[b, a]]]], "identity"},
{"id-wedge-assoc",
  HoldComplete[Module[{a = d[x] + y d[z], b = d[y] - d[x], c = z d[x] + d[z]},
    dfEqualQ[Wedge[Wedge[a, b], c], Wedge[a, Wedge[b, c]]]]], "identity"},
{"id-wedge-graded",
  HoldComplete[Module[{a = d[x] + d[w], b = d[y, z] + x d[x, w]},
    dfEqualQ[Wedge[a, b], Wedge[b, a]]]], "identity"},
{"id-wedge-self", HoldComplete[dfVanishQ[Wedge[x d[y] + d[z], x d[y] + d[z]]]], "identity"},

(* ---------------- predicates, components, bases ---------------- *)
{"formq-1",      HoldComplete[FormQ[d[x]]], "recorded"},
{"formq-0",      HoldComplete[FormQ[x^2]], "recorded"},
{"formq-t",      HoldComplete[FormQ[t[x, y]]], "recorded"},
{"component-1",  HoldComplete[Component[3 x d[x, y], d[x, y]]], "recorded"},
{"component-0",  HoldComplete[Component[3 x d[x, y], d[x, z]]], "recorded"},
{"complist-1",   HoldComplete[ComponentList[3 d[x, y] + 5 d[y, z], Basis -> {x, y, z}]], "recorded"},
{"complist-2",   HoldComplete[ComponentList[x d[x] + y d[y] + z d[z], Basis -> {x, y, z}]], "recorded"},
{"complist-auto",HoldComplete[ComponentList[x d[x] + y d[y]]], "recorded"},
{"basis-1",      HoldComplete[Basis[{x, y, z}, 1]], "recorded"},
{"basis-2",      HoldComplete[Basis[{x, y, z}, 2]], "recorded"},
{"basis-3",      HoldComplete[Basis[{x, y, z}, 3]], "recorded"},
{"basis-0",      HoldComplete[Basis[{x, y, z}, 0]], "recorded"},
{"basis-big",    HoldComplete[Basis[{x, y}, 3]], "recorded"},
{"basis-4-2",    HoldComplete[Basis[{w, x, y, z}, 2]], "recorded"},
{"basis-unsorted", HoldComplete[Basis[{z, x, y}, 2]], "recorded"},
{"standard-1",   HoldComplete[Standard[d[x, y]]], "recorded"},
{"standard-sum", HoldComplete[Standard[3 d[x] + x d[y, z]]], "recorded"},
{"standard-list",HoldComplete[Standard[{d[x], d[y]}]], "recorded"},
{"restrict-usage", HoldComplete[Restrict[(x^2) d[x, y], {x -> a, y -> b}]], "recorded"},
{"restrict-list",  HoldComplete[Restrict[{x d[x], y d[y]}, {x -> 1, y -> 2}]], "recorded"},
{"restrict-vec",   HoldComplete[Restrict[x X[y], {x -> 3}]], "recorded"},

(* ---------------- vectors, linear functionals, Pair ---------------- *)
{"X-1",          HoldComplete[X[x]], "recorded"},
{"X-2",          HoldComplete[X[x, y]], "recorded"},
{"X-lin",        HoldComplete[a X[x] + b X[x] + c X[y]], "recorded"},
{"t-1",          HoldComplete[t[x]], "recorded"},
{"t-2",          HoldComplete[t[x, y]], "recorded"},
{"t-lin",        HoldComplete[t[x, x] + t[x, x] + r^2 t[y, y]], "recorded"},
{"pair-usage-1", HoldComplete[Pair[a X[x] + b X[y], c t[x] + e t[y]]], "recorded"},
{"pair-usage-2", HoldComplete[Pair[a X[x, y], b t[x, y]]], "recorded"},
{"pair-usage-3", HoldComplete[Pair[a X[x, y], b t[x, z]]], "recorded"},
{"pair-usage-4", HoldComplete[Pair[a X[x] + b X[y],
    gxx t[x, x] + gxy t[x, y] + gxy t[y, x] + gyy t[y, y]]], "recorded"},
{"pair-dform",   HoldComplete[Pair[a X[x] + b X[y], f[x, y] d[x] + g[x, y] d[y]]], "recorded"},
{"pair-2form",   HoldComplete[Pair[X[x, y], d[x, y]]], "recorded"},
{"pair-gradient",
  HoldComplete[Pair[InvertMetric[g3], ExteriorDerivative[f[x, y, z]]]], "expect",
  HoldComplete[Derivative[1, 0, 0][f][x, y, z] X[x] + Derivative[0, 1, 0][f][x, y, z] X[y] +
    Derivative[0, 0, 1][f][x, y, z] X[z]]},

(* ---------------- InteriorProduct ---------------- *)
{"ip-1",         HoldComplete[InteriorProduct[X[x], d[x, y]]], "recorded"},
{"ip-2",         HoldComplete[InteriorProduct[X[y], d[x, y]]], "recorded"},
{"ip-0",         HoldComplete[InteriorProduct[X[z], d[x, y]]], "recorded"},
{"ip-1form",     HoldComplete[InteriorProduct[a X[x] + b X[y], f[x, y] d[x] + g[x, y] d[y]]], "recorded"},
{"ip-3form",     HoldComplete[InteriorProduct[X[y], x d[x, y, z]]], "recorded"},
{"ip-vector-sum",HoldComplete[InteriorProduct[a X[x] + b X[y] + c X[z], d[x, y, z]]], "recorded"},
{"id-ip-twice",
  HoldComplete[Module[{v = a X[x] + b X[y] + c X[z]},
    dfVanishQ[InteriorProduct[v, InteriorProduct[v, f[x, y, z] d[x, y] + g[x, y, z] d[y, z]]]]]], "identity"},

(* ---------------- Pullback ---------------- *)
{"pb-usage-polar",   HoldComplete[Pullback[d[x, y], polar]], "recorded"},
{"pb-usage-cyl",     HoldComplete[Pullback[d[x, y, z], cylindrical]], "recorded"},
{"pb-usage-sph",     HoldComplete[Pullback[d[x, y, z], spherical]], "recorded"},
{"pb-usage-generic", HoldComplete[Pullback[f[x, y] d[x] \[Wedge] d[y], {x -> x[u, v], y -> y[u, v]}]], "recorded"},
{"pb-function",      HoldComplete[Pullback[x^2 + y^2, polar]], "recorded"},
{"pb-1form",         HoldComplete[Pullback[x d[y] - y d[x], polar]], "recorded"},
{"pb-list",          HoldComplete[Pullback[{d[x], d[y]}, polar]], "recorded"},
{"pb-curve",         HoldComplete[Pullback[x d[y], {x -> Cos[s], y -> Sin[s]}]], "recorded"},
{"pb-metric-polar",  HoldComplete[Pullback[g2, polar]], "expect",
  HoldComplete[t[r, r] + r^2 t[theta, theta]]},
{"pb-metric-sph",    HoldComplete[Simple[Pullback[g3, spherical]]], "expect",
  HoldComplete[t[r, r] + r^2 t[phi, phi] + r^2 Sin[phi]^2 t[theta, theta]]},
{"id-pb-d-commute-0",
  HoldComplete[dfEqualQ[Pullback[d[f[x, y]], polar], d[Pullback[f[x, y], polar]]]], "identity"},
{"id-pb-d-commute-1",
  HoldComplete[dfEqualQ[Pullback[d[x^2 y d[x] + Sin[y] d[y]], polar],
    d[Pullback[x^2 y d[x] + Sin[y] d[y], polar]]]], "identity"},
{"id-pb-wedge",
  HoldComplete[dfEqualQ[Pullback[Wedge[x d[y], d[x] + d[y]], polar],
    Wedge[Pullback[x d[y], polar], Pullback[d[x] + d[y], polar]]]], "identity"},
{"id-pb-jacobian",
  HoldComplete[dfEqualQ[Pullback[d[x, y, z], spherical],
    Simplify[Det[D[{r Cos[theta] Sin[phi], r Sin[theta] Sin[phi], r Cos[phi]}, {{r, theta, phi}}]]] d[r, theta, phi]]], "identity"},

(* ---------------- Push ---------------- *)
{"push-r",       HoldComplete[Push[X[r], polar]], "recorded"},
{"push-theta",   HoldComplete[Push[X[theta], polar]], "recorded"},
{"push-list",    HoldComplete[Push[{X[r], X[theta]}, polar]], "recorded"},
{"id-push-pair",
  HoldComplete[dfEqualQ[Pair[Push[X[theta], polar], d[x]] /. polar,
    Pair[X[theta], Pullback[d[x], polar]]]], "identity"},

(* ---------------- Chain, Boundary, Integral ---------------- *)
{"chain-usage",      HoldComplete[unitSquare], "recorded"},
{"chain-orient",     HoldComplete[Chain[-1, {x -> s, y -> t}, {s, 0, 1}, {t, 0, 1}]], "recorded"},
{"boundary-square",  HoldComplete[Boundary[unitSquare]], "recorded"},
{"boundary-boundary",HoldComplete[Boundary[Boundary[unitSquare]]], "recorded"},
{"boundary-cube",    HoldComplete[Boundary[Chain[{x -> u, y -> v, z -> w}, {u, 0, 1}, {v, 0, 1}, {w, 0, 1}]]], "recorded"},
{"boundary-bb-cube", HoldComplete[Boundary[Boundary[Chain[{x -> u, y -> v, z -> w}, {u, 0, 1}, {v, 0, 1}, {w, 0, 1}]]]], "recorded"},
{"boundary-disk",    HoldComplete[Boundary[Chain[polar /. r -> s, {s, 0, 1}, {theta, 0, 2 Pi}]]], "recorded"},
{"boundary-segment", HoldComplete[Boundary[Chain[{x -> s}, {s, 0, 1}]]], "recorded"},
{"int-usage-square", HoldComplete[Integral[d[x, y], unitSquare]], "recorded"},
{"int-usage-circle", HoldComplete[SetAttributes[R, Constant];
    Integral[d[x, y], Chain[{x -> r Cos[theta], y -> r Sin[theta]}, {r, 0, R}, {theta, 0, 2 Pi}]]], "recorded"},
{"int-usage-sphere", HoldComplete[Integral[(x) d[y, z] + (-y) d[x, z] + (z) d[x, y],
    Chain[{x -> Cos[theta] Sin[phi], y -> Sin[theta] Sin[phi], z -> Cos[phi]}, {theta, 0, 2 Pi}, {phi, 0, Pi}]]], "recorded"},
{"int-usage-ball",   HoldComplete[SetAttributes[R, Constant];
    Integral[d[x, y, z], Chain[-1, {x -> r Cos[theta] Sin[phi], y -> r Sin[theta] Sin[phi], z -> r Cos[phi]},
      {r, 0, R}, {theta, 0, 2 Pi}, {phi, 0, Pi}]]], "recorded"},
{"int-usage-stokes2d", HoldComplete[Integral[d @ ((x/2) d[y] - (y/2) d[y]), unitSquare] ==
    Integral[((x/2) d[y] - (y/2) d[y]), Boundary @ unitSquare]], "recorded"},
{"int-usage-stokes3d", HoldComplete[SetAttributes[R, Constant];
    Integral[d @ (x d[y, z]), Chain[-1, {x -> r Cos[theta] Sin[phi], y -> r Sin[theta] Sin[phi], z -> r Cos[phi]},
      {r, 0, R}, {theta, 0, 2 Pi}, {phi, 0, Pi}]] ==
    Integral[x d[y, z], Boundary @ Chain[-1, {x -> r Cos[theta] Sin[phi], y -> r Sin[theta] Sin[phi], z -> r Cos[phi]},
      {r, 0, R}, {theta, 0, 2 Pi}, {phi, 0, Pi}]]], "recorded"},
{"int-line",         HoldComplete[Integral[x d[y] - y d[x], Chain[{x -> Cos[s], y -> Sin[s]}, {s, 0, 2 Pi}]]], "recorded"},
{"int-function-0",   HoldComplete[Integral[x^2 y d[x, y], Chain[{x -> u, y -> v}, {u, 0, 2}, {v, 0, 3}]]], "recorded"},
{"int-sum-chains",   HoldComplete[Integral[d[x], Boundary[Chain[{x -> s}, {s, 0, 5}]]]], "recorded"},
{"id-stokes-square",
  HoldComplete[Integral[d[x^2 y d[x] + x y^3 d[y]], unitSquare] ===
    Integral[x^2 y d[x] + x y^3 d[y], Boundary[unitSquare]]], "identity"},
{"id-stokes-disk",
  HoldComplete[Module[{c = Chain[{x -> s Cos[theta], y -> s Sin[theta]}, {s, 0, 1}, {theta, 0, 2 Pi}]},
    Simplify[Integral[d[x d[y] - y d[x]], c] - Integral[x d[y] - y d[x], Boundary[c]]] === 0]], "identity"},
{"id-stokes-cube",
  HoldComplete[Module[{c = Chain[{x -> u, y -> v, z -> w}, {u, 0, 1}, {v, 0, 1}, {w, 0, 1}],
      om = x^2 d[y, z] + y z d[z, x] + x y z d[x, y]},
    Simplify[Integral[d[om], c] - Integral[om, Boundary[c]]] === 0]], "identity"},

(* ---------------- HomotopyOperator ---------------- *)
{"homotopy-usage",   HoldComplete[HomotopyOperator[d[x, y]]], "recorded"},
{"homotopy-1form",   HoldComplete[HomotopyOperator[x d[y] + y d[x]]], "recorded"},
{"homotopy-3form",   HoldComplete[HomotopyOperator[d[x, y, z]]], "recorded"},
{"homotopy-coords",  HoldComplete[HomotopyOperator[d[u, v], Coordinates -> {u, v}]], "recorded"},
{"id-homotopy-1",
  HoldComplete[Module[{om = x^2 y d[z] + z d[x]},
    dfEqualQ[HomotopyOperator[d[om]] + d[HomotopyOperator[om]], om]]], "identity"},
{"id-homotopy-2",
  HoldComplete[Module[{om = x y z d[x, y] + z^2 d[y, z]},
    dfEqualQ[HomotopyOperator[d[om]] + d[HomotopyOperator[om]], om]]], "identity"},
{"id-homotopy-closed",
  HoldComplete[Module[{om = d[x^2 y z + Sin[x]]},
    dfEqualQ[d[HomotopyOperator[om]], om]]], "identity"},

(* ---------------- LieDerivative ---------------- *)
{"lie-usage-fn",     HoldComplete[LieDerivative[a[x, y] X[x] + b[x, y] X[y], f[x, y]]], "recorded"},
{"lie-usage-1form",  HoldComplete[LieDerivative[a[x, y] X[x] + b[x, y] X[y], f[x, y] d[x] + g[x, y] d[y]]], "recorded"},
{"lie-usage-2form",  HoldComplete[LieDerivative[a[x, y] X[x] + b[x, y] X[y], f[x, y] d[x, y]]], "recorded"},
{"lie-usage-vec",    HoldComplete[LieDerivative[a[x, y] X[x] + b[x, y] X[y], f[x, y] X[x] + g[x, y] X[y]]], "recorded"},
{"lie-poly",         HoldComplete[LieDerivative[x X[y] - y X[x], x^2 + y^2]], "recorded"},
{"id-lie-usage-fn",
  HoldComplete[LieDerivative[a[x, y] X[x] + b[x, y] X[y], f[x, y]] ==
    b[x, y]*Derivative[0, 1][f][x, y] + a[x, y]*Derivative[1, 0][f][x, y]], "identity"},
{"id-lie-antisym",
  HoldComplete[Module[{v = x y X[x] + z X[y], u = y^2 X[z] + x X[x]},
    dfVanishQ[LieDerivative[v, u] + LieDerivative[u, v]]]], "identity"},
{"id-lie-jacobi",
  HoldComplete[Module[{u = x y X[x] + z X[y], v = y^2 X[z] + x X[x], w0 = z X[x] + x^2 X[y]},
    dfVanishQ[LieDerivative[u, LieDerivative[v, w0]] + LieDerivative[v, LieDerivative[w0, u]] +
      LieDerivative[w0, LieDerivative[u, v]]]]], "identity"},
{"id-cartan-1form",
  HoldComplete[Module[{v = x y X[x] + z X[y] + X[z], om = x^2 d[y] + y z d[z]},
    dfEqualQ[LieDerivative[v, om], InteriorProduct[v, d[om]] + d[InteriorProduct[v, om]]]]], "identity"},
{"id-cartan-2form",
  HoldComplete[Module[{v = y X[x] - x X[y] + z X[z], om = x d[y, z] + y^2 d[x, y]},
    dfEqualQ[LieDerivative[v, om], InteriorProduct[v, d[om]] + d[InteriorProduct[v, om]]]]], "identity"},
{"id-lie-commutes-d",
  HoldComplete[Module[{v = x y X[x] + z X[y], om = x^2 y d[z]},
    dfEqualQ[LieDerivative[v, d[om]], d[LieDerivative[v, om]]]]], "identity"},

(* ---------------- metrics: Orientation, HodgeStar, cod, Laplace ---------------- *)
{"orient-usage-2",   HoldComplete[Orientation[g2]], "recorded"},
{"orient-usage-cyl", HoldComplete[Orientation[t[r, r] + r^2 t[theta, theta] + t[z, z]]], "recorded"},
{"orient-3",         HoldComplete[Orientation[g3]], "recorded"},
{"orient-polar",     HoldComplete[Orientation[gPolar]], "recorded"},
{"orient-lorentz",   HoldComplete[Orientation[t[x, x] - t[y, y], Riemannian -> False]], "recorded"},
{"orient-scaled",    HoldComplete[Orientation[4 t[x, x] + 9 t[y, y]]], "recorded"},
{"hodge-0",          HoldComplete[HodgeStar[1, g3]], "recorded"},
{"hodge-fn",         HoldComplete[HodgeStar[f[x, y, z], g3]], "recorded"},
{"hodge-dx",         HoldComplete[HodgeStar[d[x], g3]], "recorded"},
{"hodge-dy",         HoldComplete[HodgeStar[d[y], g3]], "recorded"},
{"hodge-dz",         HoldComplete[HodgeStar[d[z], g3]], "recorded"},
{"hodge-dxy",        HoldComplete[HodgeStar[d[x, y], g3]], "recorded"},
{"hodge-dxz",        HoldComplete[HodgeStar[d[x, z], g3]], "recorded"},
{"hodge-dxyz",       HoldComplete[HodgeStar[d[x, y, z], g3]], "recorded"},
{"hodge-sum",        HoldComplete[HodgeStar[f[x, y, z] d[x] + g[x, y, z] d[y], g3]], "recorded"},
{"hodge-2d",         HoldComplete[HodgeStar[d[x], g2]], "recorded"},
{"hodge-polar-dr",   HoldComplete[HodgeStar[d[r], gPolar]], "recorded"},
{"hodge-polar-dth",  HoldComplete[HodgeStar[d[theta], gPolar]], "recorded"},
{"hodge-scaled",     HoldComplete[HodgeStar[d[x], 4 t[x, x] + 9 t[y, y]]], "recorded"},
{"id-hodge-twice-1",
  HoldComplete[dfEqualQ[HodgeStar[HodgeStar[f[x, y, z] d[x] + g[x, y, z] d[z], g3], g3],
    f[x, y, z] d[x] + g[x, y, z] d[z]]], "identity"},
{"id-hodge-twice-2",
  HoldComplete[dfEqualQ[HodgeStar[HodgeStar[x d[x, y] + y d[y, z], g3], g3], x d[x, y] + y d[y, z]]], "identity"},
{"id-hodge-twice-2d",
  HoldComplete[dfEqualQ[HodgeStar[HodgeStar[x d[x] + y d[y], g2], g2], -(x d[x] + y d[y])]], "identity"},
{"id-hodge-volume",
  HoldComplete[dfEqualQ[Wedge[d[x] + 2 d[y], HodgeStar[d[x] + 2 d[y], g3]],
    InnerProduct[d[x] + 2 d[y], d[x] + 2 d[y], g3] d[x, y, z]]], "identity"},
{"cod-usage-1",      HoldComplete[cod[f[x, y, z] d[x], g3]], "recorded"},
{"cod-fn",           HoldComplete[cod[f[x, y, z], g3]], "recorded"},
{"cod-2form",        HoldComplete[cod[f[x, y, z] d[x, y], g3]], "recorded"},
{"cod-div",          HoldComplete[cod[f[x, y, z] d[x] + g[x, y, z] d[y] + h[x, y, z] d[z], g3]], "recorded"},
{"laplace-fn-3",     HoldComplete[Laplace[f[x, y, z], g3]], "recorded"},
{"laplace-fn-2",     HoldComplete[Laplace[f[x, y], g2]], "recorded"},
{"laplace-1form",    HoldComplete[Laplace[f[x, y, z] d[x], g3]], "recorded"},
{"laplace-polar",    HoldComplete[Laplace[f[r, theta], gPolar]], "recorded"},
{"id-laplace-polar",
  HoldComplete[Simplify[PowerExpand[Laplace[f[r, theta], gPolar] +
    (D[f[r, theta], r, r] + D[f[r, theta], r]/r + D[f[r, theta], theta, theta]/r^2)]] === 0], "identity"},
{"id-laplace-sph",
  HoldComplete[Module[{g = t[r, r] + r^2 t[phi, phi] + r^2 Sin[phi]^2 t[theta, theta]},
    Simplify[PowerExpand[Laplace[r^2, g] + 6]] === 0]], "identity"},
{"id-laplace-harmonic",
  HoldComplete[Laplace[x^2 - y^2 + x y z, g3] === 0], "identity"},
{"id-laplace-polar-harmonic",
  HoldComplete[Simplify[PowerExpand[Laplace[Log[r], gPolar]]] === 0], "identity"},

(* ---------------- InnerProduct, InvertMetric, frames ---------------- *)
{"inner-1forms",     HoldComplete[InnerProduct[d[x], d[x], g3]], "recorded"},
{"inner-1forms-0",   HoldComplete[InnerProduct[d[x], d[y], g3]], "recorded"},
{"inner-2forms",     HoldComplete[InnerProduct[d[x, y], d[x, y], g3]], "recorded"},
{"inner-2forms-0",   HoldComplete[InnerProduct[d[x, y], d[x, z], g3]], "recorded"},
{"inner-vec",        HoldComplete[InnerProduct[X[x], X[x], g3]], "recorded"},
{"inner-scalars",    HoldComplete[InnerProduct[a, b, g3]], "recorded"},
{"inner-polar",      HoldComplete[InnerProduct[d[theta], d[theta], gPolar]], "recorded"},
{"inner-sum",        HoldComplete[InnerProduct[a d[x] + b d[y], c d[x] + e d[y], g2]], "recorded"},
{"invmetric-3",      HoldComplete[InvertMetric[g3]], "expect",
  HoldComplete[X[x, x] + X[y, y] + X[z, z]]},
{"invmetric-polar",  HoldComplete[InvertMetric[gPolar]], "expect",
  HoldComplete[X[r, r] + X[theta, theta]/r^2]},
{"invmetric-scaled", HoldComplete[InvertMetric[4 t[x, x] + 9 t[y, y]]], "expect",
  HoldComplete[X[x, x]/4 + X[y, y]/9]},
{"unitframe-3",      HoldComplete[UnitFrame[{x, y, z}, g3]], "recorded"},
{"unitframe-polar",  HoldComplete[UnitFrame[{r, theta}, gPolar]], "recorded"},
{"simple-1",         HoldComplete[Simple[(x + y)^2 - x^2]], "recorded"},

(* ---------------- tensors ---------------- *)
{"tensor-tt",        HoldComplete[Tensor[t[x], t[y]]], "expect", HoldComplete[t[x, y]]},
{"tensor-lists",     HoldComplete[Tensor[{t[x], t[y]}, {t[x], t[y]}]], "expect", HoldComplete[t[x, x] + t[y, y]]},
{"tensorproduct-XX", HoldComplete[TensorProduct[X[x], X[y]]], "expect", HoldComplete[X[x, y]]},
{"tensorproduct-sc", HoldComplete[TensorProduct[a, t[x], t[y]]], "expect", HoldComplete[a t[x, y]]},
{"tensorproduct-lin",HoldComplete[TensorProduct[t[x] + t[y], t[x]]], "expect", HoldComplete[t[x, x] + t[y, x]]},
{"pair-tensor",      HoldComplete[Pair[X[x, y], TensorProduct[t[x], t[y]]]], "expect", HoldComplete[1]},

(* ---------------- display (Format values) ---------------- *)
{"fmt-d-x",          HoldComplete[ToString[d[x]]], "recorded"},
{"fmt-d-xy",         HoldComplete[ToString[d[x, y]]], "recorded"},
{"fmt-sum",          HoldComplete[ToString[3 d[x] + x d[y]]], "recorded"},
{"fmt-input-d-x",    HoldComplete[ToString[d[x], InputForm]], "recorded"},
{"fmt-input-sum",    HoldComplete[ToString[x d[x, y], InputForm]], "recorded"},
{"fmt-X",            HoldComplete[ToString[X[x]]], "recorded"},
{"fmt-t",            HoldComplete[ToString[t[x, y]]], "recorded"},
(* TeXForm applies the package's Format[..., TeXForm] rule; how its
   "{\\wedge }" strings are escaped is not fixed, so check content only. *)
{"fmt-tex",          HoldComplete[StringContainsQ[ToString[TeXForm[d[x, y]]], "wedge"] &&
                       !StringContainsQ[ToString[TeXForm[d[x, y]]], "altform"]], "recorded"},

(* ---------------- usage text ---------------- *)
{"usage-d",          HoldComplete[StringQ[d::usage]], "recorded"},
{"usage-top",        HoldComplete[StringLength[DifferentialForms::usage] > 100], "recorded"},
{"protected-d",      HoldComplete[MemberQ[Attributes[d], Protected]], "recorded"},
{"listable-d",       HoldComplete[MemberQ[Attributes[d], Listable]], "recorded"}
};
