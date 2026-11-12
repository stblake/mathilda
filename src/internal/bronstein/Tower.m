(* ::Package:: *)

(* ========================================================================== *)
(* Tower.m -- the surface integrand as a function on a differential tower.    *)
(*                                                                            *)
(* brTower[f, x] replaces the Log/Exp kernels of f, innermost first, by        *)
(* symbols t_1, t_2, ... (monomials over Q(x)), detects the radical, and      *)
(* classifies the integrand ([B] Part 2 needs K = k(theta)(y), theta the top  *)
(* monomial and y^n in k(theta); BRONSTEIN_PLAN.md "Front end").              *)
(*                                                                            *)
(* Every rewrite is an IDENTITY on every branch, because the answer is        *)
(* verified against the user's integrand:                                     *)
(*   - an exponential reducible over the tower, E^b = E^c Prod t_i^r_i        *)
(*     Prod a_j^s_j (Risch`ExpReducible, r_i, s_j integers, c constant), is    *)
(*     rewritten (E^(u+v) = E^u E^v and E^(n Log a) = a^n hold everywhere);   *)
(*     a fractional r_i re-bases t_i = u^q, u = E^(b_i/q), the same identity;  *)
(*   - a REDUCIBLE logarithm (Log[x^2] next to Log[x]) is declined: Log[x^2]  *)
(*     = 2 Log[x] fails off the positive axis, so the rewrite is not one.     *)
(*   - Log[a]/E^b with a radical inside is admitted only as the single       *)
(*     monomial of Shape B's curve-monomial case ("BCurve"); every head other *)
(*     than rational operations, Log, Exp and rational powers is declined.    *)
(*                                                                            *)
(* The record (an Association):                                               *)
(*   "x"        the integration variable                                      *)
(*   "monos"    the monomial Associations (sym, kind, arg, D, kernel), in     *)
(*              tower order                                                   *)
(*   "syms"     {t_1, ..., t_r}, in tower order (theta last for Shape A)      *)
(*   "kinds"    {"Log" | "Exp", ...}                                          *)
(*   "args"     {a_i}: t_i = Log[a_i] or E^a_i, a_i in the symbols before t_i *)
(*   "rules"    {x -> 1, t_i -> D t_i}: the derivation, for Risch`Derivation  *)
(*   "back"     {t_i -> surface kernel}                                       *)
(*   "kernels"  the surface kernels, in tower order                           *)
(*   "f"        the integrand in the symbols (radicals still as Power)        *)
(*   "radical"  <|"b" -> b, "n" -> n|> in the symbols, or None                *)
(*   "shape"    "Rational" | "Transcendental" | "Algebraic" | "A" | "B" |     *)
(*              "BCurve" (B: radicand in Q(x), kernel arguments free of the   *)
(*              radical; BCurve: one kernel whose argument contains it)       *)
(*   "theta"    the top monomial (Shape A), else None                         *)
(*   "kind"     "Log" | "Exp" for theta, else None                            *)
(*   "lower"    the monomials below theta (Shape A), else "syms"              *)
(* ========================================================================== *)

brUnsupportedHeads::usage = "brUnsupportedHeads[] is the pattern of heads the Bronstein front end declines (trigonometric, hyperbolic, inverse and special functions, Abs/Sign, derivatives, integrals, and approximate or complex numbers). Phase 12 may add a TrigToExp front end.";
brUnsupportedHeads[] := Sin | Cos | Tan | Cot | Sec | Csc | ArcSin | ArcCos | ArcTan | ArcCot |
  ArcSec | ArcCsc | Sinh | Cosh | Tanh | Coth | Sech | Csch | ArcSinh | ArcCosh | ArcTanh |
  ArcCoth | ArcSech | ArcCsch | Erf | Erfc | Erfi | ExpIntegralEi | LogIntegral | PolyLog |
  Gamma | Abs | Sign | Derivative | Integrate | RootSum | Root | Piecewise | Floor | Ceiling;

brCheckHeads::usage = "brCheckHeads[f, x] throws a Decline when f contains a head outside the scope (brUnsupportedHeads), an approximate number, or a complex constant; otherwise returns True.";
brCheckHeads[f_, x_Symbol] := (
  If[!FreeQ[f, brUnsupportedHeads[]], throwDecline["unsupported function in the integrand"]];
  If[!FreeQ[f, _Real], throwDecline["approximate number in the integrand"]];
  If[!FreeQ[f, _Complex], throwDecline["complex constant in the integrand"]];
  True);

brKernelPattern::usage = "brKernelPattern[] matches a Log or Exp kernel, Log[a] or Power[E, b].";
brKernelPattern[] := Log[_] | Power[E, _];

brDependsQ::usage = "brDependsQ[e, vars] is True when e involves one of the symbols vars.";
brDependsQ[e_, vars_List] := !FreeQ[e, Alternatives @@ vars];

brInnermostKernels::usage = "brInnermostKernels[f, vars] lists, in canonical order, the Log/Exp kernels of f whose argument depends on vars and contains no further Log/Exp kernel.";
brInnermostKernels[f_, vars_List] := Sort[DeleteDuplicates[Cases[f,
  k : (Log[a_] | Power[E, a_]) /; brDependsQ[a, vars] && FreeQ[a, brKernelPattern[]], {0, Infinity}]]];

brRadicalQ::usage = "brRadicalQ[e, vars] is True when e contains a rational power b^r (r non-integer) of a base depending on vars.";
brRadicalQ[e_, vars_List] := !FreeQ[e, Power[b_, _Rational] /; brDependsQ[b, vars]];

brMonoSpec::usage = "brMonoSpec[monos] is the monomial list in the form Risch`LogReducible / Risch`ExpReducible take: {{t, kind, Dt}, ...}.";
brMonoSpec[monos_List] := {#["sym"], #["kind"], #["D"]} & /@ monos;

brRules::usage = "brRules[x, monos] is the derivation {x -> 1, t_i -> D t_i, ...} of the monomial list.";
brRules[x_Symbol, monos_List] := Prepend[(#["sym"] -> #["D"]) & /@ monos, x -> 1];

brNewSymbol::usage = "brNewSymbol[i] is the i-th tower symbol brT<i> (in Bronstein`Private`), a deterministic name so tests can name the monomials.";
brNewSymbol[i_Integer] := Symbol["Bronstein`Private`brT" <> ToString[i]];

brAddMonomial::usage = "brAddMonomial[x, monos, kind, arg, kernel] appends the monomial t = Log[arg] (kind \"Log\", D t = D arg/arg) or E^arg (kind \"Exp\", D t = t D arg) to monos, arg in the symbols of monos. It returns {monos', t}.";
brAddMonomial[x_Symbol, monos_List, kind_String, arg_, kernel_] := Module[{t, Da, Dt},
  t = brNewSymbol[Length[monos] + 1];
  Da = Risch`Derivation[arg, brRules[x, monos]];
  Dt = Together[If[kind === "Log", Da/arg, t Da]];
  If[Dt === 0, throwInternal["brAddMonomialConstant", <|"kernel" -> kernel|>]];
  {Append[monos, <|"sym" -> t, "kind" -> kind, "arg" -> arg, "D" -> Dt, "kernel" -> kernel|>], t}];

brRebaseExp::usage = "brRebaseExp[x, monos, i, q] replaces the exponential monomial t_i = E^b by u^q with u = E^(b/q) a fresh monomial in the same position (D u = u D b/q), substituting t_i -> u^q in every later monomial. It returns {monos', t_i -> u^q}. The identity E^b = (E^(b/q))^q holds everywhere.";
brRebaseExp[x_Symbol, monos_List, i_Integer, q_Integer] := Module[{m = monos, old, u, b, rule, Du},
  old = m[[i]]["sym"]; b = m[[i]]["arg"];
  u = Symbol["Bronstein`Private`brT" <> ToString[i] <> "r" <> ToString[q]];
  rule = old -> u^q;
  Du = Together[u Risch`Derivation[b, brRules[x, Take[m, i - 1]]]/q];
  m[[i]] = <|"sym" -> u, "kind" -> "Exp", "arg" -> Together[b/q], "D" -> Du,
             "kernel" -> E^(Together[(m[[i]]["kernel"] /. Power[E, e_] :> e)/q])|>;
  Do[m[[j]] = Join[m[[j]], <|"arg" -> Together[m[[j]]["arg"] /. rule], "D" -> Together[m[[j]]["D"] /. rule]|>],
    {j, i + 1, Length[m]}];
  {m, rule}];

brExpRewrite::usage = "brExpRewrite[x, monos, b] decides E^b over the tower monos. It returns {\"new\"} when E^b is a new monomial (Risch`ExpReducible: False), {\"rewrite\", e} with E^b = e an identity in the symbols, or {\"rebase\", i, q} when the relation has a fractional exponent on the exponential monomial i (re-base it as (E^(b_i/q))^q first). A fractional power of a logarithm's argument is a radical: Decline.";
brExpRewrite[x_Symbol, monos_List, b_] := Module[{r, c, e, frac, i0},
  (* No Return inside With/If here: a nested Return only leaves the inner scope. *)
  r = If[monos === {}, False, Risch`ExpReducible[b, x, brMonoSpec[monos]]];
  If[r === False, {"new"},
    If[!ListQ[r] || Length[r] =!= Length[monos], throwDecline["exponential relation not decided"]];
    frac = Select[Range[Length[monos]], !IntegerQ[r[[#]]] &];
    i0 = If[frac === {}, None, First[frac]];
    Which[
      i0 =!= None && monos[[i0]]["kind"] === "Exp", {"rebase", i0, Denominator[r[[i0]]]},
      i0 =!= None, throwDecline["exponential of a fractional multiple of a logarithm (a radical)"],
      True,
        c = Together[b - Sum[r[[i]] If[monos[[i]]["kind"] === "Exp", monos[[i]]["arg"], monos[[i]]["sym"]],
                            {i, Length[monos]}]];
        If[brDependsQ[c, Prepend[#["sym"] & /@ monos, x]],
          throwInternal["brExpRewriteNonConstant", <|"b" -> b, "c" -> c|>]];
        e = E^c Product[If[monos[[i]]["kind"] === "Exp", monos[[i]]["sym"]^r[[i]], monos[[i]]["arg"]^r[[i]]],
                        {i, Length[monos]}];
        {"rewrite", e}]]];

brKernelsToSymbols::usage = "brKernelsToSymbols[f, x] replaces the Log/Exp kernels of f innermost first by tower monomials. It returns {fs, monos}: fs the integrand in x and the symbols, monos the monomial list (Associations with \"sym\", \"kind\", \"arg\", \"D\", \"kernel\"). Throws a Decline for a reducible logarithm, a radical inside a kernel, or a fractional power of a logarithm's argument.";
brKernelsToSymbols[f0_, x_Symbol] := Module[{f = f0, monos = {}, vars, ks, k, arg, kind, step, t, res, rule, guard = 0},
  While[True,
    vars = Prepend[#["sym"] & /@ monos, x];
    ks = brInnermostKernels[f, vars];
    If[ks === {}, Break[]];
    If[++guard > 200, throwDecline["tower too large"]];
    k = First[ks];
    {kind, arg} = Switch[k, Log[_], {"Log", First[k]}, _, {"Exp", Last[k]}];
    (* A radical inside a kernel: Shape B's curve-monomial case (BaseCurve.m)
       takes exactly one such kernel, its argument in Q(x)(y) with y^n = b(x).
       The reducibility tests (Risch`LogReducible / ExpReducible) do not know
       the relation y^n = b, so a radical argument is only admitted as the
       FIRST and ONLY monomial, and no monomial may follow it. *)
    If[AnyTrue[monos, TrueQ[#["algebraic"]] &],
      throwDecline["a further monomial over a kernel with a radical argument (Shape B, curve-monomial case, takes one kernel)"]];
    If[brRadicalQ[arg, vars],
      If[monos =!= {} || !brRadicalOverXQ[arg, x],
        throwDecline["radical inside a logarithm or exponential"]];
      {monos, t} = brAddMonomial[x, monos, kind, Together[arg], k];
      monos[[-1]] = Append[monos[[-1]], "algebraic" -> True];
      f = f /. k -> t;
      Continue[]];
    If[kind === "Log",
      If[monos =!= {} && Risch`LogReducible[arg, x, brMonoSpec[monos]] =!= False,
        throwDecline["reducible logarithm (not an identity on every branch)", <|"kernel" -> k|>]];
      {monos, t} = brAddMonomial[x, monos, "Log", Together[arg], k /. brBack[monos]];
      f = f /. k -> t,
      (* Exp *)
      step = brExpRewrite[x, monos, arg];
      Switch[First[step],
        "new",
          {monos, t} = brAddMonomial[x, monos, "Exp", Together[arg], k /. brBack[monos]];
          f = f /. k -> t,
        "rewrite", f = f /. k -> step[[2]],
        "rebase",
          {monos, rule} = brRebaseExp[x, monos, step[[2]], step[[3]]];
          f = f /. rule]]];
  {f, monos}];

brBack::usage = "brBack[monos] is the list of back-substitution rules t_i -> surface kernel (each surface kernel is written in x alone, so one ReplaceAll pass suffices).";
brBack[monos_List] := (#["sym"] -> #["kernel"]) & /@ monos;

brToSurface::usage = "brToSurface[T, e] is the expression e in the tower symbols of T written back in the surface kernels.";
brToSurface[T_Association, e_] := e /. T["back"];

brDetectRadical::usage = "brDetectRadical[fs, vars] finds the single radical of fs (in the symbols vars). It returns <|\"b\" -> b, \"n\" -> n|> with every radical subexpression a power of b^(1/n), or None; throws a Decline on nested radicals, several radicands (a compositum; a Kummer split is a later phase), or a radicand that is not rational in vars.";
brDetectRadical[fs_, vars_List] := Module[{pows, bases, b0, exps},
  pows = DeleteDuplicates[Cases[fs, Power[b_, r_Rational] /; brDependsQ[b, vars], {0, Infinity}]];
  If[pows === {}, Return[None]];
  If[!FreeQ[First /@ pows, Power[_, _Rational]], throwDecline["nested radicals"]];
  If[!FreeQ[fs, Power[c_, _Rational] /; !brDependsQ[c, vars]],
    throwDecline["algebraic constant in the integrand"]];
  bases = Together[First[#]] & /@ pows;
  b0 = First[bases];
  If[!brRationalInQ[b0, vars], throwDecline["radicand is not rational in the tower"]];
  exps = MapThread[
    Which[
      Together[#1 - b0] === 0, #2,
      Together[#1 b0 - 1] === 0, -#2,
      True, throwDecline["more than one radicand (a radical compositum)"]] &,
    {bases, Last /@ pows}];
  <|"b" -> b0, "n" -> LCM @@ (Denominator /@ exps)|>];

brRadicalOverXQ::usage = "brRadicalOverXQ[e, x] is True when every radical b^r (r non-integer) in e has a radicand b rational in x alone (no kernel, no other symbol, no nested radical): e is then in Q(x)(y) for the radicals of e.";
brRadicalOverXQ[e_, x_Symbol] := AllTrue[
  Cases[e, Power[b_, _Rational] /; !FreeQ[b, x], {0, Infinity}], brRationalInQ[First[#], {x}] &];

brRationalInQ::usage = "brRationalInQ[e, vars] is True when e is a rational function of the symbols vars with rational coefficients (no other symbol, no radical, no kernel).";
brRationalInQ[e_, vars_List] := Module[{t = Together[e], num, den},
  num = Expand[Numerator[t]]; den = Expand[Denominator[t]];
  PolynomialQ[num, vars] && PolynomialQ[den, vars] &&
    Complement[Variables[{num, den}], vars] === {} &&
    FreeQ[{num, den}, Power[_, _Rational] | _Complex | _Real]];

brTower::usage = "brTower[f, x] is the tower record of the integrand f (see the header of Tower.m): kernels replaced by monomials, the radical detected, the shape classified and, for Shape A, theta moved to the top. Throws a Decline outside the scope.";
brTower[f0_, x_Symbol] := Module[{f, fs, monos, vars, rad, shape, theta = None, kind = None, lower, bsyms, pos, above},
  f = f0 /. Exp[a_] :> E^a;
  brCheckHeads[f, x];
  {fs, monos} = brKernelsToSymbols[f, x];
  vars = Prepend[#["sym"] & /@ monos, x];
  If[!FreeQ[fs, brKernelPattern[]] && brDependsQ[Cases[fs, brKernelPattern[], {0, Infinity}], vars],
    throwInternal["brTowerKernelLeft", <|"f" -> fs|>]];
  rad = brDetectRadical[{fs, Select[monos, TrueQ[#["algebraic"]] &] /. a_Association :> a["arg"]}, vars];
  If[!brRationalInQ[fs /. Power[b_, r_Rational] /; brDependsQ[b, vars] :> brRadicalPlaceholder, Append[vars, brRadicalPlaceholder]],
    throwDecline["integrand is not algebraic over the tower"]];
  shape = Which[
    rad === None && monos === {}, "Rational",
    rad === None, "Transcendental",
    monos === {}, "Algebraic",
    AnyTrue[monos, TrueQ[#["algebraic"]] &],
      (* brKernelsToSymbols admitted the radical-argument kernel as the only
         monomial, its radicand rational in x *)
      If[Length[monos] =!= 1 || brDependsQ[rad["b"], #["sym"] & /@ monos],
        throwInternal["brTowerCurveMonomialShape", <|"f" -> fs|>]];
      "BCurve",
    True,
      bsyms = Select[#["sym"] & /@ monos, !FreeQ[rad["b"], #] &];
      If[bsyms === {}, "B", "A"]];
  lower = #["sym"] & /@ monos;
  If[shape === "A",
    theta = Last[bsyms];
    pos = First[FirstPosition[#["sym"] & /@ monos, theta]];
    above = Drop[monos, pos];
    (* a monomial above theta that depends on it: the radical sits in the   *)
    (* middle of the tower (theta, y, then something over both): out of scope *)
    If[AnyTrue[above, !FreeQ[{#["D"], #["arg"]}, theta] &],
      throwDecline["a monomial above the radicand's top monomial depends on it"]];
    monos = Join[Delete[monos, pos], {monos[[pos]]}];
    kind = Last[monos]["kind"];
    lower = Most[#["sym"] & /@ monos]];
  <|"x" -> x, "monos" -> monos, "syms" -> (#["sym"] & /@ monos), "kinds" -> (#["kind"] & /@ monos),
    "args" -> (#["arg"] & /@ monos), "rules" -> brRules[x, monos],
    "back" -> brBack[monos], "kernels" -> (#["kernel"] & /@ monos),
    "f" -> fs, "radical" -> rad, "shape" -> shape, "theta" -> theta, "kind" -> kind,
    "lower" -> lower|>];

brRadicalPlaceholder::usage = "brRadicalPlaceholder is the inert symbol standing for the radical while brTower checks that the integrand is rational in the tower and the radical.";

brFromSurface::usage = "brFromSurface[T, e] writes the surface expression e (in x and Log/Exp kernels, e.g. an answer of the C integrator) in the tower symbols of T: every kernel of T is replaced by its symbol and every other exponential reducible over the tower is rewritten (brExpRewrite). Kernels that are not in the tower (a new logarithm) are left in place.";
brFromSurface[T_Association, e_] := Module[{r = e /. Exp[a_] :> E^a, monos = T["monos"], ks},
  r = r /. ((#["kernel"] -> #["sym"]) & /@ monos);
  r = r /. Power[E, b_] /; brDependsQ[b, Prepend[T["syms"], T["x"]]] :>
    With[{st = brExpRewrite[T["x"], monos, b /. ((#["kernel"] -> #["sym"]) & /@ monos)]},
      If[First[st] === "rewrite", st[[2]], Power[E, b]]];
  r];

brLowerTower::usage = "brLowerTower[T] is the tower record of k alone (T without theta): the monomials below theta, with their rules, symbols and kernels; shape \"Lower\".";
brLowerTower[T_Association] := Module[{monos = If[T["theta"] === None, T["monos"], Most[T["monos"]]]},
  <|"x" -> T["x"], "monos" -> monos, "syms" -> (#["sym"] & /@ monos), "kinds" -> (#["kind"] & /@ monos),
    "args" -> (#["arg"] & /@ monos), "rules" -> brRules[T["x"], monos], "back" -> brBack[monos],
    "kernels" -> (#["kernel"] & /@ monos), "f" -> None, "radical" -> None, "shape" -> "Lower",
    "theta" -> None, "kind" -> None, "lower" -> (#["sym"] & /@ monos)|>];

brTowerD::usage = "brTowerD[T, e] is the tower derivation of e (in x and the symbols of T): Risch`Derivation[e, T[\"rules\"]].";
brTowerD[T_Association, e_] := Risch`Derivation[e, T["rules"]];

brEta::usage = "brEta[T] is, for Shape A, the derivative data of theta: Dtheta for a logarithm (an element of k), eta' = Dtheta/theta for an exponential (an element of k).";
brEta[T_Association] := With[{d = Last[T["rules"]][[2]], th = T["theta"]},
  If[T["kind"] === "Exp", Together[d/th], d]];

(* -------------------------------------------------------------------------- *)
(* Linear (Moebius) radicand in x: a change of variable that removes the     *)
(* radical entirely.                                                          *)
(* -------------------------------------------------------------------------- *)

brMoebiusRadicandQ::usage = "brMoebiusRadicandQ[b, x] is True when b = (alpha x + beta)/(gamma x + delta) with rational alpha, beta, gamma, delta and alpha delta - beta gamma != 0 (a linear radicand is the case gamma = 0).";
brMoebiusRadicandQ[b_, x_Symbol] := Module[{t = Together[b], num, den},
  num = Numerator[t]; den = Denominator[t];
  PolynomialQ[num, x] && PolynomialQ[den, x] && Exponent[num, x] <= 1 && Exponent[den, x] <= 1 &&
    Complement[Variables[{num, den}], {x}] === {} &&
    (Exponent[num, x] === 1 || Exponent[den, x] === 1) &&
    Expand[Coefficient[num, x, 1] Coefficient[den, x, 0] - Coefficient[num, x, 0] Coefficient[den, x, 1]] =!= 0];

brFlattenLinear::usage = "brFlattenLinear[f, x, b, n, u] changes variable u = b^(1/n) for a Moebius radicand b(x) (brMoebiusRadicandQ): x = (delta u^n - beta)/(alpha - gamma u^n), every power b^r becomes u^(r n), and the result is g(u) = f(x(u)) x'(u), a function of u with no radical in u. It returns <|\"g\" -> g, \"u\" -> u, \"x\" -> x(u), \"back\" -> (u -> b^(1/n))|>, or None when b is not Moebius. Integrating g du and substituting back is an antiderivative of f dx wherever u = b^(1/n) is the principal branch (the final verification decides).";
brFlattenLinear[f_, x_Symbol, b_, n_Integer, u_Symbol] := Module[{t = Together[b], num, den, al, be, ga, de, xu, g},
  If[!brMoebiusRadicandQ[b, x], Return[None]];
  num = Numerator[t]; den = Denominator[t];
  al = Coefficient[num, x, 1]; be = Coefficient[num, x, 0];
  ga = Coefficient[den, x, 1]; de = Coefficient[den, x, 0];
  xu = Together[(de u^n - be)/(al - ga u^n)];
  g = f /. Power[c_, r_Rational] /; !FreeQ[c, x] :>
        Which[Together[c - b] === 0, u^(r n), Together[c b - 1] === 0, u^(-r n), True, Power[c, r]];
  If[brRadicalQ[g, {x}], Return[None]];
  g = (g /. x -> xu) D[xu, u];
  <|"g" -> g, "u" -> u, "x" -> xu, "back" -> (u -> b^(1/n))|>];

(* -------------------------------------------------------------------------- *)
(* Quadratic radicand: the conic y^2 = Q(x) is rational when it has a rational *)
(* point, and the parametrisation removes the radical (Euler's substitutions). *)
(* -------------------------------------------------------------------------- *)

brConicPoint::usage = "brConicPoint[Q, x] is a rational point of the conic y^2 = Q(x) (Q a squarefree quadratic over Q) usable by brEulerFlatten: {\"Infinity\", alpha} when the leading coefficient is alpha^2 (the points at infinity are rational), else {x0, y0} with Q(x0) = y0^2 for the first such x0 in a small search (x0 = 0, +-1, ..., +-12, +-1/2, ..., +-12/7), else None. y0 may be 0 (a rational root). Not finding a point is not a proof that none exists: the caller then takes the general method.";
brConicPoint[Q_, x_Symbol] := Module[{a = Coefficient[Q, x, 2], v, cands, sq},
  sq[r_] := With[{w = Sqrt[r]}, If[r >= 0 && Element[w, Rationals] === True, w, None]];
  If[sq[a] =!= None, Return[{"Infinity", sq[a]}]];
  cands = DeleteDuplicates[Flatten[Table[{p/q, -p/q}, {q, 1, 7}, {p, 0, 12}]]];
  Do[v = Expand[Q /. x -> c];
    If[sq[v] =!= None, Return[{c, sq[v]}, Module]], {c, cands}];
  None];

brEulerFlatten::usage = "brEulerFlatten[f, x, b, u] removes the square root of a squarefree quadratic radicand b(x) over Q by Euler's substitution through a rational point of y^2 = b (brConicPoint): with the point at infinity (b = alpha^2 x^2 + ...) u = y + alpha x; with a point (x0, y0), u = (y - y0)/(x - x0). Then x and y = b^(1/2) are rational in u, every power b^r (r with denominator 2) becomes y(u)^(2r), and the result is g(u) = f(x(u)) x'(u). It returns <|\"g\", \"u\", \"x\" -> x(u), \"back\" -> (u -> u(x, Sqrt[b]))|> like brFlattenLinear, or None (b not such a quadratic, no rational point found, or a radical of x left). Q(x)(y) = Q(u), so integrating g du decides f dx, and u is defined from the user's own Sqrt[b], so x(u(x)) = x and y(u(x)) = Sqrt[b] hold identically.";
brEulerFlatten[f_, x_Symbol, b_, u_Symbol] := Module[{Q = Expand[b], pt, al, x0, y0, d1, xu, yu, back, g},
  If[!PolynomialQ[Q, x] || Exponent[Q, x] =!= 2 || Complement[Variables[Q], {x}] =!= {} ||
     Exponent[PolynomialGCD[Q, D[Q, x]], x] > 0, Return[None]];
  pt = brConicPoint[Q, x];
  If[pt === None, Return[None]];
  If[First[pt] === "Infinity",
    al = pt[[2]];
    (* y = u - al x:  u^2 - 2 al x u = b1 x + b0 *)
    xu = Together[(u^2 - Coefficient[Q, x, 0])/(2 al u + Coefficient[Q, x, 1])];
    yu = Together[u - al xu];
    back = u -> Sqrt[b] + al x,
    {x0, y0} = pt;
    (* Q = y0^2 + Q'(x0) t + a t^2, t = x - x0;  y = y0 + u t *)
    d1 = D[Q, x] /. x -> x0;
    xu = Together[x0 + (d1 - 2 y0 u)/(u^2 - Coefficient[Q, x, 2])];
    yu = Together[y0 + u (xu - x0)];
    back = u -> (Sqrt[b] - y0)/(x - x0)];
  g = f /. Power[c_, r_Rational] /; !FreeQ[c, x] && Denominator[r] === 2 :>
        Which[Together[c - b] === 0, yu^(2 r), Together[c b - 1] === 0, yu^(-2 r), True, Power[c, r]];
  If[brRadicalQ[g, {x}], Return[None]];
  g = Together[(g /. x -> xu) D[xu, u]];
  <|"g" -> g, "u" -> u, "x" -> xu, "back" -> back|>];

brGenus0Flatten::usage = "brGenus0Flatten[f, x, u] is the genus-0 pre-pass of Shape B on the SURFACE integrand f: when f has Log/Exp kernels and exactly one radicand b, rational in x alone (radicals inside kernels included), and the curve y^n = b is rational by brFlattenLinear (Moebius b) or brEulerFlatten (squarefree quadratic b with a rational point, n = 2), it returns that change-of-variable record (with \"b\" -> b and \"n\" -> n added); otherwise None. Q(x)(y) = Q(u), so the transcendental Risch algorithm in u then decides the integral, whatever the kernels contain.";
brGenus0Flatten[f_, x_Symbol, u_Symbol] := Module[{pows, bases, b0, exps, n, r},
  If[FreeQ[f, (Log[a_] | Power[E, a_]) /; !FreeQ[a, x]], Return[None]];
  pows = DeleteDuplicates[Cases[f, Power[c_, _Rational] /; !FreeQ[c, x], {0, Infinity}]];
  If[pows === {}, Return[None]];
  bases = Together[First[#]] & /@ pows;
  b0 = First[bases];
  If[!brRationalInQ[b0, {x}], Return[None]];
  exps = MapThread[Which[Together[#1 - b0] === 0, #2, Together[#1 b0 - 1] === 0, -#2, True, None] &,
    {bases, Last /@ pows}];
  If[MemberQ[exps, None], Return[None]];
  n = LCM @@ (Denominator /@ exps);
  r = brFlattenLinear[f, x, b0, n, u];
  If[r === None && n === 2, r = brEulerFlatten[f, x, b0, u]];
  If[r === None, None, Join[r, <|"b" -> b0, "n" -> n|>]]];

brCurveCanon::usage = "brCurveCanon[e, x, b, n] writes e, an expression in x, the radical y = b^(1/n) (b rational in x) and Log/Exp/ArcTan/ArcTanh kernels, in the canonical form of Q(x)(y): every kernel argument first (innermost first), then the whole expression with the kernels as opaque symbols, are reduced to Sum_i c_i y^i (i < n), c_i rational, the denominator rationalised in Q(x)[y]/(y^n - b). Since y = b^(1/n) satisfies y^n = b pointwise this is an identity of functions; it is the presentation of an answer found in a rational parametrisation (brGenus0Flatten), whose logarithms otherwise carry the parametrisation's quotients (Log[(x^2 + x Sqrt[1 + x^2])/(x + Sqrt[1 + x^2])] for Log[x]).";
brCurveCanon[e_, x_Symbol, b_, n_Integer] := Module[{Y = Unique["Bronstein`Private`brCY"], toY, canon, ks, atoms, r},
  toY[ex_] := ex /. Power[c_, r0_Rational] /; !FreeQ[c, x] :>
    Which[Together[c - b] === 0, Y^(r0 n), Together[c b - 1] === 0, Y^(-r0 n), True, Power[c, r0]];
  canon[ex0_] := Module[{t = Together[ex0], num, den, inv, out},
    num = Numerator[t]; den = Denominator[t];
    If[!FreeQ[den, Y],
      inv = powerBasisInverse[reduceYDegree[den, Y, n, b], Y, n, b, x];
      num = Expand[num (inv . (Y^Range[0, n - 1]))]; den = 1];
    out = PadRight[CoefficientList[reduceYDegree[num, Y, n, b], Y], n];
    Sum[Together[out[[i + 1]]/den] Y^i, {i, 0, n - 1}]];
  r = toY[e];
  (* kernel arguments, innermost first: a kernel's argument has no kernel left *)
  r = r //. (h : (Log | ArcTan | ArcTanh))[a_] /; FreeQ[a, Log | ArcTan | ArcTanh | Power[E, _] | brCanonDone] :>
        brCanonDone[h, canon[a]];
  r = r //. Power[E, a_] /; FreeQ[a, Log | ArcTan | ArcTanh | Power[E, _] | brCanonDone] && !FreeQ[a, x] :>
        brCanonDone[Exp, canon[a]];
  ks = DeleteDuplicates[Cases[r, _brCanonDone, {0, Infinity}]];
  atoms = Table[Unique["Bronstein`Private`brCK"], {Length[ks]}];
  r = canon[r /. Thread[ks -> atoms]] /. Thread[atoms -> ks];
  r /. brCanonDone[h_, a_] :> h[a] /. Y -> b^(1/n)];

brCanonDone::usage = "brCanonDone[h, a] is the inert head brCurveCanon gives a kernel h[a] whose argument a is already canonical.";
