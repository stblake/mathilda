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

brUnsupportedHeads::usage = "brUnsupportedHeads[] is the pattern of every head the Bronstein front end does not take: brElementaryUnimplementedHeads[] (elementary, not yet supported) | brAlgebraicGeneralHeads[] (a general algebraic function) | brNonElementaryHeads[] (outside elementary functions).";
brUnsupportedHeads[] := brElementaryUnimplementedHeads[] | brAlgebraicGeneralHeads[] | brNonElementaryInputHeads[];

brElementaryUnimplementedHeads::usage = "brElementaryUnimplementedHeads[] are the elementary heads the front end does not yet rewrite (trigonometric, hyperbolic and their inverses; the TrigToExp front end over Q(i) of the completion campaign). An integrand with one of them is inside the decided class: implementation incomplete (Unimplemented).";
brElementaryUnimplementedHeads[] := Sin | Cos | Tan | Cot | Sec | Csc | ArcSin | ArcCos | ArcTan | ArcCot |
  ArcSec | ArcCsc | Sinh | Cosh | Tanh | Coth | Sech | Csch | ArcSinh | ArcCosh | ArcTanh |
  ArcCoth | ArcSech | ArcCsch;

brAlgebraicGeneralHeads::usage = "brAlgebraicGeneralHeads[] are the heads of a general algebraic function (Root, RootSum): outside one simple radical extension, so implementation incomplete (Incomplete).";
brAlgebraicGeneralHeads[] := Root | RootSum;

brNonElementaryInputHeads::usage = "brNonElementaryInputHeads[] are heads outside the elementary functions (special functions, Abs/Sign, piecewise and integer-valued functions, unevaluated derivatives and integrals): no elementary decision applies, a silent Decline.";
brNonElementaryInputHeads[] := Erf | Erfc | Erfi | ExpIntegralEi | LogIntegral | PolyLog |
  Gamma | Abs | Sign | Derivative | Integrate | Piecewise | Floor | Ceiling;

brCheckHeads::usage = "brCheckHeads[f, x] classifies the heads and constants of f: a non-elementary head or an approximate number is a Decline (no elementary decision applies), a Root/RootSum is Incomplete (a general algebraic function), an elementary head not yet rewritten (trigonometric, hyperbolic, inverse) or a complex constant is Unimplemented (inside the decided class). Otherwise returns True.";
brCheckHeads[f_, x_Symbol] := (
  If[!FreeQ[f, brNonElementaryInputHeads[]], throwDecline["non-elementary function in the integrand"]];
  If[!FreeQ[f, _Real], throwDecline["approximate number in the integrand"]];
  If[!FreeQ[f, brAlgebraicGeneralHeads[]], throwIncomplete["a general algebraic function (Root/RootSum) in the integrand"]];
  If[!FreeQ[f, brElementaryUnimplementedHeads[]],
    throwUnimplemented["Tower:trig", "trigonometric or hyperbolic function (TrigToExp front end over Q(i))"]];
  (* Gaussian-rational constants: the constant field is Q(i), carried with a
     literal I (Gaussian arithmetic reduces I^2 = -1 everywhere; residue
     fields express I in their primitive element). Anything else complex
     (an inexact or symbolic part) stays unimplemented. *)
  If[!FreeQ[f, c_Complex /; !brGaussianRationalQ[c]],
    throwUnimplemented["Tower:complex", "complex constant that is not a Gaussian rational"]];
  True);

brGaussianRationalQ::usage = "brGaussianRationalQ[c] is True when c is a Gaussian rational a + b I (a, b rational).";
brGaussianRationalQ[c_] := MatchQ[c, Complex[(_Integer | _Rational), (_Integer | _Rational)]];

brConstantField::usage = "brConstantField[e] is the constant field the expression e needs: \"Q\", or \"QI\" when a Gaussian rational occurs (Phase 4.1; bronstein_extensions.tex E4.1).";
brConstantField[e_] := If[FreeQ[e, _Complex], "Q", "QI"];

brAlgConstantsOf::usage = "brAlgConstantsOf[e] is the list of radical constants of e (tragerConstantsOf): the generators of its number-field constant field besides I (Phase 4.3).";
brAlgConstantsOf[e_] := tragerConstantsOf[e];

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

brExpRewrite::usage = "brExpRewrite[x, monos, b] decides E^b over the tower monos. It returns {\"new\"} when E^b is a new monomial (Risch`ExpReducible: False), {\"rewrite\", e} with E^b = e an identity in the symbols, or {\"rebase\", i, q} when the relation has a fractional exponent on the exponential monomial i (re-base it as (E^(b_i/q))^q first). A fractional multiple r of a logarithm Log[a] becomes the radical a^r (E^(r Log[a]) = a^r pointwise).";
brExpRewrite[x_Symbol, monos_List, b_] := Module[{r, c, e, frac, i0},
  (* No Return inside With/If here: a nested Return only leaves the inner scope. *)
  r = If[monos === {}, False, Risch`ExpReducible[b, x, brMonoSpec[monos]]];
  If[r === False, {"new"},
    If[!ListQ[r] || Length[r] =!= Length[monos], throwUnimplemented["Tower:exprel", "exponential relation not decided"]];
    (* a fractional exponent on an exponential monomial re-bases it; on a
       logarithm Log[a] it is a radical: E^(r Log[a]) = a^r holds pointwise for
       the principal branches (exp(r log a) is the definition of a^r) *)
    frac = Select[Range[Length[monos]], !IntegerQ[r[[#]]] && monos[[#]]["kind"] === "Exp" &];
    i0 = If[frac === {}, None, First[frac]];
    Which[
      i0 =!= None, {"rebase", i0, Denominator[r[[i0]]]},
      True,
        c = Together[b - Sum[r[[i]] If[monos[[i]]["kind"] === "Exp", monos[[i]]["arg"], monos[[i]]["sym"]],
                            {i, Length[monos]}]];
        If[brDependsQ[c, Prepend[#["sym"] & /@ monos, x]],
          throwInternal["brExpRewriteNonConstant", <|"b" -> b, "c" -> c|>]];
        e = E^c Product[If[monos[[i]]["kind"] === "Exp", monos[[i]]["sym"]^r[[i]], monos[[i]]["arg"]^r[[i]]],
                        {i, Length[monos]}];
        {"rewrite", e}]]];

brSplitLogExponent::usage = "brSplitLogExponent[x, monos, b] splits the exponent b (in the tower symbols) as b = u + Sum r_i t_i, t_i the logarithm monomials of monos (t_i = Log[a_i], not radical-argument ones) with r_i rational and b linear in each such t_i, u the rest. It returns {u, Prod a_i^r_i} when some r_i is nonzero and u depends on x or a tower symbol, so that E^b = E^u Prod a_i^r_i holds pointwise (exp(v + w) = exp(v) exp(w), exp(r log a) = a^r for the principal branches); otherwise None (no log part, or E^b is a pure product the structure theorem handles).";
brSplitLogExponent[x_Symbol, monos_List, b_] := Module[{be = Expand[b], logs, rs, u},
  logs = Select[monos, #["kind"] === "Log" && !TrueQ[#["algebraic"]] &];
  If[logs === {}, Return[None]];
  rs = Function[m, With[{t = m["sym"]},
      If[PolynomialQ[be, t] && Exponent[be, t] === 1 && MatchQ[Coefficient[be, t, 1], _Integer | _Rational],
        Coefficient[be, t, 1], 0]]] /@ logs;
  If[AllTrue[rs, # === 0 &], Return[None]];
  u = Expand[be - rs . (#["sym"] & /@ logs)];
  If[!brDependsQ[u, Prepend[#["sym"] & /@ monos, x]], Return[None]];
  {u, Times @@ MapThread[Power[#1["arg"], #2] &, {logs, rs}]}];

brKappaMode::usage = "brKappaMode is True inside brReducibleLogSplit: brKernelsToSymbols then replaces a reducible logarithm Log[a] = L0 + kappa by L0 + K (K a fresh constant symbol) and records {K, kappa} in brKappas, instead of throwing Tower:reducible-log.";
brKappaMode = False;
brKappas::usage = "brKappas is the list {{K, kappa}, ...} brKernelsToSymbols fills in brKappaMode: K the constant symbol, kappa = Log[a] - L0 its surface value (locally constant: D kappa = 0).";
brKappas = {};

brReducibleLogSplit::usage = "brReducibleLogSplit[f, x] splits the surface integrand f along its reducible logarithms (Tower:reducible-log). A logarithm Log[a] reducible over the tower, D Log[a] = D L0 with L0 = Sum r_i t_i + Sum s_j b_j (Risch`LogReducible), is Log[a] = L0 + kappa with kappa = Log[a] - L0 LOCALLY CONSTANT (a branch difference: 0 or a transcendental constant on each region). Writing kappa as a constant symbol K, f is a polynomial in the K's (else None: K in a denominator is a symbolic parameter, Phase 4). It returns <|\"kappas\" -> {{K, kappa}, ...}, \"terms\" -> {{K^alpha, g_alpha}, ...}|> with f = Sum K^alpha g_alpha, every g_alpha free of the reducible logarithms (surface form), or None when f has no reducible logarithm or a later monomial depends on one (K would enter the tower). Sum kappa^alpha Int g_alpha is an antiderivative of f (D kappa = 0); see brIntegrateBody for the certificates.";
brReducibleLogSplit[f_, x_Symbol] := Module[{r, fs, monos, Ks, poly, rules, terms},
  r = Block[{brKappaMode = True, brKappas = {}}, With[{kk = brCatch[brKernelsToSymbols[f, x]]}, {kk, brKappas}]];
  If[brFailQ[r[[1]]] || r[[2]] === {}, Return[None]];
  {fs, monos} = r[[1]];
  Ks = First /@ r[[2]];
  (* E6a needs k' free of Log[a]: nothing above the reducible logarithm may
     depend on it (Log[Log[x^2]], E^(Log[x^2]^2) would put K in the tower) *)
  If[AnyTrue[monos, !FreeQ[{#["arg"], #["D"], #["kernel"]}, Alternatives @@ Ks] &], Return[None]];
  poly = Together[fs];
  If[!FreeQ[Denominator[poly], Alternatives @@ Ks], Return[None]];
  poly = Expand[Numerator[poly]];
  If[!PolynomialQ[poly, Ks], Return[None]];
  rules = CoefficientRules[poly, Ks];
  terms = ({Times @@ (Ks^#[[1]]), Together[#[[2]]/Denominator[Together[fs]]] /. brBack[monos]} &) /@ rules;
  <|"kappas" -> r[[2]], "terms" -> terms|>];

brKernelsToSymbols::usage = "brKernelsToSymbols[f, x] replaces the Log/Exp kernels of f innermost first by tower monomials. It returns {fs, monos}: fs the integrand in x and the symbols, monos the monomial list (Associations with \"sym\", \"kind\", \"arg\", \"D\", \"kernel\"). Throws a Decline for a reducible logarithm, a radical inside a kernel, or a fractional power of a logarithm's argument.";
brKernelsToSymbols[f0_, x_Symbol] := Module[{f = f0, monos = {}, vars, ks, k, arg, kind, step, t, res, rule, guard = 0},
  While[True,
    vars = Prepend[#["sym"] & /@ monos, x];
    ks = brInnermostKernels[f, vars];
    If[ks === {}, Break[]];
    If[++guard > 200, throwDecline["tower too large"]];  (* a resource bound *)
    k = First[ks];
    {kind, arg} = Switch[k, Log[_], {"Log", First[k]}, _, {"Exp", Last[k]}];
    (* A radical inside a kernel: Shape B's curve-monomial case (BaseCurve.m)
       takes exactly one such kernel, its argument in Q(x)(y) with y^n = b(x).
       The reducibility tests (Risch`LogReducible / ExpReducible) do not know
       the relation y^n = b, so a radical argument is only admitted as the
       FIRST and ONLY monomial, and no monomial may follow it. *)
    If[AnyTrue[monos, TrueQ[#["algebraic"]] &],
      throwUnimplemented["Tower:curve-monomial-tower", "a further monomial over a kernel with a radical argument (transcendental Risch over an algebraic base)"]];
    If[brRadicalQ[arg, vars],
      If[monos =!= {} || !brRadicalOverXQ[arg, x],
        throwUnimplemented["Tower:radical-in-kernel", "radical inside a logarithm or exponential above other monomials (transcendental Risch over an algebraic base)"]];
      {monos, t} = brAddMonomial[x, monos, kind, Together[arg], k];
      monos[[-1]] = Append[monos[[-1]], "algebraic" -> True];
      f = f /. k -> t;
      Continue[]];
    If[kind === "Log",
      (* LogReducible is tri-state: False is a proof of a new monomial, a list
         a relation, $Failed an undecided relation (never read as "new") *)
      With[{lr = If[monos === {}, False, Risch`LogReducible[arg, x, brMonoSpec[monos]]]},
        If[lr =!= False && !ListQ[lr],
          throwUnimplemented["Tower:logrel", "logarithm relation not decided", <|"kernel" -> k|>]];
        If[ListQ[lr],
          (* Log[a] = L0 + kappa, L0 = Sum r_i (t_i | b_i), kappa locally
             constant: in the kappa mode of brReducibleLogSplit the kernel is
             replaced by L0 + K, K a fresh constant symbol standing for kappa *)
          If[!TrueQ[brKappaMode],
            throwUnimplemented["Tower:reducible-log", "reducible logarithm (rewrite with an explicit constant)", <|"kernel" -> k|>]];
          With[{L0 = Sum[lr[[i]] If[monos[[i]]["kind"] === "Exp", monos[[i]]["arg"], monos[[i]]["sym"]], {i, Length[monos]}],
                K = Symbol["Bronstein`Private`brKappa" <> ToString[Length[brKappas] + 1]]},
            AppendTo[brKappas, {K, (k /. brBack[monos]) - (L0 /. brBack[monos])}];
            f = f /. k -> L0 + K];
          Continue[]]];
      {monos, t} = brAddMonomial[x, monos, "Log", Together[arg], k /. brBack[monos]];
      f = f /. k -> t,
      (* Exp: first split off the part of the exponent that is a rational
         combination of logarithm monomials, E^(u + r Log[a]) = E^u a^r
         (pointwise), so E^(x + Log[x]/2) is Sqrt[x] E^x, not a new monomial *)
      With[{sp = brSplitLogExponent[x, monos, arg]},
        If[sp =!= None, f = f /. k -> E^sp[[1]] sp[[2]]; Continue[]]];
      step = brExpRewrite[x, monos, arg];
      Switch[First[step],
        "new",
          {monos, t} = brAddMonomial[x, monos, "Exp", Together[arg], k /. brBack[monos]];
          f = f /. k -> t,
        "rewrite", f = f /. k -> step[[2]],
        "rebase",
          {monos, rule} = brRebaseExp[x, monos, step[[2]], step[[3]]];
          f = f /. rule]]];
  {f, brPruneMonomials[f, monos]}];

brPruneMonomials::usage = "brPruneMonomials[f, monos] drops the monomials that neither f nor a kept monomial (through its argument or derivative) involves: a kernel rewritten away after it entered the tower (Log[Log[x]] once E^(x + Log[Log[x]]/2) became Sqrt[Log[x]] E^x). The order of the rest is kept.";
brPruneMonomials[f_, monos_List] := Module[{need, grow = True, syms = #["sym"] & /@ monos},
  need = Select[syms, !FreeQ[f, #] &];
  While[grow,
    With[{more = Select[syms, Function[s, !MemberQ[need, s] &&
        AnyTrue[Select[monos, MemberQ[need, #["sym"]] &], !FreeQ[{#["arg"], #["D"]}, s] &]]]},
      grow = more =!= {}; need = Join[need, more]]];
  Select[monos, MemberQ[need, #["sym"]] &]];

brBack::usage = "brBack[monos] is the list of back-substitution rules t_i -> surface kernel (each surface kernel is written in x alone, so one ReplaceAll pass suffices).";
brBack[monos_List] := (#["sym"] -> #["kernel"]) & /@ monos;

brToSurface::usage = "brToSurface[T, e] is the expression e in the tower symbols of T written back in the surface kernels.";
brToSurface[T_Association, e_] := e /. T["back"];

brDetectRadical::usage = "brDetectRadical[fs, vars] finds the single radical of fs (in the symbols vars). It returns <|\"b\" -> b, \"n\" -> n|> with every radical subexpression a power of b^(1/n), or None. Nested radicals, and several radicands whose Kummer group is not cyclic (brKummerCyclicQ), are outside one simple radical extension: Incomplete. A cyclic compositum (one simple radical in disguise), an algebraic constant, and a radicand that is not rational in vars are Unimplemented.";
brDetectRadical[fs_, vars_List] := Module[{pows, bases, b0, exps, cyc},
  pows = DeleteDuplicates[Cases[fs, Power[b_, r_Rational] /; brDependsQ[b, vars], {0, Infinity}]];
  If[pows === {}, Return[None]];
  (* a radical of a constant inside a radicand is a number-field constant, not
     nesting (Phase 4.3) *)
  If[!FreeQ[First /@ pows, Power[b_, _Rational] /; brDependsQ[b, vars]], throwIncomplete["nested radicals"]];
  If[!FreeQ[fs, Power[c_, _Rational] /; !brDependsQ[c, vars] && !(tragerAlgConstants =!= {} && rationalNumberQ[c])],
    throwUnimplemented["Tower:algebraic-constant", "algebraic constant in the integrand (a number-field constant field)"]];
  bases = Together[First[#]] & /@ pows;
  b0 = First[bases];
  If[!brRationalInQ[b0, vars],
    throwUnimplemented["Tower:radicand-not-rational", "radicand is not rational in the tower (symbolic parameter)"]];
  exps = MapThread[
    Which[
      Together[#1 - b0] === 0, #2,
      Together[#1 b0 - 1] === 0, -#2,
      True, None] &,
    {bases, Last /@ pows}];
  If[MemberQ[exps, None],
    cyc = brKummerCyclicQ[Transpose[{bases, Last /@ pows}], vars];
    Which[
      cyc === False, throwIncomplete["several radicals generating a non-cyclic Kummer extension (not one simple radical)"],
      (* brCompositumRewrite (run first by the pipeline) found no pointwise
         single-radical presentation: the radicals agree with one radical only
         up to a locally constant root of unity, or over a larger constant field *)
      cyc === True, throwUnimplemented["Tower:cyclic-compositum", "several radicands generating one simple radical extension up to a locally constant unit or a constant field extension"],
      True, throwUnimplemented["Tower:compositum-undecided", "radical compositum whose Kummer group was not decided"]]];
  <|"b" -> b0, "n" -> LCM @@ (Denominator /@ exps)|>];

brKummerCyclicQ::usage = "brKummerCyclicQ[{{b1, r1}, ...}, vars] decides whether the radicals b_i^r_i (b_i rational in vars, r_i rational) generate a CYCLIC Kummer extension of the field of vars over its constants, i.e. one simple radical extension: with the b_i factored into irreducibles p_j of Q[vars], the exponent vectors r_i e_ij mod 1 generate a subgroup of (Q/Z)^J, cyclic iff its order equals its exponent. Returns True, False, or $Failed (a non-rational constant factor, or a group too large to enumerate). Constant factors are ignored here: they change the constant field (Unimplemented elsewhere), not the Kummer group over it.";
brKummerVectors::usage = "brKummerVectors[{{b1, r1}, ...}, vars] is the Kummer data of the radicals b_i^r_i (b_i rational in vars, r_i rational): with the b_i factored into the irreducibles p_j of Q[vars] (associates identified, constant factors dropped), it returns {vecs, nn}: vecs[[i]] the exponent vector r_i e_ij reduced mod 1, an element of (Q/Z)^J, and nn the lcm of their denominators. Constant factors change the constant field, not the Kummer group over it.";
brKummerVectors[rads_List, vars_List] := Module[{facs, ps, vecs},
  facs = Map[Function[rb, Join[tragerFactorList[Numerator[Together[rb[[1]]]]],
       Map[Function[fe, {fe[[1]], -fe[[2]]}], tragerFactorList[Denominator[Together[rb[[1]]]]]]]], rads];
  facs = Select[#, brDependsQ[First[#], vars] &] & /@ facs;
  ps = DeleteDuplicates[Flatten[Map[First, facs, {2}]]];
  (* identify associate factors (p and -p, or rational multiples) *)
  ps = DeleteDuplicates[ps, NumberQ[Together[#1/#2]] &];
  vecs = MapThread[
    Function[{fl, r}, Table[r Total[Cases[fl, {q_, e_} /; NumberQ[Together[q/p]] :> e]], {p, ps}]],
    {facs, Last /@ rads}];
  vecs = Mod[#, 1] & /@ vecs;
  {vecs, LCM @@ Denominator[Flatten[{vecs, 1}]]}];

brKummerCyclicQ[rads_List, vars_List] := Module[{vecs, nn, gens, seen, frontier, new, ordG, expG},
  {vecs, nn} = brKummerVectors[rads, vars];
  gens = Round[nn #] & /@ vecs;
  expG = LCM @@ (nn/GCD[nn, GCD @@ Append[#, nn]] & /@ gens);
  seen = <|ConstantArray[0, Length[First[vecs]]] -> True|>; frontier = Keys[seen];
  While[frontier =!= {},
    new = DeleteDuplicates[Flatten[Table[Mod[f + g, nn], {f, frontier}, {g, gens}], 1]];
    new = Select[new, !KeyExistsQ[seen, #] &];
    Scan[(seen[#] = True) &, new];
    If[Length[seen] > 100000, Return[$Failed, Module]];
    frontier = new];
  ordG = Length[seen];
  ordG === expG];

brCompositumPlan::usage = "brCompositumPlan[{{b1, d1}, ...}, vars] rewrites the radical atoms rho_i = b_i^(1/d_i) (b_i rational in vars, distinct pairs) as ONE radical, using nothing but the identities rho_i^d_i = b_i and (z^r)^m = z^(r m) (m an integer), so that every step holds pointwise on every branch. It searches integers c_i in [0, d_i) and k_j in [0, N), N the order of the Kummer group, with (a) c_i N/d_i an integer, (b) (c_i k_j - delta_ij)/d_i an integer for all i, j, and (c) the Kummer class of y = Prod rho_i^c_i of order N. Then y^N = B = Prod b_i^(c_i N/d_i) is rational, Y^N - B is irreducible, and rho_j = R_j y^k_j with R_j = Prod_i b_i^(-(c_i k_j - delta_ij)/d_i) rational: Q(vars)(rho_1, ...) = Q(vars)(y) as functions. It returns <|\"B\" -> B, \"N\" -> N, \"c\" -> {c_i}, \"k\" -> {k_j}, \"R\" -> {R_j}|>, or None when no such monomial exists (the radicals agree with one radical only up to a locally constant root of unity, e.g. Sqrt[x^2 (x + 1)] and Sqrt[x + 1], or a constant field extension, e.g. Sqrt[2 x] and Sqrt[x]) or the search is too large.";
brCompositumPlan[atoms_List, vars_List, halves_ : False] := Module[{ds = Last /@ atoms, bs = First /@ atoms, vecs, nn, ordN, cs, ks, ok, res = None, n = Length[atoms], eij, sq, eps},
  {vecs, nn} = brKummerVectors[{#[[1]], 1/#[[2]]} & /@ atoms, vars];
  If[Times @@ ds > 20000, Return[None]];
  (* the group order N: the subgroup generated by the vecs, enumerated *)
  ordN = Module[{gens = Round[nn #] & /@ vecs, seen, frontier, new},
    seen = <|ConstantArray[0, Length[First[vecs]]] -> True|>; frontier = Keys[seen];
    While[frontier =!= {},
      new = Select[DeleteDuplicates[Flatten[Table[Mod[f + g, nn], {f, frontier}, {g, gens}], 1]], !KeyExistsQ[seen, #] &];
      Scan[(seen[#] = True) &, new]; frontier = new];
    Length[seen]];
  Do[
    If[!AllTrue[Range[n], IntegerQ[cs[[#]] ordN/ds[[#]]] &], Continue[]];
    If[LCM @@ Denominator[Append[Mod[cs . vecs, 1], 1]] =!= ordN, Continue[]];
    If[!halves,
      ks = Table[SelectFirst[Range[0, ordN - 1],
              Function[k, AllTrue[Range[n], IntegerQ[(cs[[#]] k - If[# === j, 1, 0])/ds[[#]]] &]], None], {j, n}];
      If[MemberQ[ks, None], Continue[]];
      res = <|"B" -> Together[Product[bs[[i]]^(cs[[i]] ordN/ds[[i]]), {i, n}]], "N" -> ordN, "c" -> cs, "k" -> ks,
        "R" -> Table[Together[Product[bs[[i]]^(-(cs[[i]] ks[[j]] - If[i === j, 1, 0])/ds[[i]]), {i, n}]], {j, n}],
        "eps" -> Table[None, {n}]|>;
      Break[]];
    (* E6d: up to a sign.  e_ij = (c_i k_j - delta_ij)/d_i in Z/2; the product
       of the b_i with a half-integer e_ij must be a square s_j^2 in Q(vars), so
       Prod_(i in H_j) Sqrt[b_i] = eps_j s_j with eps_j = +-1 locally constant *)
    ks = Table[SelectFirst[Range[0, ordN - 1],
            Function[k, AllTrue[Range[n], IntegerQ[2 (cs[[#]] k - If[# === j, 1, 0])/ds[[#]]] &] &&
              brRationalSquareRoot[Product[If[IntegerQ[(cs[[i]] k - If[i === j, 1, 0])/ds[[i]]], 1, bs[[i]]], {i, n}], vars] =!= None],
            None], {j, n}];
    If[MemberQ[ks, None], Continue[]];
    eij = Table[(cs[[i]] ks[[j]] - If[i === j, 1, 0])/ds[[i]], {j, n}, {i, n}];
    eps = Table[With[{H = Select[Range[n], !IntegerQ[eij[[j, #]]] &]},
        If[H === {}, None, {H, brRationalSquareRoot[Product[bs[[i]], {i, H}], vars]}]], {j, n}];
    res = <|"B" -> Together[Product[bs[[i]]^(cs[[i]] ordN/ds[[i]]), {i, n}]], "N" -> ordN, "c" -> cs, "k" -> ks,
      "R" -> Table[Together[1/(Product[bs[[i]]^Floor[eij[[j, i]]], {i, n}] If[eps[[j]] === None, 1, eps[[j, 2]]])], {j, n}],
      "eps" -> eps|>;
    Break[],
    {cs, Tuples[Range[0, # - 1] & /@ ds]}];
  res];

brRationalSquareRoot::usage = "brRationalSquareRoot[e, vars] is s with s^2 = e, s in Q(vars), when e (rational in vars) is the square of one up to nothing: every irreducible factor of its numerator and denominator has an even exponent and the constant is a positive rational square; otherwise None (a constant like 2 or -1 would need a constant field extension).";
brRationalSquareRoot[e_, vars_List] := Module[{t = Together[e], fl, c = 1, s = 1, ok = True},
  fl = Join[tragerFactorList[Numerator[t]], {#[[1]], -#[[2]]} & /@ tragerFactorList[Denominator[t]]];
  Do[If[brDependsQ[q[[1]], vars],
      If[OddQ[q[[2]]], ok = False, s = s q[[1]]^(q[[2]]/2)],
      c = c q[[1]]^q[[2]]], {q, fl}];
  If[!ok || !(NumberQ[c] && c > 0) || !(IntegerQ[Sqrt[Numerator[c]]] && IntegerQ[Sqrt[Denominator[c]]]), Return[None]];
  Together[Sqrt[c] s]];

brCompositumRewrite::usage = "brCompositumRewrite[f, x] is the front-end rewrite of a radical compositum that is one simple radical in disguise (Tower:cyclic-compositum): when the radicals of the surface integrand f (radicands rational in x and the tower's Log/Exp kernels, none inside a kernel) generate a cyclic Kummer extension and brCompositumPlan finds a pointwise single-radical presentation, it returns <|\"g\" -> g, \"B\" -> B, \"N\" -> N, \"y\" -> y|>: g is f with every radical written as R_j B^(k_j/N) (an algebraic identity of the field Q(...)(y), y^N = B), and y = Prod rho_i^c_i is the product of the user's own radicals that B^(1/N) stands for. An antiderivative G of g maps to one of f by brCompositumBack (B^(r) -> y^(r N)); a certificate for g is one for f (the same differential field). Otherwise None.";
brCompositumRewrite[f_, x_Symbol] := Module[{ks, fs, monos, vars, pows, atoms, plan, Ysurf, g, y, toSurf, epsKeys, epsSyms, epsOf},
  If[!FreeQ[f, (Log[a_] | Power[E, a_]) /; brRadicalQ[a, {x}]], Return[None]];
  ks = brCatch[brKernelsToSymbols[f, x]];
  If[brFailQ[ks], Return[None]];
  {fs, monos} = ks;
  vars = Prepend[#["sym"] & /@ monos, x];
  pows = DeleteDuplicates[Cases[fs, Power[b_, r_Rational] /; brDependsQ[b, vars], {0, Infinity}]];
  If[Length[pows] < 2 || !FreeQ[First /@ pows, Power[_, _Rational]] ||
     !AllTrue[First /@ pows, brRationalInQ[#, vars] &], Return[None]];
  atoms = DeleteDuplicates[{Together[First[#]], Denominator[Last[#]]} & /@ pows,
    (Together[#1[[1]] - #2[[1]]] === 0 && #1[[2]] === #2[[2]]) &];
  If[Length[DeleteDuplicates[First /@ atoms, Together[#1 - #2] === 0 &]] < 2, Return[None]];
  If[brKummerCyclicQ[{#[[1]], 1/#[[2]]} & /@ atoms, vars] =!= True, Return[None]];
  plan = brCompositumPlan[atoms, vars];
  If[plan === None, plan = brCompositumPlan[atoms, vars, True]];
  If[plan === None, Return[None]];
  toSurf = brBack[monos];
  (* E6d: the sign of radical j is the locally constant eps_j = Prod_(i in H_j)
     Sqrt[b_i]/s_j, written as a symbol E_j (one per distinct H_j, s_j) *)
  epsKeys = DeleteDuplicates[DeleteCases[plan["eps"], None]];
  epsSyms = Table[Symbol["Bronstein`Private`brEps" <> ToString[i]], {i, Length[epsKeys]}];
  epsOf = Table[If[plan["eps"][[j]] === None, 1, epsSyms[[First[FirstPosition[epsKeys, plan["eps"][[j]]]]]]], {j, Length[atoms]}];
  Ysurf = Power[plan["B"] /. toSurf, 1/plan["N"]];
  g = fs /. Power[b_, r_Rational] /; brDependsQ[b, vars] :>
    With[{j = SelectFirst[Range[Length[atoms]],
               Together[atoms[[#, 1]] - b] === 0 && atoms[[#, 2]] === Denominator[r] &, None]},
      If[j === None, Power[b, r],
        (epsOf[[j]] (plan["R"][[j]] /. toSurf) Ysurf^plan["k"][[j]])^Numerator[r]]];
  y = Product[Power[atoms[[i, 1]] /. toSurf, 1/atoms[[i, 2]]]^plan["c"][[i]], {i, Length[atoms]}];
  <|"g" -> (g /. toSurf), "B" -> (plan["B"] /. toSurf), "N" -> plan["N"], "y" -> y,
    "signs" -> MapThread[{#1, (Product[Sqrt[atoms[[i, 1]]], {i, #2[[1]]}]/#2[[2]]) /. toSurf} &, {epsSyms, epsKeys}]|>];

brCompositumBack::usage = "brCompositumBack[rw, G] maps an expression G in the single radical B^(1/N) of the rewrite rw (brCompositumRewrite) back to the user's radicals: every B^r and (1/B)^r becomes y^(r N), y the product of the user's radicals with y^N = B. Since Y -> y is a differential embedding of Q(...)(Y), Y^N = B, the image of an antiderivative of g is an antiderivative of f.";
brCompositumBack[rw_Association, G_] := Module[{B = rw["B"], N0 = rw["N"], y = rw["y"]},
  G /. {Power[b_, r_Rational] /; Together[b - B] === 0 :> y^(r N0),
        Power[b_, r_Rational] /; Together[b B - 1] === 0 :> y^(-r N0)}];

brRadicalOverXQ::usage = "brRadicalOverXQ[e, x] is True when every radical b^r (r non-integer) in e has a radicand b rational in x alone (no kernel, no other symbol, no nested radical): e is then in Q(x)(y) for the radicals of e.";
brRadicalOverXQ[e_, x_Symbol] := AllTrue[
  Cases[e, Power[b_, _Rational] /; !FreeQ[b, x], {0, Infinity}], brRationalInQ[First[#], {x}] &];

brRationalInQ::usage = "brRationalInQ[e, vars] is True when e is a rational function of the symbols vars with coefficients in the constant field (rational, or Gaussian rational for Q(i); no other symbol, no radical, no kernel).";
brRationalInQ[e_, vars_List] := Module[{t = Together[e], num, den},
  num = Expand[Numerator[t]]; den = Expand[Denominator[t]];
  PolynomialQ[num, vars] && PolynomialQ[den, vars] &&
    Complement[Variables[{num, den}], vars] === {} &&
    FreeQ[{num, den}, Power[b_, _Rational] /; !(tragerAlgConstants =!= {} && rationalNumberQ[b])] &&
    FreeQ[{num, den}, _Real] && FreeQ[{num, den}, c_Complex /; !brGaussianRationalQ[c]]];

brTower::usage = "brTower[f, x] is the tower record of the integrand f (see the header of Tower.m): kernels replaced by monomials, the radical detected, the shape classified and, for Shape A, theta moved to the top. Throws a Decline outside the scope.";
brTower[f0_, x_Symbol] := Module[{f, fs, monos, vars, rad, shape, theta = None, kind = None, lower, bsyms, pos, above},
  f = f0 /. Exp[a_] :> E^a;
  brCheckHeads[f, x];
  {fs, monos} = brKernelsToSymbols[f, x];
  vars = Prepend[#["sym"] & /@ monos, x];
  If[!FreeQ[fs, brKernelPattern[]] && brDependsQ[Cases[fs, brKernelPattern[], {0, Infinity}], vars],
    throwInternal["brTowerKernelLeft", <|"f" -> fs|>]];
  rad = brDetectRadical[{fs, Select[monos, TrueQ[#["algebraic"]] &] /. a_Association :> a["arg"]}, vars];
  (* every radical power by its OWN placeholder: one symbol for Sqrt[t] and
     t^(3/2) collapsed -2 t Sqrt[t] + 2 t^(3/2) to 0 in a denominator
     (Power::infy, a wrong "not algebraic" decline; review 2026-11-14) *)
  If[!With[{rp = DeleteDuplicates[Cases[fs, Power[b_, r_Rational] /; brDependsQ[b, vars], {0, Infinity}]]},
       With[{ph = Table[Unique["Bronstein`Private`brRadPh"], {Length[rp]}]},
         brRationalInQ[fs /. Thread[rp -> ph], Join[vars, ph]]]],
    throwUnimplemented["Tower:not-algebraic", "integrand is not algebraic over the tower (irrational exponent or symbolic parameter)"]];
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
      throwUnimplemented["Tower:radical-mid-tower", "a monomial above the radicand's top monomial depends on it (radical in the middle of the tower)"]];
    monos = Join[Delete[monos, pos], {monos[[pos]]}];
    kind = Last[monos]["kind"];
    lower = Most[#["sym"] & /@ monos]];
  <|"x" -> x, "monos" -> monos, "syms" -> (#["sym"] & /@ monos), "kinds" -> (#["kind"] & /@ monos),
    "args" -> (#["arg"] & /@ monos), "rules" -> brRules[x, monos],
    "back" -> brBack[monos], "kernels" -> (#["kernel"] & /@ monos),
    "f" -> fs, "radical" -> rad, "shape" -> shape, "theta" -> theta, "kind" -> kind,
    "lower" -> lower|>];


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
