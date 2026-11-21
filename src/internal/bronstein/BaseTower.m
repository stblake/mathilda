(* ::Package:: *)

(* ========================================================================== *)
(* BaseTower.m -- Shape B over the curve with TWO monomials (Phase 5, M1):     *)
(* K = k0(t1)(t2), k0 = Q(x)(Y), Y^m = G(x), each t_i = Log[u_i] or E^(w_i)    *)
(* with u_i, w_i in k_(i-1) (arguments may involve Y).  The specification is  *)
(* docs/notes/bronstein_extensions.tex, sections E2 and E1, and its M1 scope. *)
(*                                                                            *)
(* Admission (E2.2; brDTAdmit, called by Tower.m's brKernelsToSymbols): a      *)
(* kernel is a new monomial iff Du/u (logarithm) or Dw (exponential) is NOT   *)
(* in the Q-span of the structure generators G_(i-1), proved by an            *)
(* inconsistent exact linear system on normal forms (brDFSpan; Prop E2span).  *)
(* A dependent kernel is eliminated ONLY with a zero offset (A3): a           *)
(* logarithm when kappa = 1 exactly (normal form) and the branch multiple m = *)
(* 0 at every sample point; an exponential when the constant c = 0 exactly    *)
(* and every coefficient is an integer.  Everything else declines.            *)
(*                                                                            *)
(* One level (Thm E2level): Hermite and Rothstein-Trager in t_j over          *)
(* k_(j-1) are BaseCurve.m's routines on a LEVEL RECORD (brDTLevel: the       *)
(* derivation of every monomial of the level, the constant test against all  *)
(* of them); the polynomial part is                                           *)
(*   exponential t_j: RDE_(j-1)(i Dw_j, p_i) for every i != 0 ([SI] 5.9.1),   *)
(*     RDE_0 = brAlgRDE, RDE_1 = brDTRDE ([SI] Ch. 6 over k0, E1.1);          *)
(*   logarithmic t_j: LI_(j-1)(p_i; Dt_j) top down ([SI] 5.8.1),             *)
(*     LI_0 = brBCLimited (E1.3), LI_1 = brDTLimited (integrate, then the    *)
(*     structure theorem over k1: Cor 9.3.1);                                 *)
(* and the t_j^0 coefficient recurses to level j - 1, ending in Trager.       *)
(*                                                                            *)
(* The level-0 oracles LD_0 / PLD_0 (Prop E1ld) are brDTPLD0: Trager's        *)
(* integral, its Liouville data (rational coefficients, logands in k0; else  *)
(* Decline: I2), then the span.  ROOT_0 (OL5) is not built: a witness that   *)
(* needs it declines.  Every "no" is one of [SI] Thms 5.6.1, 5.8.1, 5.9.1,   *)
(* Ch. 6, a Trager certificate, or a decided auxiliary (E2.4); every positive *)
(* auxiliary answer is checked exactly in normal form before it is used, and *)
(* the final antiderivative is verified by the pipeline.                     *)
(*                                                                            *)
(* Elements are expressions in x, Y and the monomials, in BaseCurve.m's       *)
(* normal form Sum c_i Y^i (brBCNorm); a polynomial in t over k0 is such an   *)
(* expression polynomial in t.                                                *)
(* ========================================================================== *)

(* -------------------------------------------------------------------------- *)
(* The curve, normal forms, admission.                                         *)
(* -------------------------------------------------------------------------- *)

brDTCurve::usage = "brDTCurve[b, n, x] is the curve record of the radical b^(1/n), b in Q(x): <|\"x\", \"Y\", \"m\" -> n, \"G\", \"b\", \"phi\", \"sY\", \"Ysurf\"|> with Y^m = G the reduced, absolutely irreducible curve (componentCurve of q^(n-1) p, b = p/q), the user's b^(1/n) = phi Y, D Y = sY Y, and Ysurf the surface expression of Y. Declines (BaseTower:curve-not-reduced) when Y^n - b is not absolutely irreducible over Q(x) (A0; a constant-field extension, OL2).";
brDTCurve[b_, n_Integer, x_Symbol] := Module[{bt = Together[b], p, q, g, cc, G, Y = Unique["Bronstein`Private`brYt"]},
  p = Numerator[bt]; q = Denominator[bt];
  g = Expand[q^(n - 1) p];
  cc = componentCurve[n, g, 1, x];
  If[!cc["exact"] || cc["m"] =!= n || !reducedCurveQ[n, cc["G"], x],
    throwUnimplemented["BaseTower:curve-not-reduced", "tower over the curve: Y^n - b is not absolutely irreducible over Q(x) (a constant-field extension, OL2)"]];
  G = Expand[cc["G"]];
  <|"x" -> x, "Y" -> Y, "m" -> n, "G" -> G, "b" -> bt, "phi" -> Together[cc["s"]/q],
    "sY" -> Together[D[G, x]/(n G)], "Ysurf" -> Together[q bt^(1/n)/cc["s"]]|>];

brDTCurveOf::usage = "brDTCurveOf[es, x] is the curve record (brDTCurve) of the radicals of the expressions es (kernel arguments): every power c^r (r non-integer, c depending on x) must have the SAME radicand c = b, rational in x over Q, and n is the lcm of the denominators of the r. Another spelling ((1/b)^r, several radicands: A1, OL6) declines as BaseTower:radical-spelling.";
brDTCurveOf[es_List, x_Symbol] := Module[{pows, bases, b},
  pows = DeleteDuplicates[Cases[es, Power[c_, r_Rational] /; !FreeQ[c, x], {0, Infinity}]];
  If[pows === {}, throwInternal["brDTCurveNoRadical", <|"es" -> es|>]];
  bases = Together[First[#]] & /@ pows;
  b = First[bases];
  If[!AllTrue[bases, Together[# - b] === 0 &],
    throwUnimplemented["BaseTower:radical-spelling", "a tower over the curve whose kernels spell the radical other than as powers of one radicand b ((1/b)^r, several radicands; A1, OL6)"]];
  If[!brDTRationalQ[b, {x}],
    throwUnimplemented["BaseTower:constants", "a constant other than a rational number in a tower over the curve (hypothesis H_Q; OL2)"]];
  brDTCurve[b, LCM @@ (Denominator[Last[#]] & /@ pows), x]];

brDTRationalQ::usage = "brDTRationalQ[e, vars] is True when e is a rational function of the symbols vars with RATIONAL coefficients: no other symbol, no complex, inexact or radical constant, no Pi, E, EulerGamma, Catalan, GoldenRatio, Degree and no Log (hypothesis H_Q of E2).";
brDTRationalQ[e_, vars_List] := Module[{t = Together[e], num, den},
  num = Expand[Numerator[t]]; den = Expand[Denominator[t]];
  PolynomialQ[num, vars] && PolynomialQ[den, vars] && Complement[Variables[{num, den}], vars] === {} &&
    FreeQ[{num, den}, _Complex | _Real | Power[_?NumericQ, _Rational] | Pi | E | EulerGamma | Catalan | GoldenRatio | Degree | Log]];

brDTToY::usage = "brDTToY[cv, e] writes the expression e (in x, the tower symbols and powers b^r of the curve's radicand) with every b^r replaced by (phi Y)^(r m): an exact identity, since b^(k/n) = (b^(1/n))^k for the principal branch (A1). A power of another radicand, or one whose exponent is not a multiple of 1/m, declines (BaseTower:radical-spelling).";
brDTToY[cv_Association, e_] := e /. Power[c_, r_Rational] /; !FreeQ[c, cv["x"]] :>
  If[Together[c - cv["b"]] === 0 && IntegerQ[r cv["m"]], (cv["phi"] cv["Y"])^(r cv["m"]),
    throwUnimplemented["BaseTower:radical-spelling", "a tower over the curve whose kernels spell the radical other than as powers of one radicand b ((1/b)^r, several radicands; A1, OL6)"]];

brDTMonoData::usage = "brDTMonoData[cv, monos] is, for every monomial of monos in tower order, <|\"argY\" -> its argument in normal form, \"DY\" -> D t in normal form, \"eta\" -> the structure generator (D t for a logarithm, D w for an exponential E^w)|>, the derivation being D x = 1, D Y = sY Y and the D t of the monomials before it (E2.1).";
brDTMonoData[cv_Association, monos_List] := Module[{data = {}, rec, aY, Da},
  Do[
    rec = Join[cv, <|"lsyms" -> (#["sym"] & /@ Take[monos, i - 1]), "lDs" -> (#["DY"] & /@ data)|>];
    aY = brBCNorm[cv, brDTToY[cv, monos[[i]]["arg"]]];
    Da = brBCD[rec, aY];
    AppendTo[data, If[monos[[i]]["kind"] === "Log",
      With[{d = brBCNorm[cv, Da/aY]}, <|"argY" -> aY, "DY" -> d, "eta" -> d|>],
      <|"argY" -> aY, "DY" -> brBCNorm[cv, monos[[i]]["sym"] Da], "eta" -> Da|>]],
    {i, Length[monos]}];
  data];

brDFSpan::usage = "brDFSpan[cv, theta, gens, syms] decides theta in span_Q gens (theta, gens in normal form over the curve cv, rational in x, Y and the monomials syms): {r_1, ...} with theta = Sum r_i g_i (checked exactly in normal form), False when the exact linear system is inconsistent (a proof, Prop E2span: a normal form of Y-degree < m is zero iff each coordinate is), or $Failed when the solve did not decide (A2: a Decline, never an admission). An empty gens spans only 0.";
brDFSpan[cv_Association, theta_, gens_List, syms_List] := Module[{th = brBCNorm[cv, theta], r},
  If[th === 0, throwInternal["brDFSpanZero", <|"theta" -> theta|>]];
  If[gens === {}, Return[False]];
  r = Risch`RationalSpan[th, brBCNorm[cv, #] & /@ gens, Join[{cv["x"], cv["Y"]}, syms]];
  Which[
    r === False, False,
    ListQ[r] && Length[r] === Length[gens] && AllTrue[r, rationalNumberQ] && brBCNorm[cv, th - r . gens] === 0, r,
    True, $Failed]];

brDTZeroOffsetQ::usage = "brDTZeroOffsetQ[delta, x] is True when the locally constant surface expression delta (a branch difference Log[u] - Sum r_g theta_g with kappa = 1, so delta is in 2 Pi I Z/d) vanishes at every sample point (three real and four complex rationals, 30 digits): at least three points give a finite value and every finite value is below 10^-10 (integer separation, A3). A point where delta is not finite is skipped.";
brDTZeroOffsetQ[delta_, x_Symbol] := Module[{pts = {3/7, 5/4, 7/3, 1/3 + 2 I/5, -3/2 + 7 I/10, -5/3 - I/4, 2 - 9 I/5}, vals},
  vals = Quiet[N[delta /. x -> #, 30] & /@ pts];
  vals = Select[vals, NumericQ[#] && FreeQ[#, Indeterminate | ComplexInfinity | DirectedInfinity] &];
  Length[vals] >= 3 && AllTrue[vals, Abs[#] < 10^-10 &]];

brDTAdmit::usage = "brDTAdmit[x, monos, kind, arg, k, cv] admits the kernel k = Log[arg] or E^arg (kind \"Log\" | \"Exp\", arg in x, the symbols of monos and powers of the curve's radicand) over the tower monos on the curve cv (E2.2). It returns {\"new\", monos'} (monos' = monos with the new monomial, flagged \"curve\"; brDFSpan proved Du/u resp. Dw outside span_Q G), or {\"rewrite\", e} with k = e an identity on the sampled regions: a dependent logarithm with kappa = 1 exactly and zero offset at every sample point (e = Sum r_g theta_g), a dependent exponential with c = 0 and integer coefficients (e = Prod u_g^r_g Prod t_g^r_g). Declines: Tower:dependent-log-offset (kappa != 1), Tower:dependent-log-branch (a nonzero 2 Pi I m/d), Tower:dependent-exp-offset (c != 0), Tower:dependent-exp-fractional (a non-integer coefficient), Tower:three-monomials, Tower:span-undecided.";
brDTAdmit[x_Symbol, monos_List, kind_String, arg_, k_, cv_Association] := Module[
  {syms = #["sym"] & /@ monos, data, rec, aY, Da, theta, sp, m2, t, d, kap, L0, delta, c},
  data = brDTMonoData[cv, monos];
  rec = Join[cv, <|"lsyms" -> syms, "lDs" -> (#["DY"] & /@ data)|>];
  aY = brBCNorm[cv, brDTToY[cv, arg]];
  Da = brBCD[rec, aY];
  theta = If[kind === "Log", brBCNorm[cv, Da/aY], Da];
  If[theta === 0, throwInternal["brDTAdmitConstant", <|"kernel" -> k|>]];
  sp = brDFSpan[cv, theta, #["eta"] & /@ data, syms];
  Which[
    sp === False,
      If[Length[monos] >= 2,
        throwUnimplemented["Tower:three-monomials", "three or more monomials over the curve (Phase 5 M1 takes two)"]];
      {m2, t} = brAddMonomial[x, monos, kind, Together[arg], k /. brBack[monos]];
      m2[[-1]] = Join[m2[[-1]], <|"algebraic" -> brRadicalQ[arg, Prepend[syms, x]], "curve" -> True|>];
      {"new", m2},
    !ListQ[sp],
      throwUnimplemented["Tower:span-undecided", "the structure-theorem span over the curve was not decided (A2)"],
    kind === "Log",
      (* t = Sum r_g theta_g + delta, delta = (Log kappa + 2 Pi I m)/d *)
      d = LCM @@ (Denominator /@ sp);
      kap = brBCNorm[cv, aY^d Product[If[monos[[g]]["kind"] === "Log", data[[g]]["argY"], syms[[g]]]^(-d sp[[g]]),
        {g, Length[monos]}]];
      If[!FreeQ[kap, Alternatives @@ Join[{x, cv["Y"]}, syms]], throwInternal["brDTKappaNotConstant", <|"kappa" -> kap|>]];
      If[kap =!= 1,
        throwUnimplemented["Tower:dependent-log-offset", "a dependent logarithm over the curve whose constant offset Log[kappa] is nonzero (kappa != 1: the constant field would leave Q; A3, OL2)"]];
      L0 = Sum[sp[[g]] If[monos[[g]]["kind"] === "Log", syms[[g]], monos[[g]]["arg"]], {g, Length[monos]}];
      delta = (k /. brBack[monos]) - (L0 /. brBack[monos]);
      If[!brDTZeroOffsetQ[delta, x],
        throwUnimplemented["Tower:dependent-log-branch", "a dependent logarithm over the curve that differs from its reduction by a nonzero multiple of 2 Pi I/d on a sampled region (A3)"]];
      {"rewrite", L0},
    True,
      c = brBCNorm[cv, aY - Sum[sp[[g]] If[monos[[g]]["kind"] === "Log", syms[[g]], data[[g]]["argY"]], {g, Length[monos]}]];
      If[!FreeQ[c, Alternatives @@ Join[{x, cv["Y"]}, syms]], throwInternal["brDTExpOffsetNotConstant", <|"c" -> c|>]];
      If[c =!= 0,
        throwUnimplemented["Tower:dependent-exp-offset", "a dependent exponential over the curve with a nonzero constant offset (a factor E^c: the constant field would leave Q; A3, OL2)"]];
      If[!AllTrue[sp, IntegerQ],
        throwUnimplemented["Tower:dependent-exp-fractional", "a dependent exponential over the curve with a fractional coefficient (an algebraic element; A3, OL6)"]];
      {"rewrite", Product[If[monos[[g]]["kind"] === "Log", monos[[g]]["arg"]^sp[[g]], syms[[g]]^sp[[g]]], {g, Length[monos]}]}]];

(* -------------------------------------------------------------------------- *)
(* The field record and its levels.                                            *)
(* -------------------------------------------------------------------------- *)

brDTContext::usage = "brDTContext[T] is the field record of the tower record T (shape \"BTower\"): the curve (brDTCurve of the radical) joined with <|\"B\" -> curveBasis, \"syms\", \"kinds\", \"argsY\", \"Ds\", \"etas\" (normal forms, brDTMonoData), \"back\" -> {t_i -> kernel}, \"kernels\", \"f\" -> the integrand in normal form|>. Declines (BaseTower:constants) outside hypothesis H_Q: a Gaussian, radical or symbolic constant, or any constant other than a rational in the integrand or a kernel argument.";
brDTContext[T_Association] := Module[{x = T["x"], monos = T["monos"], syms = T["syms"], rad = T["radical"], cv, fY, data, vars},
  If[TrueQ[tragerConstantI] || tragerAlgConstants =!= {} || tragerParams =!= {},
    throwUnimplemented["BaseTower:constants", "a constant other than a rational number in a tower over the curve (hypothesis H_Q; OL2)"]];
  cv = brDTCurve[rad["b"], rad["n"], x];
  vars = Join[{x, cv["Y"]}, syms];
  fY = brBCNorm[cv, brDTToY[cv, T["f"]]];
  data = brDTMonoData[cv, monos];
  If[!brDTRationalQ[fY, vars] || !AllTrue[data, brDTRationalQ[#["argY"], vars] &],
    throwUnimplemented["BaseTower:constants", "a constant other than a rational number in a tower over the curve (hypothesis H_Q; OL2)"]];
  Join[cv, <|"B" -> curveBasis[cv["m"], cv["G"], x], "syms" -> syms, "kinds" -> (#["kind"] & /@ monos),
    "argsY" -> (#["argY"] & /@ data), "Ds" -> (#["DY"] & /@ data), "etas" -> (#["eta"] & /@ data),
    "back" -> brBack[monos], "kernels" -> (#["kernel"] & /@ monos), "f" -> fY|>]];

brDTLevel::usage = "brDTLevel[C, j] is the level-j record of the field record C: k_j = k0(t_1..t_j) with t_j as BaseCurve.m's theta (\"th\", \"kind\", \"Dth\" = D t_j, \"Dh\" = D of its argument, \"h\", \"kernel\"), the coefficient field \"F\" -> fieldFn[{t_1..t_j}] and the level's derivation \"lsyms\" -> {t_1..t_j}, \"lDs\" -> {D t_1..D t_j} (brBCD differentiates every monomial; brBCLogPart tests constants against every one). Level 0 is k0 itself (\"th\" -> None).";
brDTLevel[C_Association, j_Integer] := Module[{syms = Take[C["syms"], j]},
  Join[C, <|"level" -> j, "lsyms" -> syms, "lDs" -> Take[C["Ds"], j], "F" -> brFn[syms],
    "th" -> If[j >= 1, C["syms"][[j]], None], "kind" -> If[j >= 1, C["kinds"][[j]], None],
    "Dth" -> If[j >= 1, C["Ds"][[j]], 0], "h" -> If[j >= 1, C["argsY"][[j]], None],
    "Dh" -> If[j >= 1, If[C["kinds"][[j]] === "Exp", C["etas"][[j]], brBCNorm[C, C["Ds"][[j]] C["argsY"][[j]]]], 0],
    "kernel" -> If[j >= 1, C["kernels"][[j]], None]|>]];

brDTSurface::usage = "brDTSurface[C, e] writes e (in x, Y and the monomials) in the user's surface form: Y -> Ysurf, every monomial -> its kernel.";
brDTSurface[C_Association, e_] := e /. C["Y"] -> C["Ysurf"] /. C["back"];

(* -------------------------------------------------------------------------- *)
(* INT_j: one level of integration (Thm E2level), recursing to level 0.         *)
(* -------------------------------------------------------------------------- *)

brDTInt::usage = "brDTInt[C, j, e] integrates e in k_j (Thm E2level): Hermite and Rothstein-Trager in t_j over k_(j-1) (BaseCurve.m on brDTLevel[C, j]), the polynomial part (exponential t_j: RDE_(j-1)(i D w_j, p_i) for i != 0; logarithmic t_j: LI_(j-1) top down) and the t_j^0 coefficient by brDTInt at level j - 1; level 0 is Trager (brDTIntegrateK0). It returns <|\"alg\" -> v in k_j, \"logs\" -> {{c, u}, ...} (c rational, u in k_j), \"rs\" -> {RootSum terms}, \"surf\" -> the level-0 integral (surface form)|>, the antiderivative being v + Sum c Log[u] + Sum rs + surf. Throws NonElementary only on a certificate (E2.4 (a)-(e)).";
brDTInt[C_Association, 0, e_] := <|"alg" -> 0, "logs" -> {}, "rs" -> {}, "surf" -> brDTIntegrateK0[C, e]|>;
brDTInt[C_Association, j_Integer, e0_] /; j >= 1 := Module[
  {Cj = brDTLevel[C, j], th, e, her, sp, lp, L, alg = 0, low, logs, rs, eta, beta = 0, es, p, vc, q, e00 = 0},
  th = Cj["th"];
  e = brBCNorm[Cj, e0];
  If[e === 0, Return[<|"alg" -> 0, "logs" -> {}, "rs" -> {}, "surf" -> 0|>]];
  her = brBCHermite[Cj, e];
  sp = brBCSplit[Cj, her["R"]];
  lp = brBCLogPart[Cj, sp["a"], sp["d"]];
  L = Normal[Merge[{Association[Rule @@@ sp["L"]], Association[Rule @@@ brBCLaurent[Cj, lp["Rem"]]]},
    brBCNorm[Cj, Total[#]] &]];
  L = Select[List @@@ L, #[[2]] =!= 0 &];
  logs = Select[lp["Logs"], Length[#] === 2 &];
  rs = #[[2]] & /@ Select[lp["Logs"], Length[#] === 3 &];
  If[Cj["kind"] === "Exp",
    Do[If[l[[1]] === 0, e00 = l[[2]],
        q = If[j === 1, brDTRDE0[C, brBCNorm[Cj, l[[1]] Cj["Dh"]], l[[2]]],
                        brDTRDE[C, brBCNorm[Cj, l[[1]] Cj["Dh"]], l[[2]]]];
        If[First[q] === "NoSolution",
          throwNonElementary["no algebraic part: the Risch DE D b + i w' b = p_i has no solution over the tower on the curve ([SI] Thm 5.9.1, E1)"]];
        alg = alg + q[[2]] th^l[[1]]],
      {l, L}];
    low = brDTInt[C, j - 1, e00],
    eta = Cj["Dth"];
    If[L =!= {} && Min[First /@ L] < 0, throwInternal["brDTLogNegativePower", <|"L" -> L|>]];
    es = Association[Rule @@@ L];
    p = If[L === {}, 0, Max[First /@ L]];
    Do[
      vc = If[j === 1, brBCLimited[Cj, brBCNorm[Cj, Lookup[es, i, 0] - (i + 1) beta eta]],
                       brDTLimited[C, brBCNorm[Cj, Lookup[es, i, 0] - (i + 1) beta eta]]];
      alg = alg + brBCNorm[Cj, beta + vc[[2]]/(i + 1)] th^(i + 1);
      beta = vc[[1]],
      {i, p, 1, -1}];
    alg = alg + beta th;
    low = brDTInt[C, j - 1, brBCNorm[Cj, Lookup[es, 0, 0] - beta eta]]];
  <|"alg" -> her["H"] + alg + low["alg"], "logs" -> Join[logs, low["logs"]], "rs" -> Join[rs, low["rs"]],
    "surf" -> low["surf"]|>];

brDTIntegrateK0::usage = "brDTIntegrateK0[C, e] is INT_0: an elementary integral of e in k0 (free of every monomial) as a surface expression, by Trager's algorithm (brBCIntegrateK0: a Trager certificate is NonElementary, its decline a Decline).";
brDTIntegrateK0[C_Association, e_] := (
  If[!FreeQ[e, Alternatives @@ C["syms"]], throwInternal["brDTIntegrateK0Tower", <|"e" -> e|>]];
  brBCIntegrateK0[brDTLevel[C, 1], e]);

brDTLiouville0::usage = "brDTLiouville0[C, F0] reads the level-0 integral F0 (surface form, Trager's or the rational integrator's answer) as Liouville data {v, {{c, u}, ...}}: v in k0 and u in k0 in normal form, c RATIONAL (Prop E1ld(ii) needs both). An irrational coefficient declines (BaseTower:irrational-log, gap I2); a RootSum, ArcTan, a logand outside k0 over Q or any other shape declines (BaseTower:rootsum-logand).";
brDTLiouville0[C_Association, F0_] := Module[{C1 = brDTLevel[C, 1], x = C["x"], terms, v = 0, logs = {}, vars = {C["x"], C["Y"]}, u, vY},
  terms = Expand[F0];
  terms = If[Head[terms] === Plus, List @@ terms, {terms}];
  Do[Which[
      FreeQ[tm, Log | ArcTan | ArcTanh | ArcCot | RootSum | Root], v = v + tm,
      MatchQ[tm, Log[_]], AppendTo[logs, {1, First[tm]}],
      MatchQ[tm, Times[c_, Log[_]] /; FreeQ[c, x]],
        With[{c = tm /. Times[c0_, Log[_]] :> c0, a = tm /. Times[_, Log[a0_]] :> a0},
          If[!rationalNumberQ[c],
            throwUnimplemented["BaseTower:irrational-log", "an irrational logarithm coefficient in the integral over k0 inside LD/PLD/LI (gap I2)"]];
          AppendTo[logs, {c, a}]],
      True, throwUnimplemented["BaseTower:rootsum-logand", "a logarithmic part of an integral over the curve not given as rational multiples of logarithms of elements of the field (RootSum, ArcTan)"]],
    {tm, terms}];
  vY = brBCFromSurface[C1, v];
  u = {#[[1]], brBCFromSurface[C1, #[[2]]]} & /@ logs;
  If[!brDTRationalQ[vY, vars] || !AllTrue[u, brDTRationalQ[#[[2]], vars] && #[[2]] =!= 0 &],
    throwUnimplemented["BaseTower:rootsum-logand", "a logarithmic part of an integral over the curve not given as rational multiples of logarithms of elements of the field (RootSum, ArcTan)"]];
  {vY, u}];

(* -------------------------------------------------------------------------- *)
(* Level-0 oracles.                                                            *)
(* -------------------------------------------------------------------------- *)

brDTPLD0::usage = "brDTPLD0[C, alpha, eta] decides PLD_0 (Prop E1ld, j = 0) for alpha, eta in k0 (eta = D w, or 0 for LD_0): {\"none\"} when there are no n >= 1, m in Z, U in k0* with n alpha = m eta + DU/U (a proof: Trager proved alpha non-integrable, or the span system with n = 1 is inconsistent), or {\"radical\", n0, m0, {{e_l, u_l}, ...}} with U = Prod u_l^e_l, checked exactly in normal form. Declines when Trager declines, on Liouville data outside I2 (brDTLiouville0) or an undecided span (BaseTower:pld-undecided).";
brDTPLD0[C_Association, alpha_, eta_] := Module[{C0 = brDTLevel[C, 0], a = brBCNorm[C, alpha], F0, v, logs, gens, sp, mq, rq, n0},
  If[a === 0, Return[{"radical", 1, 0, {}}]];
  F0 = brCatch[brDTIntegrateK0[C, a]];
  If[brFailQ[F0], If[First[brFailTag[F0]] === "NonElementary", Return[{"none"}, Module], Throw[F0, "TragerFail"]]];
  {v, logs} = brDTLiouville0[C, F0];
  (* the Liouville data must reproduce alpha exactly before a "none" rests on it *)
  If[brBCNorm[C, a - brBCD[C0, v] - Sum[l[[1]] brBCD[C0, l[[2]]]/l[[2]], {l, logs}]] =!= 0,
    throwInternal["brDTLiouvilleCheck", <|"alpha" -> alpha|>]];
  gens = Join[If[eta === 0, {}, {brBCNorm[C, eta]}], brBCNorm[C, brBCD[C0, #[[2]]]/#[[2]]] & /@ logs];
  If[gens === {}, Return[{"none"}]];
  sp = Risch`RationalSpan[a, gens, {C["x"], C["Y"]}];
  If[sp === False, Return[{"none"}]];
  If[!ListQ[sp] || Length[sp] =!= Length[gens] || !AllTrue[sp, rationalNumberQ],
    throwUnimplemented["BaseTower:pld-undecided", "the span system of a (parametric) logarithmic-derivative problem over k0 was not decided"]];
  mq = If[eta === 0, 0, First[sp]]; rq = If[eta === 0, sp, Rest[sp]];
  n0 = LCM @@ (Denominator /@ Join[{mq, 1}, rq]);
  If[brBCNorm[C, n0 a - n0 (sp . gens)] =!= 0, throwInternal["brDTPLDCheck", <|"alpha" -> alpha|>]];
  {"radical", n0, n0 mq, Transpose[{n0 rq, #[[2]] & /@ logs}]}];

brDTIsDerivK0::usage = "brDTIsDerivK0[C, e] decides whether e in k0 is a derivative D v, v in k0: {\"yes\", v} or {\"no\"}. Trager's integral with its logarithmic part dropped gives v; e - D v = Sum c_l D v_l/v_l is then a derivative only if it is 0 (E1.3: reduce the c_l to a Q-basis, each component has zero residues everywhere), so the test e = D v decides. A Trager certificate is a \"no\"; a decline propagates.";
brDTIsDerivK0[C_Association, e_] := Module[{C1 = brDTLevel[C, 1], C0 = brDTLevel[C, 0], F0, v},
  If[brBCNorm[C, e] === 0, Return[{"yes", 0}]];
  F0 = brCatch[brDTIntegrateK0[C, e]];
  If[brFailQ[F0], If[First[brFailTag[F0]] === "NonElementary", Return[{"no"}, Module], Throw[F0, "TragerFail"]]];
  v = brBCFromSurface[C1, brBCLogFree[F0]];
  If[brBCNorm[C, e - brBCD[C0, v]] === 0, {"yes", v}, {"no"}]];

brDTRDE0::usage = "brDTRDE0[C, f, g] is RDE_0: D q + f q = g over k0 = Q(x)(Y) (f, g free of the monomials) by brAlgRDE ([B] Part 3; Lemmas E1choice, E1inf): {\"Solution\", q} (q in normal form), {\"NoSolution\"} (a proof) or a Decline (RDEAlgebraic:undecided, gap I1).";
brDTRDE0[C_Association, f_, g_] := Module[{B = C["B"], Y = C["Y"], fn = brBCNorm[C, f], gn = brBCNorm[C, g], r},
  If[!FreeQ[{fn, gn}, Alternatives @@ C["syms"]], throwInternal["brDTRDE0Tower", <|"f" -> f, "g" -> g|>]];
  If[gn === 0, Return[{"Solution", 0}]];
  r = brAlgRDE[B, afFromY[fieldQ[], B, fn, Y], afFromY[fieldQ[], B, gn, Y]];
  Switch[First[r],
    "Solution", {"Solution", brBCNorm[C, afToY[B, r[[2]], Y]]},
    "NoSolution", {"NoSolution"},
    _, throwUnimplemented["RDEAlgebraic:undecided", "Risch DE on the curve undecided: " <> ToString[Last[r]]]]];

(* -------------------------------------------------------------------------- *)
(* LI_1: limited integration over k1 against D t2, t2 a logarithm.             *)
(* -------------------------------------------------------------------------- *)

brDTLimited::usage = "brDTLimited[C, p] solves LI_1(p; eta), eta = D t_2 (t_2 = Log[u_2], u_2 in k1): {w, c} with p = D w + c eta, w in k1, c rational. INT_1 gives p = D V + Sum c_l D u_l/u_l (rational c_l, u_l in k1; else Decline); a solution exists iff Sum c_l D u_l/u_l - c eta is in span_Q G_1 for some c ([SI] Cor 9.3.1 over k1, constants Q), decided by Risch`RationalSpan; then w = V + r theta_1 (theta_1 = t_1 or w_1), checked exactly. Throws NonElementary when there is none (p not elementary over k1, or the span system inconsistent): [SI] Thm 5.8.1.";
brDTLimited[C_Association, p_] := Module[{C1 = brDTLevel[C, 1], C2 = brDTLevel[C, 2], eta, r, v0, logs0, logs, theta, sp, w, th1},
  eta = C2["Dth"];
  If[brBCNorm[C1, p] === 0, Return[{0, 0}]];
  r = brCatch[brDTInt[C, 1, p]];
  If[brFailQ[r], If[First[brFailTag[r]] === "NonElementary",
    throwNonElementary["no algebraic part: the limited integration problem has no solution over the tower on the curve (its integrand is not elementary over k1; [SI] Thm 5.8.1)"],
    Throw[r, "TragerFail"]]];
  If[r["rs"] =!= {}, throwUnimplemented["BaseTower:rootsum-logand", "a logarithmic part of an integral over the curve not given as rational multiples of logarithms of elements of the field (RootSum, ArcTan)"]];
  {v0, logs0} = brDTLiouville0[C, r["surf"]];
  logs = Join[r["logs"], logs0];
  theta = brBCNorm[C1, Sum[l[[1]] brBCD[C1, l[[2]]]/l[[2]], {l, logs}]];
  (* the Liouville data must reproduce p exactly before a "no" rests on it *)
  If[brBCNorm[C1, p - brBCD[C1, r["alg"] + v0] - theta] =!= 0, throwInternal["brDTLiouvilleCheck", <|"p" -> p|>]];
  th1 = If[C["kinds"][[1]] === "Log", C["syms"][[1]], C["argsY"][[1]]];
  sp = If[theta === 0, {0, 0}, Risch`RationalSpan[theta, {eta, C["etas"][[1]]}, {C["x"], C["Y"], C["syms"][[1]]}]];
  If[sp === False,
    throwNonElementary["no algebraic part: the limited integration problem has no solution over the tower on the curve (the logarithmic part of the integral over k1 is not in the span of the structure generators; [SI] Thm 5.8.1, Cor 9.3.1)"]];
  If[!ListQ[sp] || Length[sp] =!= 2 || !AllTrue[sp, rationalNumberQ],
    throwUnimplemented["BaseTower:li-undecided", "the span system of a limited integration problem over k1 was not decided"]];
  w = brBCNorm[C1, r["alg"] + v0 + sp[[2]] th1];
  If[brBCNorm[C1, p - brBCD[C1, w] - sp[[1]] eta] =!= 0, throwInternal["brDTLimitedCheck", <|"p" -> p|>]];
  {w, sp[[1]]}];

(* -------------------------------------------------------------------------- *)
(* RDE_1: [SI] Ch. 6 over k0, t = t1 primitive or exponential (E1.1).          *)
(* Polynomials in t over k0 are normal-form expressions polynomial in t.       *)
(* -------------------------------------------------------------------------- *)

brDTDeg::usage = "brDTDeg[C1, e] is the degree in t = C1[\"th\"] of e in k0[t] (normal form), -Infinity for 0; an internal error when a coordinate denominator involves t.";
brDTDeg[C1_Association, e_] := Module[{t = C1["th"], n = brBCNorm[C1, e], cs},
  If[n === 0, Return[-Infinity]];
  cs = Together[Coefficient[n, C1["Y"], #]] & /@ Range[0, C1["m"] - 1];
  If[!FreeQ[Denominator /@ cs, t], throwInternal["brDTNotPolynomial", <|"e" -> e|>]];
  Max[If[# === 0, -Infinity, Exponent[Numerator[#], t]] & /@ cs]];

brDTCoeff::usage = "brDTCoeff[C1, e, k] is the coefficient of t^k of e in k0[t] (an element of k0 in normal form).";
brDTCoeff[C1_Association, e_, k_Integer] := Module[{t = C1["th"], Y = C1["Y"], n = brBCNorm[C1, e]},
  brBCNorm[C1, Sum[With[{c = Together[Coefficient[n, Y, i]]},
      Coefficient[Expand[Numerator[c]], t, k]/Denominator[c]] Y^i, {i, 0, C1["m"] - 1}]]];

brDTDivide::usage = "brDTDivide[C1, a, b] is {q, r} with a = q b + r, deg r < deg b, in k0[t] (long division, leading coefficients inverted in k0).";
brDTDivide[C1_Association, a_, b_] := Module[{t = C1["th"], db = brDTDeg[C1, b], lcInv, q = 0, r = brBCNorm[C1, a], dr, term, guard = 0},
  If[db === -Infinity, throwInternal["brDTDivideByZero"]];
  lcInv = brBCNorm[C1, 1/brDTCoeff[C1, b, db]];
  While[(dr = brDTDeg[C1, r]) >= db,
    If[++guard > 500, throwInternal["brDTDivideNoProgress"]];
    term = brBCNorm[C1, brDTCoeff[C1, r, dr] lcInv] t^(dr - db);
    q = q + term;
    r = brBCNorm[C1, r - term b];
    If[brDTDeg[C1, r] >= dr, throwInternal["brDTDivideLeadingTerm", <|"r" -> r|>]]];
  {brBCNorm[C1, q], r}];

brDTNormalDen::usage = "brDTNormalDen[C1, e] is the normal part of the t-denominator of e (a polynomial in x and t over Q, free of Y): the t-dependent factors of its normal-form denominator, without the special factor t for an exponential t. Since 1, Y, ..., Y^(m-1) is a basis of k0[t] over Q(x)[t], an irreducible rho of Q(x)[t] divides it with multiplicity e exactly when max over the places pi | rho of -nu_pi(e) is e (rho is squarefree over k0).";
brDTNormalDen[C1_Association, e_] := Module[{t = C1["th"], d},
  d = Last[brBCSplitDen[C1, e]];
  Expand[Times @@ (#[[1]]^#[[2]] & /@ Select[tragerFactorList[d],
    !FreeQ[#[[1]], t] && !(C1["kind"] === "Exp" && NumberQ[Together[#[[1]]/t]]) &])]];

brDTDegreeBound::usage = "brDTDegreeBound[C, a, b, c] is a valid bound n on deg q for a D q + b q = c in k0[t] ([SI] Sec. 6.3, E1.1 table): primitive t: max(0, dc - db) when db > da, else max(0, dc - da + 1), raised to c when db = da - 1 and LI_0(-lc b/lc a; Dt) has c a positive integer, or when db = da, LD_0(alpha) has a witness z (n0 = 1; n0 > 1 declines, ROOT_0) and LI_0(beta; Dt) has c a positive integer (beta = -(a_(da-1) Dz + b_(da-1) z)/(z lc a)); exponential t: max(0, dc - max(da, db)), raised to m0/n0 when da = db and PLD_0(-lc b/lc a; Dw) has a radical solution with m0/n0 an integer (remark (2) of E1.1). A non-integer c or ratio is discarded, never rounded.";
brDTDegreeBound[C_Association, a_, b_, c_] := Module[{C1 = brDTLevel[C, 1], da, db, dc, n, alpha, vc, ld, z, beta, pld},
  da = brDTDeg[C1, a]; db = brDTDeg[C1, b]; dc = brDTDeg[C1, c];
  alpha = If[db === -Infinity, 0, brBCNorm[C1, -brDTCoeff[C1, b, db]/brDTCoeff[C1, a, da]]];
  If[C1["kind"] === "Log",
    n = If[db > da, Max[0, dc - db], Max[0, dc - da + 1]];
    If[db === da - 1,
      vc = brCatch[brBCLimited[C1, alpha]];
      If[brFailQ[vc], If[First[brFailTag[vc]] =!= "NonElementary", Throw[vc, "TragerFail"]],
        If[IntegerQ[vc[[2]]] && vc[[2]] > 0, n = Max[n, vc[[2]]]]]];
    If[db === da,
      ld = brDTPLD0[C, alpha, 0];
      If[First[ld] === "radical",
        If[ld[[2]] =!= 1, throwUnimplemented["BaseTower:root0", "a logarithmic-derivative witness that needs an n-th root in k0 (ROOT_0, OL5)"]];
        z = brBCNorm[C1, Product[u[[2]]^u[[1]], {u, ld[[4]]}]];
        beta = brBCNorm[C1, -(brDTCoeff[C1, a, da - 1] brBCD[brDTLevel[C, 0], z] + brDTCoeff[C1, b, da - 1] z)/(z brDTCoeff[C1, a, da])];
        vc = brCatch[brBCLimited[C1, beta]];
        If[brFailQ[vc], If[First[brFailTag[vc]] =!= "NonElementary", Throw[vc, "TragerFail"]],
          If[IntegerQ[vc[[2]]] && vc[[2]] > 0, n = Max[n, vc[[2]]]]]]],
    n = Max[0, dc - Max[db, da]];
    If[da === db,
      pld = brDTPLD0[C, alpha, C1["Dh"]];
      If[First[pld] === "radical" && IntegerQ[pld[[3]]/pld[[2]]], n = Max[n, pld[[3]]/pld[[2]]]]]];
  n];

brDTSPDE::usage = "brDTSPDE[C, a, b, c, n] is Rothstein's SPDE ([SI] Sec. 6.4) for a D q + b q = c, deg q <= n, in k0[t]: {\"NoSolution\"} (a proof: n < 0 with c != 0, or gcd(a, b) not dividing c), {\"Zero\", alpha, beta} (q = beta solves it), or {\"Reduced\", b', c', m, alpha, beta}: every solution is q = alpha h + beta with D h + b' h = c', deg h <= m.";
brDTSPDE[C_Association, a0_, b0_, c0_, n0_] := Module[{C1 = brDTLevel[C, 1], a = a0, b = b0, c = c0, n = n0, alpha = 1, beta = 0, g, dv, eg, r, z, guard = 0},
  While[True,
    If[++guard > 100, throwInternal["brDTSPDENoProgress"]];
    c = brBCNorm[C1, c];
    If[c === 0, Return[{"Zero", alpha, beta}, Module]];
    If[n < 0, Return[{"NoSolution"}, Module]];
    g = brDTGCD[C1, a, b];
    dv = brDTDivide[C1, c, g];
    If[dv[[2]] =!= 0, Return[{"NoSolution"}, Module]];
    c = dv[[1]]; a = First[brDTDivide[C1, a, g]]; b = First[brDTDivide[C1, b, g]];
    If[brDTDeg[C1, a] === 0,
      Return[{"Reduced", brBCNorm[C1, b/a], brBCNorm[C1, c/a], n, alpha, beta}, Module]];
    (* b r + a z = c, deg r < deg a *)
    eg = brK0PolyExtGCD[C1["F"], C1["B"], brBCList[C1, b], brBCList[C1, a], C1["th"]];
    If[brK0Degree[eg[[1]], C1["th"]] =!= 0, throwInternal["brDTSPDENotCoprime"]];
    r = Last[brDTDivide[C1, brBCNorm[C1, brBCExpr[C1, eg[[2]]] c], a]];
    dv = brDTDivide[C1, brBCNorm[C1, c - b r], a];
    If[dv[[2]] =!= 0, throwInternal["brDTSPDEDiophantine"]];
    z = dv[[1]];
    b = brBCNorm[C1, b + brBCD[C1, a]];
    c = brBCNorm[C1, z - brBCD[C1, r]];
    n = n - brDTDeg[C1, a];
    beta = brBCNorm[C1, beta + alpha r];
    alpha = brBCNorm[C1, alpha a]]];

brDTGCD::usage = "brDTGCD[C1, a, b] is the monic gcd of a and b in k0[t] (brK0PolyGCD over the level-1 coefficient field), as a normal-form expression.";
brDTGCD[C1_Association, a_, b_] := brBCExpr[C1, brK0PolyGCD[C1["F"], C1["B"], brBCList[C1, a], brBCList[C1, b], C1["th"]]];

brDTPolyDeriv::usage = "brDTPolyDeriv[C, c, n] decides D q = c with q in k0[t], deg q <= n ([SI] IntegratePrimitivePolynomial / the Laurent coefficients, with an empty logarithmic part): primitive t top down through LI_0 (brBCLimited, every step's c forced), exponential t coefficientwise through RDE_0(i Dw, c_i) and brDTIsDerivK0 at t^0. {\"Solution\", q} (checked exactly) or {\"NoSolution\"} (a proof).";
brDTPolyDeriv[C_Association, c0_, n_] := Module[{C1 = brDTLevel[C, 1], t, c, d, eta, beta = 0, alg = 0, vc, s, isd},
  t = C1["th"]; c = brBCNorm[C1, c0];
  If[c === 0, Return[{"Solution", 0}]];
  d = brDTDeg[C1, c];
  If[d > n, Return[{"NoSolution"}]];
  If[C1["kind"] === "Log",
    eta = C1["Dth"];
    Do[
      vc = brCatch[brBCLimited[C1, brBCNorm[C1, brDTCoeff[C1, c, i] - (i + 1) beta eta]]];
      If[brFailQ[vc], If[First[brFailTag[vc]] === "NonElementary", Return[{"NoSolution"}, Module], Throw[vc, "TragerFail"]]];
      alg = alg + brBCNorm[C1, beta + vc[[2]]/(i + 1)] t^(i + 1);
      beta = vc[[1]],
      {i, d, 0, -1}];
    alg = brBCNorm[C1, alg + beta];
    If[brDTDeg[C1, alg] > n, Return[{"NoSolution"}]],
    Do[
      If[i === 0,
        isd = brDTIsDerivK0[C, brDTCoeff[C1, c, 0]];
        If[First[isd] === "no", Return[{"NoSolution"}, Module]];
        alg = alg + isd[[2]],
        s = brDTRDE0[C, brBCNorm[C1, i C1["Dh"]], brDTCoeff[C1, c, i]];
        If[First[s] === "NoSolution", Return[{"NoSolution"}, Module]];
        alg = alg + s[[2]] t^i],
      {i, 0, d}];
    alg = brBCNorm[C1, alg]];
  If[brBCNorm[C1, brBCD[C1, alg] - c] =!= 0, throwInternal["brDTPolyDerivCheck", <|"c" -> c0|>]];
  {"Solution", alg}];

brDTPolyRDE::usage = "brDTPolyRDE[C, b, c, n] solves D q + b q = c, deg q <= n, in k0[t] after SPDE ([SI] Sec. 6.5-6.6): b = 0 is brDTPolyDeriv; deg b > 0 the no-cancellation case (top down, each leading coefficient forced); b in k0* the cancellation case, top down through RDE_0(b + m Dw, c_m) (exponential) or RDE_0(b, c_m) (primitive) once LD_0(-b) resp. PLD_0(-b; Dw) shows the homogeneous equation has no polynomial solution (a homogeneous solution declines: BaseTower:cancellation-homogeneous). {\"Solution\", q} or {\"NoSolution\"} (a proof).";
brDTPolyRDE[C_Association, b0_, c0_, n0_] := Module[{C1 = brDTLevel[C, 1], t, b, c, n = n0, db, q = 0, m, p, s, stm, ld, guard = 0},
  t = C1["th"]; b = brBCNorm[C1, b0]; c = brBCNorm[C1, c0];
  If[c === 0, Return[{"Solution", 0}]];
  If[b === 0, Return[brDTPolyDeriv[C, c, n]]];
  db = brDTDeg[C1, b];
  If[db >= 1,
    While[c =!= 0,
      If[++guard > 500, throwInternal["brDTPolyRDENoProgress"]];
      m = brDTDeg[C1, c] - db;
      If[!(0 <= m <= n), Return[{"NoSolution"}, Module]];
      p = brBCNorm[C1, brDTCoeff[C1, c, brDTDeg[C1, c]]/brDTCoeff[C1, b, db]] t^m;
      q = q + p; n = m - 1;
      c = brBCNorm[C1, c - brBCD[C1, p] - b p]];
    Return[{"Solution", brBCNorm[C1, q]}]];
  (* b in k0*: the cancellation case *)
  ld = If[C1["kind"] === "Log", brDTPLD0[C, -b, 0], brDTPLD0[C, -b, C1["Dh"]]];
  If[First[ld] === "radical" && (C1["kind"] === "Log" || IntegerQ[ld[[3]]/ld[[2]]]),
    throwUnimplemented["BaseTower:cancellation-homogeneous", "the cancellation case of the Risch DE over k1 with a homogeneous solution (or one that needs ROOT_0)"]];
  While[c =!= 0,
    If[++guard > 500, throwInternal["brDTPolyRDENoProgress"]];
    m = brDTDeg[C1, c];
    If[m > n, Return[{"NoSolution"}, Module]];
    s = brDTRDE0[C, If[C1["kind"] === "Log", b, brBCNorm[C1, b + m C1["Dh"]]], brDTCoeff[C1, c, m]];
    If[First[s] === "NoSolution", Return[{"NoSolution"}, Module]];
    stm = s[[2]] t^m;
    q = q + stm; n = m - 1;
    c = brBCNorm[C1, c - b stm - brBCD[C1, stm]]];
  {"Solution", brBCNorm[C1, q]}];

brDTRDE::usage = "brDTRDE[C, f, g] is RDE_1: D q + f q = g over k1 = k0(t1) ([SI] Ch. 6, E1.1), for f = i D w with w in k1 (the exponential polynomial part of level 2): such an f is weakly normalised (its normal poles have order >= 2), so q has a normal pole only where g has one, of order at most its order there minus 1 (the bound h); the special denominator (exponential t1) by the order at t1 = 0 and PLD_0; the degree bound brDTDegreeBound; SPDE; brDTPolyRDE. {\"Solution\", q} (checked exactly in normal form) or {\"NoSolution\"} (a proof: valid bounds and decided auxiliaries, E2.4).";
brDTRDE[C_Association, f0_, g0_] := Module[{C1 = brDTLevel[C, 1], t, eta, f, g, E1, hN, ft, gt, M, a, b, c, ns = 0, nb, nc, N, alpha, pld, n, sp, hs, qb, q, Lb},
  t = C1["th"]; eta = If[C1["kind"] === "Exp", C1["Dh"], C1["Dth"]];
  f = brBCNorm[C1, f0]; g = brBCNorm[C1, g0];
  If[g === 0, Return[{"Solution", 0}]];
  (* normal denominator: -nu_pi(q) <= max(0, -nu_pi(g) - 1) at every normal pi *)
  E1 = brDTNormalDen[C1, g];
  hN = Expand[Times @@ (#[[1]]^(#[[2]] - 1) & /@ Select[tragerFactorList[E1], !FreeQ[#[[1]], t] &])];
  ft = brBCNorm[C1, f - brBCD[C1, hN]/hN]; gt = brBCNorm[C1, g hN];
  M = Expand[PolynomialLCM[brDTNormalDen[C1, ft], brDTNormalDen[C1, gt]]];
  a = M; b = brBCNorm[C1, M ft]; c = brBCNorm[C1, M gt];
  If[C1["kind"] === "Exp",
    (* special denominator: q = t^ns qbar, the order at t = 0 bounded below *)
    Lb = brBCLaurent[C1, b];
    nb = If[Lb === {}, Infinity, Min[First /@ Lb]];
    nc = Min[First /@ brBCLaurent[C1, c]];
    ns = Min[0, nc - Min[0, nb]];
    If[nb === 0,
      alpha = brBCNorm[C1, -(Lookup[Association[Rule @@@ Lb], 0, 0])/(a /. t -> 0)];
      pld = brDTPLD0[C, alpha, eta];
      If[First[pld] === "radical" && IntegerQ[pld[[3]]/pld[[2]]], ns = Min[ns, pld[[3]]/pld[[2]]]]];
    N = Max[0, -nb, ns - nc];
    b = brBCNorm[C1, (b + ns a eta) t^N];
    c = brBCNorm[C1, c t^(N - ns)];
    a = Expand[a t^N]];
  n = brDTDegreeBound[C, a, b, c];
  sp = brDTSPDE[C, a, b, c, n];
  Switch[First[sp],
    "NoSolution", Return[{"NoSolution"}],
    "Zero", qb = sp[[3]],
    _, hs = brDTPolyRDE[C, sp[[2]], sp[[3]], sp[[4]]];
       If[First[hs] === "NoSolution", Return[{"NoSolution"}]];
       qb = brBCNorm[C1, sp[[5]] hs[[2]] + sp[[6]]]];
  q = brBCNorm[C1, qb t^ns/hN];
  If[brBCNorm[C1, brBCD[C1, q] + f q - g] =!= 0, throwInternal["brDTRDECheck", <|"f" -> f0, "g" -> g0|>]];
  {"Solution", q}];

(* -------------------------------------------------------------------------- *)
(* Entry point.                                                                *)
(* -------------------------------------------------------------------------- *)

brShapeBTower::usage = "brShapeBTower[T, f, x] integrates the integrand of the tower record T (shape \"BTower\": two monomials over the curve, Phase 5 M1) by brDTInt at the top level. It returns <|\"F\" -> antiderivative (surface form), \"verified\" -> False|> (the pipeline verifies it) or throws: NonElementary only on a certificate of Thm E2cert.";
brShapeBTower[T_Association, f_, x_Symbol] := Module[{C, r},
  C = brDTContext[T];
  r = brDTInt[C, Length[C["syms"]], C["f"]];
  <|"F" -> brDTSurface[C, r["alg"] + Total[#[[1]] Log[#[[2]]] & /@ r["logs"]] + Total[r["rs"]]] + r["surf"],
    "verified" -> False|>];
