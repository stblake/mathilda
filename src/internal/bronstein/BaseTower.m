(* ::Package:: *)

(* ========================================================================== *)
(* BaseTower.m -- Shape B over the curve with ANY NUMBER of monomials (Phase   *)
(* 5, M1 two monomials, M2 the general tower): K = k0(t1)..(ts), k0 =         *)
(* Q(x)(Y), Y^m = G(x), each t_i = Log[u_i] or E^(w_i) with u_i, w_i in        *)
(* k_(i-1) (arguments may involve Y).  The specification is                   *)
(* docs/notes/bronstein_extensions.tex, sections E2 and E1, M1 and M2.        *)
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
(*     RDE_0 = brAlgRDE, RDE_i = brDTRDE ([SI] Ch. 6 over k_(i-1), Thm E1rde, *)
(*     with the weak normaliser of [SI] 6.1);                                 *)
(*   logarithmic t_j: LI_(j-1)(p_i; Dt_j) top down ([SI] 5.8.1),             *)
(*     LI_0 = brBCLimited (E1.3), LI_i = brDTLimitedAt (integrate over k_i,   *)
(*     then the structure theorem over k_i: Cor 9.3.1);                       *)
(* and the t_j^0 coefficient recurses to level j - 1, ending in Trager.       *)
(* Every routine takes the level as an argument; the oracles of a level i    *)
(* (brDTLI, brDTIsDeriv, brDTPLD, brDTRDEAt) dispatch to the level-0 ones at *)
(* i = 0.                                                                     *)
(*                                                                            *)
(* LD_i / PLD_i (Prop E1ld) are brDTPLD: INT_i, its Liouville data (rational  *)
(* coefficients, logands in k_i; else Decline: I2 at level 0, OL3 above),    *)
(* then the span over G_i and the logands.  ROOT_0 (OL5) is not built: a     *)
(* witness with n0 > 1 declines.  The cancellation case with an n = 1        *)
(* witness uses the homogeneous solution ([SI] 6.6, brDTCancelHom).  Every   *)
(* "no" is one of [SI] Thms 5.6.1, 5.8.1, 5.9.1, Ch. 6, a Trager certificate,*)
(* or a decided auxiliary (E2.4); every positive auxiliary answer is checked *)
(* exactly in normal form before it is used, and the final antiderivative is *)
(* verified by the pipeline.                                                 *)
(*                                                                            *)
(* Elements are expressions in x, Y and the monomials, in BaseCurve.m's       *)
(* normal form Sum c_i Y^i (brBCNorm); a polynomial in t over k_(j-1) is such *)
(* an expression polynomial in t.                                             *)
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

brDTAdmit::usage = "brDTAdmit[x, monos, kind, arg, k, cv] admits the kernel k = Log[arg] or E^arg (kind \"Log\" | \"Exp\", arg in x, the symbols of monos and powers of the curve's radicand) over the tower monos on the curve cv (E2.2). It returns {\"new\", monos'} (monos' = monos with the new monomial, flagged \"curve\"; brDFSpan proved Du/u resp. Dw outside span_Q G), or {\"rewrite\", e} with k = e an identity on the sampled regions: a dependent logarithm with kappa = 1 exactly and zero offset at every sample point (e = Sum r_g theta_g), a dependent exponential with c = 0 and integer coefficients (e = Prod u_g^r_g Prod t_g^r_g). Declines: Tower:dependent-log-offset (kappa != 1), Tower:dependent-log-branch (a nonzero 2 Pi I m/d), Tower:dependent-exp-offset (c != 0), Tower:dependent-exp-fractional (a non-integer coefficient), Tower:span-undecided. Any number of monomials is admitted (Phase 5 M2).";
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

brDTInt::usage = "brDTInt[C, j, e] integrates e in k_j (Thm E2level): Hermite and Rothstein-Trager in t_j over k_(j-1) (BaseCurve.m on brDTLevel[C, j]), the polynomial part (exponential t_j: RDE_(j-1)(i D w_j, p_i) for i != 0, brDTRDEAt; logarithmic t_j: LI_(j-1)(p_i; D t_j) top down, brDTLI) and the t_j^0 coefficient by brDTInt at level j - 1; level 0 is Trager (brDTIntegrateK0). It returns <|\"alg\" -> v in k_j, \"logs\" -> {{c, u}, ...} (u in k_j), \"rs\" -> {RootSum terms}, \"surf\" -> the level-0 integral (surface form)|>, the antiderivative being v + Sum c Log[u] + Sum rs + surf. Throws NonElementary only on a certificate (E2.4 (a)-(e)).";
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
        q = brDTRDEAt[C, j - 1, brBCNorm[Cj, l[[1]] Cj["Dh"]], l[[2]]];
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
      vc = brDTLI[C, j, brBCNorm[Cj, Lookup[es, i, 0] - (i + 1) beta eta]];
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

brDTRadicalSolve::usage = "brDTRadicalSolve[C, a, gens, vars] finds the least n >= 1 and an INTEGER vector e with n a = Sum e_k gens_k (a, gens in normal form over the level record C, rational in vars): {n, e}, or None when a is not in span_Q gens (a proof, Prop E2span), or $Failed when a span system was not decided. The gens may be dependent (Trager's logands beside the structure units): a greedy independent subset I gives the unique coordinates a_I and the relations q_d of the others; the solutions are a_I - Sum nu_d q_d (dependent coordinates nu_d), so the least n is the order of a_I in Q^I modulo the lattice Z^I + Sum Z q_d (Hermite normal form of its generators), and the integer combination is read off the HNF transformation. A minimal n matters: a consumer that needs an n = 1 witness declines (ROOT_0) on any larger one.";
brDTRadicalSolve[C_Association, a_, gens0_List, vars_List] := Module[{gens = gens0, K, ind = {}, rel = <||>, sp, aI, dl, L, M, U, H, B, c, n, w, e},
  K = Length[gens];
  Do[
    If[gens[[k]] =!= 0,
      sp = If[ind === {}, False, Risch`RationalSpan[gens[[k]], gens[[ind]], vars]];
      Which[
        sp === False, AppendTo[ind, k],
        ListQ[sp] && Length[sp] === Length[ind] && AllTrue[sp, rationalNumberQ] && brBCNorm[C, gens[[k]] - sp . gens[[ind]]] === 0,
          rel[k] = sp,
        True, Return[$Failed, Module]]],
    {k, K}];
  If[ind === {}, Return[If[brBCNorm[C, a] === 0, {1, ConstantArray[0, K]}, None]]];
  aI = Risch`RationalSpan[a, gens[[ind]], vars];
  If[aI === False, Return[None]];
  If[!ListQ[aI] || Length[aI] =!= Length[ind] || !AllTrue[aI, rationalNumberQ] || brBCNorm[C, a - aI . gens[[ind]]] =!= 0,
    Return[$Failed]];
  dl = Keys[rel];
  (* a relation found while the independent set was smaller: its later coordinates are 0 *)
  rel = PadRight[#, Length[ind]] & /@ rel;
  L = LCM @@ (Denominator /@ Flatten[Join[{aI, 1}, Values[rel]]]);
  M = Join[L IdentityMatrix[Length[ind]], L (rel[#] & /@ dl)];
  {U, H} = HermiteDecomposition[M];
  B = Take[H, Length[ind]];
  c = LinearSolve[Transpose[B], L aI];
  n = LCM @@ (Denominator /@ c);
  w = (n c) . Take[U, Length[ind]];
  e = ConstantArray[0, K];
  Do[e[[ind[[i]]]] = w[[i]], {i, Length[ind]}];
  Do[e[[dl[[d]]]] = w[[Length[ind] + d]], {d, Length[dl]}];
  If[!AllTrue[e, IntegerQ] || brBCNorm[C, n a - e . gens] =!= 0, throwInternal["brDTRadicalSolveCheck", <|"a" -> a|>]];
  {n, e}];

brDTPLD0::usage = "brDTPLD0[C, alpha, eta] decides PLD_0 (Prop E1ld, j = 0) for alpha, eta in k0 (eta = D w, or 0 for LD_0): {\"none\"} when there are no n >= 1, m in Z, U in k0* with n alpha = m eta + DU/U (a proof: Trager proved alpha non-integrable, or the span system with n = 1 is inconsistent), or {\"radical\", n0, m0, {{e_l, u_l}, ...}} with U = Prod u_l^e_l, checked exactly in normal form. Declines when Trager declines, on Liouville data outside I2 (brDTLiouville0) or an undecided span (BaseTower:pld-undecided).";
brDTPLD0[C_Association, alpha_, eta_] := Module[{C0 = brDTLevel[C, 0], a = brBCNorm[C, alpha], F0, v, logs, gens, sp, mq, rq, n0, ev},
  If[a === 0, Return[{"radical", 1, 0, {}}]];
  F0 = brCatch[brDTIntegrateK0[C, a]];
  If[brFailQ[F0], If[First[brFailTag[F0]] === "NonElementary", Return[{"none"}, Module], Throw[F0, "TragerFail"]]];
  {v, logs} = brDTLiouville0[C, F0];
  (* the Liouville data must reproduce alpha exactly before a "none" rests on it *)
  If[brBCNorm[C, a - brBCD[C0, v] - Sum[l[[1]] brBCD[C0, l[[2]]]/l[[2]], {l, logs}]] =!= 0,
    throwInternal["brDTLiouvilleCheck", <|"alpha" -> alpha|>]];
  gens = Join[If[eta === 0, {}, {brBCNorm[C, eta]}], brBCNorm[C, brBCD[C0, #[[2]]]/#[[2]]] & /@ logs];
  If[gens === {}, Return[{"none"}]];
  sp = brDTRadicalSolve[C0, a, gens, {C["x"], C["Y"]}];
  If[sp === None, Return[{"none"}]];
  If[sp === $Failed,
    throwUnimplemented["BaseTower:pld-undecided", "the span system of a (parametric) logarithmic-derivative problem over k0 was not decided"]];
  {n0, ev} = sp;
  mq = If[eta === 0, 0, First[ev]]; rq = If[eta === 0, ev, Rest[ev]];
  {"radical", n0, mq, Select[Transpose[{rq, #[[2]] & /@ logs}], #[[1]] =!= 0 &]}];

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
(* The oracles of level i >= 0 (Thm E1rde's table, Prop E1ld), each one       *)
(* dispatching to the level-0 oracle above at i = 0.                          *)
(* -------------------------------------------------------------------------- *)

brDTTheta::usage = "brDTTheta[C, g] is theta_g, the element of k_g whose derivative is the structure generator eta_g: the monomial t_g for a logarithm, its argument w_g (normal form) for an exponential E^(w_g).";
brDTTheta[C_Association, g_Integer] := If[C["kinds"][[g]] === "Log", C["syms"][[g]], C["argsY"][[g]]];

brDTUnit::usage = "brDTUnit[C, g] is the element u of k_g with D u/u = eta_g: the argument u_g (normal form) of a logarithm t_g = Log[u_g], the monomial t_g itself for an exponential.";
brDTUnit[C_Association, g_Integer] := If[C["kinds"][[g]] === "Log", C["argsY"][[g]], C["syms"][[g]]];

brDTLiouville::usage = "brDTLiouville[C, i, a] is INT_i(a) read as Liouville data over k_i: {\"none\"} when a has no elementary integral over k_i (a certificate of brDTInt), else {\"data\", V, {{c, u}, ...}} with a = D V + Sum c D u/u checked exactly in normal form, V and every u in k_i, every c RATIONAL and nonzero. An irrational c declines (BaseTower:irrational-log at level 0, gap I2; BaseTower:pld-irrational-above at a level >= 1, OL3); a RootSum or a logand outside k_i over Q declines (BaseTower:rootsum-logand); a Trager decline propagates.";
brDTLiouville[C_Association, i_Integer, a_] := Module[{Ci = brDTLevel[C, i], vars, r, v0, logs0, logs, V},
  vars = Join[{C["x"], C["Y"]}, Take[C["syms"], i]];
  r = brCatch[brDTInt[C, i, a]];
  If[brFailQ[r], If[First[brFailTag[r]] === "NonElementary", Return[{"none"}, Module], Throw[r, "TragerFail"]]];
  If[r["rs"] =!= {},
    throwUnimplemented["BaseTower:rootsum-logand", "a logarithmic part of an integral over the curve not given as rational multiples of logarithms of elements of the field (RootSum, ArcTan)"]];
  {v0, logs0} = brDTLiouville0[C, r["surf"]];
  logs = Join[r["logs"], logs0];
  If[!AllTrue[logs, rationalNumberQ[First[#]] &],
    If[i === 0,
      throwUnimplemented["BaseTower:irrational-log", "an irrational logarithm coefficient in the integral over k0 inside LD/PLD/LI (gap I2)"],
      throwUnimplemented["BaseTower:pld-irrational-above", "an irrational logarithm coefficient in the integral over k_j, j >= 1, inside LD/PLD/LI (OL3)"]]];
  logs = Select[{#[[1]], brBCNorm[Ci, #[[2]]]} & /@ logs, #[[1]] =!= 0 &];
  If[!AllTrue[logs, #[[2]] =!= 0 && brDTRationalQ[#[[2]], vars] &],
    throwUnimplemented["BaseTower:rootsum-logand", "a logarithmic part of an integral over the curve not given as rational multiples of logarithms of elements of the field (RootSum, ArcTan)"]];
  V = brBCNorm[Ci, r["alg"] + v0];
  (* the Liouville data must reproduce a exactly before a "no" rests on it *)
  If[brBCNorm[Ci, a - brBCD[Ci, V] - Sum[l[[1]] brBCD[Ci, l[[2]]]/l[[2]], {l, logs}]] =!= 0,
    throwInternal["brDTLiouvilleCheck", <|"a" -> a, "level" -> i|>]];
  {"data", V, logs}];

brDTLimitedAt::usage = "brDTLimitedAt[C, i, p, etas] solves p = D w + Sum c_k etas_k with w in k_i and RATIONAL c_k (the etas in k_i; etas = {} asks \"is p a derivative in k_i\"): {w, {c_k}}, checked exactly in normal form. INT_i gives p = D V + Sum c_l D u_l/u_l (brDTLiouville); a solution exists iff Sum c_l D u_l/u_l is in span_Q(etas, eta_1..eta_i) ([SI] Thm 9.3.1, Cor 9.3.1 over k_i, constants Q), decided by Risch`RationalSpan, and then w = V + Sum r_g theta_g. Throws NonElementary when there is none ([SI] Thm 5.8.1): p not elementary over k_i, or the span system inconsistent.";
brDTLimitedAt[C_Association, i_Integer, p_, etas_List] := Module[{Ci = brDTLevel[C, i], ne = Length[etas], ld, V, logs, theta, gens, sp, w, cs},
  If[brBCNorm[Ci, p] === 0, Return[{0, ConstantArray[0, ne]}]];
  ld = brDTLiouville[C, i, p];
  If[First[ld] === "none",
    throwNonElementary["no algebraic part: the limited integration problem has no solution over the tower on the curve (its integrand is not elementary over k" <> ToString[i] <> "; [SI] Thm 5.8.1)"]];
  {V, logs} = Rest[ld];
  theta = brBCNorm[Ci, Sum[l[[1]] brBCD[Ci, l[[2]]]/l[[2]], {l, logs}]];
  gens = brBCNorm[Ci, #] & /@ Join[etas, Take[C["etas"], i]];
  sp = Which[
    theta === 0, ConstantArray[0, Length[gens]],
    gens === {}, False,
    True, Risch`RationalSpan[theta, gens, Join[{C["x"], C["Y"]}, Take[C["syms"], i]]]];
  If[sp === False,
    throwNonElementary["no algebraic part: the limited integration problem has no solution over the tower on the curve (the logarithmic part of the integral over k" <> ToString[i] <> " is not in the span of the structure generators; [SI] Thm 5.8.1, Cor 9.3.1)"]];
  If[!ListQ[sp] || Length[sp] =!= Length[gens] || !AllTrue[sp, rationalNumberQ],
    throwUnimplemented["BaseTower:li-undecided", "the span system of a limited integration problem over k1 was not decided"]];
  cs = Take[sp, ne];
  w = brBCNorm[Ci, V + Sum[sp[[ne + g]] brDTTheta[C, g], {g, i}]];
  If[brBCNorm[Ci, p - brBCD[Ci, w] - cs . etas] =!= 0, throwInternal["brDTLimitedCheck", <|"p" -> p, "level" -> i|>]];
  {w, cs}];

brDTLI::usage = "brDTLI[C, j, p] is LI_(j-1)(p; D t_j) for a logarithmic t_j and p in k_(j-1): {w, c} with p = D w + c D t_j, w in k_(j-1), c rational. Level 1 is brBCLimited over k0 (E1.3); above it brDTLimitedAt over k_(j-1). Throws NonElementary when there is none ([SI] Thm 5.8.1).";
brDTLI[C_Association, 1, p_] := brBCLimited[brDTLevel[C, 1], p];
brDTLI[C_Association, j_Integer, p_] /; j >= 2 := With[{r = brDTLimitedAt[C, j - 1, p, {C["Ds"][[j]]}]}, {r[[1]], r[[2, 1]]}];

brDTLimited::usage = "brDTLimited[C, p] is LI_1(p; D t_2) (t_2 a logarithm over k1): brDTLI[C, 2, p].";
brDTLimited[C_Association, p_] := brDTLI[C, 2, p];

brDTIsDeriv::usage = "brDTIsDeriv[C, i, e] decides whether e in k_i is a derivative D v, v in k_i: {\"yes\", v} (checked exactly) or {\"no\"} (a proof). Level 0 is brDTIsDerivK0; above it the limited integration problem with no eta (brDTLimitedAt[C, i, e, {}]).";
brDTIsDeriv[C_Association, 0, e_] := brDTIsDerivK0[C, e];
brDTIsDeriv[C_Association, i_Integer, e_] /; i >= 1 := Module[{r = brCatch[brDTLimitedAt[C, i, e, {}]]},
  If[brFailQ[r], If[First[brFailTag[r]] === "NonElementary", {"no"}, Throw[r, "TragerFail"]], {"yes", r[[1]]}]];

brDTPLD::usage = "brDTPLD[C, i, alpha, eta] decides PLD_i (Prop E1ld) for alpha, eta in k_i (eta = D w with w in k_i, or 0 for LD_i): {\"none\"} (a proof: alpha has no elementary integral over k_i, or the span system n alpha - m eta in span_Q(G_i, D v_l/v_l) is inconsistent, Prop E1ld (i)-(ii)), or {\"radical\", n0, m0, {{e, u}, ...}} with n0 alpha = m0 eta + DU/U, U = Prod u^e in k_i (the logands of INT_i and the units of the structure generators: u_g of a logarithm, the monomial of an exponential), checked exactly. Level 0 is brDTPLD0. Declines on Liouville data with an irrational coefficient (OL3 above level 0) or an undecided span.";
brDTPLD[C_Association, 0, alpha_, eta_] := brDTPLD0[C, alpha, eta];
brDTPLD[C_Association, i_Integer, alpha_, eta_] /; i >= 1 := Module[
  {Ci = brDTLevel[C, i], a, ld, logs, gens, sp, ne, mq, re, rl, n0, facs, ev},
  a = brBCNorm[Ci, alpha];
  If[a === 0, Return[{"radical", 1, 0, {}}]];
  ld = brDTLiouville[C, i, a];
  If[First[ld] === "none", Return[{"none"}]];
  logs = ld[[3]];
  ne = If[eta === 0, 0, 1];
  gens = Join[If[eta === 0, {}, {brBCNorm[Ci, eta]}], brBCNorm[Ci, #] & /@ Take[C["etas"], i],
    brBCNorm[Ci, brBCD[Ci, #[[2]]]/#[[2]]] & /@ logs];
  sp = brDTRadicalSolve[Ci, a, gens, Join[{C["x"], C["Y"]}, Take[C["syms"], i]]];
  If[sp === None, Return[{"none"}]];
  If[sp === $Failed,
    throwUnimplemented["BaseTower:pld-undecided", "the span system of a (parametric) logarithmic-derivative problem over k0 was not decided"]];
  {n0, ev} = sp;
  mq = If[eta === 0, 0, First[ev]];
  re = Take[ev, {ne + 1, ne + i}]; rl = Drop[ev, ne + i];
  facs = Select[Join[Table[{re[[g]], brDTUnit[C, g]}, {g, i}], Table[{rl[[l]], logs[[l, 2]]}, {l, Length[logs]}]],
    #[[1]] =!= 0 &];
  If[brBCNorm[Ci, n0 a - mq If[eta === 0, 0, eta] - Sum[u[[1]] brBCD[Ci, u[[2]]]/u[[2]], {u, facs}]] =!= 0,
    throwInternal["brDTPLDCheck", <|"alpha" -> alpha, "level" -> i|>]];
  {"radical", n0, mq, facs}];

brDTRDEAt::usage = "brDTRDEAt[C, i, f, g] is RDE_i: D q + f q = g over k_i: brDTRDE0 at level 0, brDTRDE above it. {\"Solution\", q} or {\"NoSolution\"} (a proof).";
brDTRDEAt[C_Association, 0, f_, g_] := brDTRDE0[C, f, g];
brDTRDEAt[C_Association, i_Integer, f_, g_] /; i >= 1 := brDTRDE[C, i, f, g];

(* -------------------------------------------------------------------------- *)
(* RDE_j, j >= 1: [SI] Ch. 6 over k_(j-1), t = t_j primitive or exponential   *)
(* (Thm E1rde).  Polynomials in t over k_(j-1) are normal-form expressions     *)
(* polynomial in t; every call to k_(j-1) is one of the level oracles above.  *)
(* Each routine also has its level-1 form (the M1 signature, one argument    *)
(* fewer), kept for the unit suites.                                           *)
(* -------------------------------------------------------------------------- *)

brDTDeg::usage = "brDTDeg[Cj, e] is the degree in t = Cj[\"th\"] of e in k_(j-1)[t] (normal form), -Infinity for 0; an internal error when a coordinate denominator involves t.";
brDTDeg[C1_Association, e_] := Module[{t = C1["th"], n = brBCNorm[C1, e], cs},
  If[n === 0, Return[-Infinity]];
  cs = Together[Coefficient[n, C1["Y"], #]] & /@ Range[0, C1["m"] - 1];
  If[!FreeQ[Denominator /@ cs, t], throwInternal["brDTNotPolynomial", <|"e" -> e|>]];
  Max[If[# === 0, -Infinity, Exponent[Numerator[#], t]] & /@ cs]];

brDTCoeff::usage = "brDTCoeff[Cj, e, k] is the coefficient of t^k of e in k_(j-1)[t] (an element of k_(j-1) in normal form).";
brDTCoeff[C1_Association, e_, k_Integer] := Module[{t = C1["th"], Y = C1["Y"], n = brBCNorm[C1, e]},
  brBCNorm[C1, Sum[With[{c = Together[Coefficient[n, Y, i]]},
      Coefficient[Expand[Numerator[c]], t, k]/Denominator[c]] Y^i, {i, 0, C1["m"] - 1}]]];

brDTDivide::usage = "brDTDivide[Cj, a, b] is {q, r} with a = q b + r, deg r < deg b, in k_(j-1)[t] (long division, leading coefficients inverted in k_(j-1)).";
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

brDTNormalDen::usage = "brDTNormalDen[Cj, e] is the normal part of the t-denominator of e (a polynomial in x and the monomials over Q, free of Y): the t-dependent factors of its normal-form denominator, without the special factor t for an exponential t. Since 1, Y, ..., Y^(m-1) is a basis of k_(j-1)[t] over Q(x, t_1..t_(j-1))[t], an irreducible rho of Q(x, t_1..)[t] divides it with multiplicity e exactly when max over the places pi | rho of -nu_pi(e) is e (rho is squarefree over k_(j-1)).";
brDTNormalDen[C1_Association, e_] := Module[{t = C1["th"], d},
  d = Last[brBCSplitDen[C1, e]];
  Expand[Times @@ (#[[1]]^#[[2]] & /@ Select[tragerFactorList[d],
    !FreeQ[#[[1]], t] && !(C1["kind"] === "Exp" && NumberQ[Together[#[[1]]/t]]) &])]];

brDTWeakNormalizer::usage = "brDTWeakNormalizer[C, j, f] is a weak normaliser q in k_(j-1)[t_j] of f in k_j ([SI] Sec. 6.1, Thm 6.1.1): f - Dq/q has no simple normal pole whose residue is a positive integer. For d1 the product of the normal factors of multiplicity 1 in the denominator of f = N/d (d = d1 E), the residues at the roots of d1 are the roots of R(z) = Res_t(N - z E D d1, d1) in k_(j-1)[z]; q = Prod gcd(N - n E D d1, d1)^n over the positive integer roots n. 1 when there is none (always so for f = D w: it has no simple pole).";
brDTWeakNormalizer[C_Association, j_Integer, f_] := Module[
  {Cj = brDTLevel[C, j], t, z = Unique["Bronstein`Private`brWNz"], N, d, fl, d1, E2, Dd1, P, R, Rn, polys, gz, ns, q = 1, oth},
  t = Cj["th"];
  If[brBCNorm[Cj, f] === 0, Return[1]];
  {N, d} = brBCSplitDen[Cj, f];
  fl = Select[tragerFactorList[d], !FreeQ[#[[1]], t] && !(Cj["kind"] === "Exp" && NumberQ[Together[#[[1]]/t]]) &];
  d1 = Expand[Times @@ (First /@ Select[fl, #[[2]] === 1 &])];
  If[FreeQ[d1, t], Return[1]];
  E2 = Together[d/d1];
  If[!PolynomialQ[E2, t], throwInternal["brDTWeakNormalizerCofactor", <|"f" -> f|>]];
  Dd1 = brBCD[Cj, d1];
  P = Together[N - z E2 Dd1];
  If[!FreeQ[Denominator[P], t], throwInternal["brDTWeakNormalizerDen", <|"f" -> f|>]];
  R = Resultant[Expand[Numerator[P]], d1, t];
  Rn = brBCNorm[Cj, R];
  If[Rn === 0, throwInternal["brDTWeakNormalizerResultant", <|"f" -> f|>]];
  (* an integer n is a root iff every Q-coefficient (in x, Y, the monomials) of R vanishes at n *)
  polys = Flatten[Table[
    With[{num = Expand[Numerator[Together[Coefficient[Rn, Cj["Y"], i]]]]},
      If[num === 0, {}, (oth = Complement[Variables[num], {z}];
        If[oth === {}, {num}, Last /@ CoefficientRules[num, oth]])]],
    {i, 0, Cj["m"] - 1}]];
  polys = Select[Expand /@ polys, # =!= 0 &];
  If[polys === {}, throwInternal["brDTWeakNormalizerResultant", <|"f" -> f|>]];
  gz = PolynomialGCD @@ polys;
  ns = Select[Cases[FactorList[gz], {fac_, _} /; Exponent[fac, z] === 1 :> -Coefficient[fac, z, 0]/Coefficient[fac, z, 1]],
    IntegerQ[#] && # > 0 &];
  Do[q = q brDTGCD[Cj, brBCNorm[Cj, N - n E2 Dd1], d1]^n, {n, ns}];
  brBCNorm[Cj, q]];

brDTDegreeBound::usage = "brDTDegreeBound[C, j, a, b, c] is a valid bound n on deg q for a D q + b q = c in k_(j-1)[t_j] ([SI] Sec. 6.3, Thm E1rde): primitive t: max(0, dc - db) when db > da, else max(0, dc - da + 1), raised to c when db = da - 1 and LI_(j-1)(-lc b/lc a; Dt) has c a positive integer, or when db = da, LD_(j-1)(alpha) has a witness z (n0 = 1; n0 > 1 declines, ROOT_0) and LI_(j-1)(beta; Dt) has c a positive integer (beta = -(a_(da-1) Dz + b_(da-1) z)/(z lc a)); exponential t: max(0, dc - max(da, db)), raised to m0/n0 when da = db and PLD_(j-1)(-lc b/lc a; Dw) has a radical solution with m0/n0 an integer (remark (2) of E1.1). A non-integer c or ratio is discarded, never rounded. brDTDegreeBound[C, a, b, c] is level 1.";
brDTDegreeBound[C_Association, a_, b_, c_] := brDTDegreeBound[C, 1, a, b, c];
brDTDegreeBound[C_Association, j_Integer, a_, b_, c_] := Module[{C1 = brDTLevel[C, j], Cl = brDTLevel[C, j - 1], da, db, dc, n, alpha, vc, ld, z, beta, pld},
  da = brDTDeg[C1, a]; db = brDTDeg[C1, b]; dc = brDTDeg[C1, c];
  alpha = If[db === -Infinity, 0, brBCNorm[C1, -brDTCoeff[C1, b, db]/brDTCoeff[C1, a, da]]];
  If[C1["kind"] === "Log",
    n = If[db > da, Max[0, dc - db], Max[0, dc - da + 1]];
    If[db === da - 1,
      vc = brCatch[brDTLI[C, j, alpha]];
      If[brFailQ[vc], If[First[brFailTag[vc]] =!= "NonElementary", Throw[vc, "TragerFail"]],
        If[IntegerQ[vc[[2]]] && vc[[2]] > 0, n = Max[n, vc[[2]]]]]];
    If[db === da,
      ld = brDTPLD[C, j - 1, alpha, 0];
      If[First[ld] === "radical",
        If[ld[[2]] =!= 1, throwUnimplemented["BaseTower:root0", "a logarithmic-derivative witness that needs an n-th root in k0 (ROOT_0, OL5)"]];
        z = brBCNorm[C1, Product[u[[2]]^u[[1]], {u, ld[[4]]}]];
        beta = brBCNorm[C1, -(brDTCoeff[C1, a, da - 1] brBCD[Cl, z] + brDTCoeff[C1, b, da - 1] z)/(z brDTCoeff[C1, a, da])];
        vc = brCatch[brDTLI[C, j, beta]];
        If[brFailQ[vc], If[First[brFailTag[vc]] =!= "NonElementary", Throw[vc, "TragerFail"]],
          If[IntegerQ[vc[[2]]] && vc[[2]] > 0, n = Max[n, vc[[2]]]]]]],
    n = Max[0, dc - Max[db, da]];
    If[da === db,
      pld = brDTPLD[C, j - 1, alpha, C1["Dh"]];
      If[First[pld] === "radical" && IntegerQ[pld[[3]]/pld[[2]]], n = Max[n, pld[[3]]/pld[[2]]]]]];
  n];

brDTSPDE::usage = "brDTSPDE[C, j, a, b, c, n] is Rothstein's SPDE ([SI] Sec. 6.4) for a D q + b q = c, deg q <= n, in k_(j-1)[t_j]: {\"NoSolution\"} (a proof: n < 0 with c != 0, or gcd(a, b) not dividing c), {\"Zero\", alpha, beta} (q = beta solves it), or {\"Reduced\", b', c', m, alpha, beta}: every solution is q = alpha h + beta with D h + b' h = c', deg h <= m. brDTSPDE[C, a, b, c, n] is level 1.";
brDTSPDE[C_Association, a0_, b0_, c0_, n0_] := brDTSPDE[C, 1, a0, b0, c0, n0];
brDTSPDE[C_Association, j_Integer, a0_, b0_, c0_, n0_] := Module[{C1 = brDTLevel[C, j], a = a0, b = b0, c = c0, n = n0, alpha = 1, beta = 0, g, dv, eg, r, z, guard = 0},
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

brDTGCD::usage = "brDTGCD[Cj, a, b] is the monic gcd of a and b in k_(j-1)[t] (brK0PolyGCD over the level's coefficient field), as a normal-form expression.";
brDTGCD[C1_Association, a_, b_] := brBCExpr[C1, brK0PolyGCD[C1["F"], C1["B"], brBCList[C1, a], brBCList[C1, b], C1["th"]]];

brDTPolyDeriv::usage = "brDTPolyDeriv[C, j, c, n] decides D q = c with q in k_(j-1)[t_j], deg q <= n ([SI] IntegratePrimitivePolynomial / the Laurent coefficients, with an empty logarithmic part): primitive t top down through LI_(j-1) (brDTLI, every step's c forced), exponential t coefficientwise through RDE_(j-1)(i Dw, c_i) and brDTIsDeriv at t^0. {\"Solution\", q} (checked exactly) or {\"NoSolution\"} (a proof). brDTPolyDeriv[C, c, n] is level 1.";
brDTPolyDeriv[C_Association, c0_, n_] := brDTPolyDeriv[C, 1, c0, n];
brDTPolyDeriv[C_Association, j_Integer, c0_, n_] := Module[{C1 = brDTLevel[C, j], t, c, d, eta, beta = 0, alg = 0, vc, s, isd},
  t = C1["th"]; c = brBCNorm[C1, c0];
  If[c === 0, Return[{"Solution", 0}]];
  d = brDTDeg[C1, c];
  If[d > n, Return[{"NoSolution"}]];
  If[C1["kind"] === "Log",
    eta = C1["Dth"];
    Do[
      vc = brCatch[brDTLI[C, j, brBCNorm[C1, brDTCoeff[C1, c, i] - (i + 1) beta eta]]];
      If[brFailQ[vc], If[First[brFailTag[vc]] === "NonElementary", Return[{"NoSolution"}, Module], Throw[vc, "TragerFail"]]];
      alg = alg + brBCNorm[C1, beta + vc[[2]]/(i + 1)] t^(i + 1);
      beta = vc[[1]],
      {i, d, 0, -1}];
    alg = brBCNorm[C1, alg + beta];
    If[brDTDeg[C1, alg] > n, Return[{"NoSolution"}]],
    Do[
      If[i === 0,
        isd = brDTIsDeriv[C, j - 1, brDTCoeff[C1, c, 0]];
        If[First[isd] === "no", Return[{"NoSolution"}, Module]];
        alg = alg + isd[[2]],
        s = brDTRDEAt[C, j - 1, brBCNorm[C1, i C1["Dh"]], brDTCoeff[C1, c, i]];
        If[First[s] === "NoSolution", Return[{"NoSolution"}, Module]];
        alg = alg + s[[2]] t^i],
      {i, 0, d}];
    alg = brBCNorm[C1, alg]];
  If[brBCNorm[C1, brBCD[C1, alg] - c] =!= 0, throwInternal["brDTPolyDerivCheck", <|"c" -> c0|>]];
  {"Solution", alg}];

brDTLaurentDeriv::usage = "brDTLaurentDeriv[C, j, L] decides D p = Sum c_i t^i for an exponential t = t_j and the Laurent polynomial L = {{i, c_i}, ...} (c_i in k_(j-1)), p in k_(j-1)[t, 1/t]: D(p_i t^i) = (D p_i + i Dw p_i) t^i, so coefficientwise RDE_(j-1)(i Dw, c_i) for i != 0 and brDTIsDeriv at t^0. {\"Solution\", {{i, p_i}, ...}} (the t^0 coefficient is determined up to a constant) or {\"NoSolution\"} (a proof).";
brDTLaurentDeriv[C_Association, j_Integer, L_List] := Module[{C1 = brDTLevel[C, j], out = {}, s},
  If[C1["kind"] =!= "Exp", throwInternal["brDTLaurentDerivKind", <|"j" -> j|>]];
  Do[
    If[l[[1]] === 0,
      s = brDTIsDeriv[C, j - 1, l[[2]]];
      If[First[s] === "no", Return[{"NoSolution"}, Module]],
      s = brDTRDEAt[C, j - 1, brBCNorm[C1, l[[1]] C1["Dh"]], l[[2]]];
      If[First[s] === "NoSolution", Return[{"NoSolution"}, Module]]];
    If[s[[2]] =!= 0, AppendTo[out, {l[[1]], s[[2]]}]],
    {l, L}];
  {"Solution", out}];

brDTCancelHom::usage = "brDTCancelHom[C, j, b, c, n, U, m] is the cancellation case of D q + b q = c in k_(j-1)[t_j] (b in k_(j-1)*, c != 0) when the homogeneous equation has the solution qh = U t_j^m (-b = DU/U + m eta, an n = 1 witness of PLD_(j-1); m = 0 for a primitive t_j), verified exactly: q = qh r with D r = c/qh ([SI] PolyRischDECancelPrim, PolyRischDECancelExp; remark (3) of E1.1). Primitive t_j: D r = c/U by brDTPolyDeriv (deg r = deg q <= n), q = U r. Exponential t_j: D r = c U^-1 t^-m coefficientwise (brDTLaurentDeriv), r's t^0 coefficient free up to a constant, so q = U Sum r_i t^(i+m) is in k_(j-1)[t] with deg q <= n iff every exponent i + m lies in [0, n] once a constant t^0 coefficient is dropped when m is outside [0, n]. {\"Solution\", q} (checked exactly) or {\"NoSolution\"} (a proof).";
brDTCancelHom[C_Association, j_Integer, b_, c_, n_, U_, m_Integer] := Module[{C1 = brDTLevel[C, j], t, qh, r, L, p, p0, q},
  t = C1["th"];
  qh = brBCNorm[C1, U t^m];
  If[brBCNorm[C1, brBCD[C1, qh] + b qh] =!= 0, throwInternal["brDTCancelHomWitness", <|"b" -> b, "U" -> U, "m" -> m|>]];
  If[C1["kind"] === "Log",
    r = brDTPolyDeriv[C, j, brBCNorm[C1, c/U], n];
    If[First[r] === "NoSolution", Return[{"NoSolution"}]];
    q = brBCNorm[C1, U r[[2]]],
    L = Select[Table[{i - m, brBCNorm[C1, brDTCoeff[C1, c, i]/U]}, {i, 0, Max[0, brDTDeg[C1, c]]}], #[[2]] =!= 0 &];
    r = brDTLaurentDeriv[C, j, L];
    If[First[r] === "NoSolution", Return[{"NoSolution"}]];
    p = r[[2]];
    p0 = Lookup[Association[Rule @@@ p], 0, 0];
    If[!(0 <= m <= n) && p0 =!= 0,
      If[FreeQ[p0, Alternatives @@ Join[{C1["x"], C1["Y"]}, C["syms"]]],
        p = Select[p, #[[1]] =!= 0 &],
        Return[{"NoSolution"}]]];
    If[AnyTrue[p, !(0 <= #[[1]] + m <= n) &], Return[{"NoSolution"}]];
    q = brBCNorm[C1, U Sum[pp[[2]] t^(pp[[1]] + m), {pp, p}]]];
  If[brBCNorm[C1, brBCD[C1, q] + b q - c] =!= 0, throwInternal["brDTCancelHomCheck", <|"b" -> b, "c" -> c|>]];
  {"Solution", q}];

brDTPolyRDE::usage = "brDTPolyRDE[C, j, b, c, n] solves D q + b q = c, deg q <= n, in k_(j-1)[t_j] after SPDE ([SI] Sec. 6.5-6.6): b = 0 is brDTPolyDeriv; deg b > 0 the no-cancellation case (top down, each leading coefficient forced); b in k_(j-1)* the cancellation case: when LD_(j-1)(-b) resp. PLD_(j-1)(-b; Dw) has an n = 1 witness the homogeneous solution is used (brDTCancelHom), a witness with n0 > 1 declines (BaseTower:root0, ROOT_0), otherwise top down through RDE_(j-1)(b + m Dw, c_m) (exponential) or RDE_(j-1)(b, c_m) (primitive). {\"Solution\", q} or {\"NoSolution\"} (a proof). brDTPolyRDE[C, b, c, n] is level 1.";
brDTPolyRDE[C_Association, b0_, c0_, n0_] := brDTPolyRDE[C, 1, b0, c0, n0];
brDTPolyRDE[C_Association, j_Integer, b0_, c0_, n0_] := Module[{C1 = brDTLevel[C, j], t, b, c, n = n0, db, q = 0, m, p, s, stm, ld, guard = 0},
  t = C1["th"]; b = brBCNorm[C1, b0]; c = brBCNorm[C1, c0];
  If[c === 0, Return[{"Solution", 0}]];
  If[b === 0, Return[brDTPolyDeriv[C, j, c, n]]];
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
  (* b in k_(j-1)*: the cancellation case *)
  ld = If[C1["kind"] === "Log", brDTPLD[C, j - 1, -b, 0], brDTPLD[C, j - 1, -b, C1["Dh"]]];
  If[First[ld] === "radical" && (C1["kind"] === "Log" || IntegerQ[ld[[3]]/ld[[2]]]),
    If[ld[[2]] =!= 1,
      throwUnimplemented["BaseTower:root0", "a logarithmic-derivative witness that needs an n-th root in k0 (ROOT_0, OL5)"]];
    Return[brDTCancelHom[C, j, b, c, n, brBCNorm[C1, Product[u[[2]]^u[[1]], {u, ld[[4]]}]], ld[[3]]]]];
  While[c =!= 0,
    If[++guard > 500, throwInternal["brDTPolyRDENoProgress"]];
    m = brDTDeg[C1, c];
    If[m > n, Return[{"NoSolution"}, Module]];
    s = brDTRDEAt[C, j - 1, If[C1["kind"] === "Log", b, brBCNorm[C1, b + m C1["Dh"]]], brDTCoeff[C1, c, m]];
    If[First[s] === "NoSolution", Return[{"NoSolution"}, Module]];
    stm = s[[2]] t^m;
    q = q + stm; n = m - 1;
    c = brBCNorm[C1, c - b stm - brBCD[C1, stm]]];
  {"Solution", brBCNorm[C1, q]}];

brDTRDE::usage = "brDTRDE[C, j, f, g] is RDE_j: D q + f q = g over k_j = k_(j-1)(t_j) ([SI] Ch. 6, Thm E1rde) for any f, g in k_j: the weak normaliser w (brDTWeakNormalizer; q = r/w with r solving the RDE of f - Dw/w and w g), then brDTRDEWeak. {\"Solution\", q} (checked exactly in normal form) or {\"NoSolution\"} (a proof: valid bounds and decided auxiliaries, E2.4). brDTRDE[C, f, g] is level 1.";
brDTRDE[C_Association, f0_, g0_] := brDTRDE[C, 1, f0, g0];
brDTRDE[C_Association, j_Integer, f0_, g0_] := Module[{C1 = brDTLevel[C, j], f, g, w, r, q},
  f = brBCNorm[C1, f0]; g = brBCNorm[C1, g0];
  If[g === 0, Return[{"Solution", 0}]];
  w = brDTWeakNormalizer[C, j, f];
  If[w === 1, Return[brDTRDEWeak[C, j, f, g]]];
  r = brDTRDEWeak[C, j, brBCNorm[C1, f - brBCD[C1, w]/w], brBCNorm[C1, w g]];
  If[First[r] === "NoSolution", Return[r]];
  q = brBCNorm[C1, r[[2]]/w];
  If[brBCNorm[C1, brBCD[C1, q] + f q - g] =!= 0, throwInternal["brDTRDECheck", <|"f" -> f0, "g" -> g0|>]];
  {"Solution", q}];

brDTRDEWeak::usage = "brDTRDEWeak[C, j, f, g] is RDE_j for a WEAKLY NORMALISED f ([SI] Def 6.1.1: no simple normal pole with a positive integer residue; f = i D w always is): q has a normal pole only where g has one, of order at most its order there minus 1 (the bound h, [SI] Thm 6.1.2); the special denominator (exponential t_j) by the order at t_j = 0 and PLD_(j-1); the degree bound brDTDegreeBound; SPDE; brDTPolyRDE. {\"Solution\", q} (checked exactly) or {\"NoSolution\"} (a proof).";
brDTRDEWeak[C_Association, j_Integer, f0_, g0_] := Module[{C1 = brDTLevel[C, j], t, eta, f, g, E1, hN, ft, gt, M, a, b, c, ns = 0, nb, nc, N, alpha, pld, n, sp, hs, qb, q, Lb},
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
      pld = brDTPLD[C, j - 1, alpha, eta];
      If[First[pld] === "radical" && IntegerQ[pld[[3]]/pld[[2]]], ns = Min[ns, pld[[3]]/pld[[2]]]]];
    N = Max[0, -nb, ns - nc];
    b = brBCNorm[C1, (b + ns a eta) t^N];
    c = brBCNorm[C1, c t^(N - ns)];
    a = Expand[a t^N]];
  n = brDTDegreeBound[C, j, a, b, c];
  sp = brDTSPDE[C, j, a, b, c, n];
  Switch[First[sp],
    "NoSolution", Return[{"NoSolution"}],
    "Zero", qb = sp[[3]],
    _, hs = brDTPolyRDE[C, j, sp[[2]], sp[[3]], sp[[4]]];
       If[First[hs] === "NoSolution", Return[{"NoSolution"}]];
       qb = brBCNorm[C1, sp[[5]] hs[[2]] + sp[[6]]]];
  q = brBCNorm[C1, qb t^ns/hN];
  If[brBCNorm[C1, brBCD[C1, q] + f q - g] =!= 0, throwInternal["brDTRDECheck", <|"f" -> f0, "g" -> g0|>]];
  {"Solution", q}];

(* -------------------------------------------------------------------------- *)
(* Entry point.                                                                *)
(* -------------------------------------------------------------------------- *)

brShapeBTower::usage = "brShapeBTower[T, f, x] integrates the integrand of the tower record T (shape \"BTower\": two or more monomials over the curve, Phase 5 M1/M2) by brDTInt at the top level. It returns <|\"F\" -> antiderivative (surface form), \"verified\" -> False|> (the pipeline verifies it) or throws: NonElementary only on a certificate of Thm E2cert.";
brShapeBTower[T_Association, f_, x_Symbol] := Module[{C, r},
  C = brDTContext[T];
  r = brDTInt[C, Length[C["syms"]], C["f"]];
  <|"F" -> brDTSurface[C, r["alg"] + Total[#[[1]] Log[#[[2]]] & /@ r["logs"]] + Total[r["rs"]]] + r["surf"],
    "verified" -> False|>];
