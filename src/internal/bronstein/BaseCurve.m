(* ::Package:: *)

(* ========================================================================== *)
(* BaseCurve.m -- Shape B, the curve-monomial case: K = Q(x)(Y)(theta), with  *)
(* Y^m = G(x) and theta = Log[h] or E^h for ONE h in Q(x)(Y) that involves Y   *)
(* (BRONSTEIN_PLAN.md "Shape B"; [B] 1990 Parts 2-3 with n = 1 over the       *)
(* algebraic base k0 = Q(x)(Y); ISSAC'91 for the Risch DE on the curve).      *)
(*                                                                            *)
(* D theta = h'/h or theta h' lies in k0 and mixes the weights, so the        *)
(* integrand is a general element of k0(theta) and the algorithm is the       *)
(* transcendental one over the field k0, every coefficient operation done in  *)
(* k0 (Y^m = G):                                                              *)
(*   Hermite in theta over k0[theta] (extended Euclid in k0[theta], K0.m);    *)
(*   Rothstein-Trager: R(z) = Res_theta(a - z D d, d), formed with Y an       *)
(*     indeterminate and reduced modulo Y^m - G (the resultant is a           *)
(*     polynomial in the coefficients), made monic in k0: its coefficients    *)
(*     must be constants ([B] Prop 2.9 over k0, constants Q: the curve is     *)
(*     reduced and absolutely irreducible); logands are gcds in k0[theta];    *)
(*   polynomial part, exponential theta: b_j' + j h' b_j = e_j over Q(x)(Y)   *)
(*     for j != 0 (RDEAlgebraic.m, Part 3), the theta^0 coefficient           *)
(*     integrated in k0 (Trager);                                             *)
(*   polynomial part, logarithmic theta: top down, the limited integration    *)
(*     problem r = D w + c eta (eta = h'/h) in k0: integrate r by Trager and  *)
(*     decide the constant from r - D v = c eta exactly (a nonzero constant   *)
(*     combination of logarithmic derivatives that is an exact derivative in  *)
(*     k0 has a zero residue divisor, so it vanishes: the test decides).      *)
(*                                                                            *)
(* Elements are expressions in x, the symbol Y and theta, normalised to        *)
(* Sum c_i(x, theta) Y^i, i < m (brBCNorm), and coordinate lists on Trager's  *)
(* basis over fieldFn[{theta}] for the k0 arithmetic.  Every certificate is   *)
(* one of the theorems above; the caller verifies every answer.               *)
(* ========================================================================== *)

brBCContext::usage = "brBCContext[T] is the curve-monomial context of the tower record T (shape \"BCurve\"): <|\"x\", \"Y\", \"m\", \"G\", \"B\", \"F\", \"th\", \"kind\", \"h\", \"Dh\", \"Dth\", \"sY\", \"Ysurf\", \"f\"|> with the integrand \"f\" and the monomial argument \"h\" in x, Y and theta, Y^m = G the reduced curve of the radical (componentCurve of yh), Ysurf the surface expression of Y. Throws a Decline when the radical's curve is not reduced or not exact.";
brBCContext[T_Association] := Module[
  {x = T["x"], mono = First[T["monos"]], th = First[T["syms"]], cf, n, g, cc, Y = Unique["Bronstein`Private`brYc"],
   yh, toY, f, h, sY, C, Dh, rad = T["radical"], q, s},
  cf = brCanonicalForm[T];
  n = cf["n"]; g = cf["g"]; yh = cf["y"];
  cc = componentCurve[n, g, 1, x];
  If[!cc["exact"] || cc["m"] =!= n || !reducedCurveQ[cc["m"], cc["G"], x],
    throwUnimplemented["BaseCurve:curve-not-reduced", "curve-monomial case: the radical's curve is not reduced over Q(x) (a constant-field extension)"]];
  s = cc["s"]; q = cf["q"];
  (* yh = s Y; the user's b^(1/n) = yh/q = s Y/q *)
  f = Total[Table[cf["coeffs"][[i + 1]] (s Y)^i, {i, 0, n - 1}]];
  toY[e_] := brSubstituteRadical[e, {x}, rad["b"], n, Y] /. Y -> s Y/q;
  h = toY[mono["arg"]];
  sY = Together[D[cc["G"], x]/(n cc["G"])];
  C = <|"x" -> x, "Y" -> Y, "m" -> n, "G" -> Expand[cc["G"]], "B" -> curveBasis[n, Expand[cc["G"]], x],
        "F" -> fieldFn[{th}], "th" -> th, "kind" -> mono["kind"], "sY" -> sY,
        "Ysurf" -> Together[cf["yexpr"]/s], "kernel" -> mono["kernel"],
        "b" -> rad["b"], "phi" -> Together[s/q]|>;
  C = Append[C, "h" -> brBCNorm[C, h]];
  Dh = brBCNorm[C, D[C["h"], x] + D[C["h"], Y] sY Y];
  C = Join[C, <|"Dh" -> Dh, "Dth" -> brBCNorm[C, If[C["kind"] === "Log", Dh/C["h"], th Dh]], "f" -> brBCNorm[C, f]|>];
  C];

brBCNorm::usage = "brBCNorm[C, e] is the expression e (rational in x, Y, theta and any parameter) in the normal form Sum c_i Y^i, i < m, c_i rational in the other symbols (denominator rationalised modulo Y^m - G).";
brBCNorm[C_Association, e_] := Module[{Y = C["Y"], r},
  r = brAlgReduceY[e, Y, C["m"], C["G"]];
  Sum[Together[Coefficient[Expand[Numerator[Together[r]]], Y, i]/Denominator[Together[r]]] Y^i, {i, 0, C["m"] - 1}]];

brBCList::usage = "brBCList[C, e] is the element e (in x, Y, theta) as a coordinate list on the basis C[\"B\"] over the field F (default C[\"F\"]).";
brBCList[C_Association, e_] := brBCList[C, e, C["F"]];
brBCList[C_Association, e_, F_Association] := afFromY[F, C["B"], brBCNorm[C, e], C["Y"]];

brBCExpr::usage = "brBCExpr[C, a] is the coordinate list a as the expression Sum c_i Y^i/d_i in x, Y and theta.";
brBCExpr[C_Association, a_List] := brBCNorm[C, afToY[C["B"], a, C["Y"]]];

brBCD::usage = "brBCD[C, e] is the derivative of e (in x, Y, theta) on K: D x = 1, D Y = sY Y, D theta = C[\"Dth\"], in normal form.";
brBCD[C_Association, e_] := brBCNorm[C, D[e, C["x"]] + D[e, C["Y"]] C["sY"] C["Y"] + D[e, C["th"]] C["Dth"]];

brBCSurface::usage = "brBCSurface[C, e] writes the expression e in x, Y and theta in the user's surface form: Y -> C[\"Ysurf\"], theta -> the kernel.";
brBCSurface[C_Association, e_] := e /. {C["Y"] -> C["Ysurf"], C["th"] -> C["kernel"]};

brBCSplitDen::usage = "brBCSplitDen[C, e] is {N, d}: e = N/d with d the theta-dependent part of the denominator (a polynomial in x and theta over Q, free of Y) and N in normal form, polynomial in theta.";
brBCSplitDen[C_Association, e_] := Module[{th = C["th"], t = Together[brBCNorm[C, e]], den, dth},
  den = Denominator[t];
  dth = Times @@ (#[[1]]^#[[2]] & /@ Select[tragerFactorList[den], !FreeQ[#[[1]], th] &]);
  {brBCNorm[C, e dth], Expand[dth]}];

brBCInverseMod::usage = "brBCInverseMod[C, c, V] is the inverse of c (in k0[theta]) modulo the polynomial V (in k[theta], normal) in k0[theta], as an expression of theta-degree < deg V; an internal error when c is not invertible modulo V.";
brBCInverseMod[C_Association, c_, V_] := Module[{F = C["F"], B = C["B"], th = C["th"], eg},
  eg = brK0PolyExtGCD[F, B, brBCList[C, c], brBCList[C, V], th];
  If[brK0Degree[eg[[1]], th] =!= 0, throwInternal["brBCNotInvertibleMod", <|"c" -> c, "V" -> V|>]];
  brBCExpr[C, eg[[2]]]];

brBCRemTheta::usage = "brBCRemTheta[C, a, V] is the remainder of a (in k0[theta]) by V (in k[theta], free of Y): coordinatewise PolynomialRemainder in theta.";
brBCRemTheta[C_Association, a_, V_] := Module[{th = C["th"], Y = C["Y"], t = brBCNorm[C, a]},
  Sum[With[{c = Together[Coefficient[t, Y, i]]},
      Together[PolynomialRemainder[Expand[Numerator[c]], V, th]/Denominator[c]]] Y^i, {i, 0, C["m"] - 1}]];

brBCHermite::usage = "brBCHermite[C, e] is Hermite reduction in theta over k0 at the normal places: <|\"H\" -> H, \"R\" -> R|> with e = D H + R and the normal part of R's denominator squarefree.";
brBCHermite[C_Association, e0_] := Module[{th = C["th"], e = brBCNorm[C, e0], H = 0, N, d, fl, mult, V, j, U, Bv, it = 0},
  While[True,
    If[++it > 200, throwInternal["brBCHermiteIterationBound"]];
    {N, d} = brBCSplitDen[C, e];
    fl = Select[tragerFactorList[d], !FreeQ[#[[1]], th] && #[[2]] >= 2 &&
                 !(C["kind"] === "Exp" && Together[#[[1]]/th] === 1) &];
    If[fl === {}, Break[]];
    j = Max[Last /@ fl];
    V = Expand[Times @@ (First /@ Select[fl, #[[2]] === j &])];
    U = Together[d/V^j];
    (* B ((1 - j) U DV) == N (mod V) *)
    Bv = brBCRemTheta[C, brBCNorm[C, N brBCInverseMod[C, brBCNorm[C, (1 - j) U brBCD[C, V]], V]], V];
    H = brBCNorm[C, H + Bv/V^(j - 1)];
    e = brBCNorm[C, e - brBCD[C, Bv/V^(j - 1)]]];
  <|"H" -> H, "R" -> e|>];

brBCLaurent::usage = "brBCLaurent[C, e] reads e (in k0[theta, 1/theta]: its denominator a power of theta times an x-polynomial) as {{j, e_j}, ...}, e_j in x and Y; an internal error when e has a normal denominator.";
brBCLaurent[C_Association, e0_] := Module[{th = C["th"], e = brBCNorm[C, e0], t, q = 0, cl, den},
  If[e === 0, Return[{}]];
  t = Together[e];
  While[!FreeQ[Denominator[t], th], t = Together[t th]; q++;
    If[q > 200, throwInternal["brBCNotLaurent", <|"e" -> e0|>]]];
  den = Denominator[t];
  cl = CoefficientList[Expand[Numerator[t]], th];
  Select[Table[{k - 1 - q, brBCNorm[C, cl[[k]]/den]}, {k, Length[cl]}], #[[2]] =!= 0 &]];

brBCSplit::usage = "brBCSplit[C, R] splits the Hermite remainder R as L + a/d: <|\"L\" -> {{j, e_j}}, \"a\" -> a, \"d\" -> d|>, d the normal (squarefree) part of the theta-denominator, deg a < deg d, L a Laurent polynomial (exponential theta: the power of theta is special).";
brBCSplit[C_Association, R_] := Module[{th = C["th"], N, d, q = 0, dn, a, inv, L},
  {N, d} = brBCSplitDen[C, R];
  dn = d;
  If[C["kind"] === "Exp",
    While[!FreeQ[dn, th] && Expand[dn /. th -> 0] === 0, dn = Expand[Together[dn/th]]; q++]];
  If[FreeQ[dn, th], Return[<|"L" -> brBCLaurent[C, R], "a" -> 0, "d" -> 1|>]];
  (* a th^q == N (mod dn): th^q is invertible modulo dn in k[theta] *)
  inv = PolynomialExtendedGCD[th^q, dn, th];
  a = brBCRemTheta[C, brBCNorm[C, N inv[[2, 1]]/inv[[1]]], dn];
  L = brBCLaurent[C, brBCNorm[C, R - a/dn]];
  <|"L" -> L, "a" -> a, "d" -> dn|>];

brBCLogPart::usage = "brBCLogPart[C, a, d] is the Rothstein-Trager log part of a/d (d squarefree normal in k[theta], deg a < deg d) over k0: <|\"Logs\" -> {{c, u}, ..., {1, rootsum, \"RootSum\"}}, \"Rem\" -> a/d - Sum c Du/u|> (Rem in k0[theta, 1/theta]). Throws NonElementary when the monic R(z) has a non-constant coefficient.";
brBCLogPart[C_Association, a_, d0_] := Module[
  {th = C["th"], Y = C["Y"], m = C["m"], B = C["B"], F = C["F"], z = Unique["Bronstein`Private`brZc"],
   d, lcd, Dd, R, cl, lc, lcInv, P, facs, logs = {}, sumDl, v, mu, t, Ft, val, w, dL},
  If[Together[a] === 0, Return[<|"Logs" -> {}, "Rem" -> 0|>]];
  lcd = Coefficient[d0, th, Exponent[d0, th]];
  d = Expand[Together[d0/lcd]];
  Dd = brBCD[C, d];
  R = Expand[Resultant[Expand[Numerator[Together[a - z Dd]]], d, th]/
             Denominator[Together[a - z Dd]]^Exponent[d, th]];
  cl = brBCNorm[C, #] & /@ CoefficientList[R, z];
  While[cl =!= {} && Last[cl] === 0, cl = Most[cl]];
  If[Length[cl] <= 1, throwInternal["brBCResultantConstant", <|"a" -> a, "d" -> d|>]];
  lc = brBCList[C, Last[cl]];
  lcInv = brK0Inverse[F, B, lc];
  P = 0;
  Do[With[{cc = brK0Times[F, B, brBCList[C, cl[[k]]], lcInv]},
      If[AnyTrue[Rest[cc], Together[#] =!= 0 &] || !FreeQ[cc[[1]], C["x"] | th | Y],
        throwNonElementary["non-constant residues (resultant criterion over the algebraic base, [B] Prop 2.9)"]];
      P = P + Together[cc[[1]]] z^(k - 1)],
    {k, Length[cl]}];
  facs = brShapeBResidueRoots[P, z];
  dL = brBCList[C, d];
  sumDl = afZero[B];
  Do[
    mu = fac[[1]];
    If[fac[[2]] === 1,
      val = Together[-Coefficient[mu, z, 0]/Coefficient[mu, z, 1]];
      v = brK0PolyGCD[F, B, dL, brBCList[C, a - val Dd], th];
      If[brK0Degree[v, th] > 0,
        AppendTo[logs, {val, brBCExpr[C, v]}];
        sumDl = afPlus[F, sumDl, afScale[F, val,
          brK0Times[F, B, brBCList[C, brBCD[C, brBCExpr[C, v]]], brK0Inverse[F, B, v]]]]],
      t = Unique["Bronstein`Private`brTc"];
      Ft = fieldFn[{th}, {t, mu /. z -> t, Root[Function[Evaluate[mu /. z -> #]], 1]}];
      v = brK0PolyGCD[Ft, B, dL, brBCList[C, a - t Dd, Ft], th];
      If[brK0Degree[v, th] > 0,
        w = afScale[Ft, t, brK0Times[Ft, B, brBCList[C, brBCD[C, afToY[B, v, Y]], Ft], brK0Inverse[Ft, B, v]]];
        sumDl = afPlus[F, sumDl, brK0Trace[Ft, #] & /@ w];
        With[{vt = brBCNorm[C, afToY[B, v, Y]], mf = Function[Evaluate[mu /. z -> #]]},
          AppendTo[logs, {1, RootSum[mf, Function[Evaluate[# Log[vt /. t -> #]]]], "RootSum"}]]]],
    {fac, facs}];
  <|"Logs" -> logs, "Rem" -> brBCNorm[C, a/d - brBCExpr[C, sumDl]]|>];

brBCIntegrateK0::usage = "brBCIntegrateK0[C, e] is an elementary integral of e in k0 = Q(x)(Y) (free of theta), as a SURFACE expression: Trager's algorithm on the user's radical (the rational part by the rational integrator). Trager's certificate becomes NonElementary, its decline a Decline.";
brBCIntegrateK0[C_Association, e_] := Module[{x = C["x"], s, r},
  If[!FreeQ[e, C["th"]], throwInternal["brBCIntegrateK0Theta", <|"e" -> e|>]];
  If[Together[e] === 0, Return[0]];
  s = Together[brBCSurface[C, e]];
  If[FreeQ[e, C["Y"]], Return[brIntegrateTranscendental[s, x]]];
  r = Trager`TragerIntegrate[s, x, "Trager"];
  Which[
    MatchQ[r, {"not elementary", _}], throwNonElementary["algebraic base (Trager): " <> r[[2]]],
    MatchQ[r, {"implementation incomplete", _}], throwIncomplete["algebraic base (Trager): " <> r[[2]]],
    MatchQ[r, {"failed", _}], throwUnimplemented["Trager", "algebraic base (Trager): " <> r[[2]]],
    True, r]];

brBCFromSurface::usage = "brBCFromSurface[C, e] writes the surface expression e (in x and powers of the user's radicand b) in normal form in x and Y: b^(k/n) = (phi Y)^k with phi = s/q (b^(1/n) = yh/q, yh = s Y).";
brBCFromSurface[C_Association, e_] := Module[{x = C["x"], b = C["b"], n = C["m"], phi = C["phi"], Y = C["Y"]},
  brBCNorm[C, e /. Power[c_, r_Rational] /; !FreeQ[c, x] :>
    Which[Together[c - b] === 0, (phi Y)^(r n), Together[c b - 1] === 0, (phi Y)^(-r n),
      True, throwInternal["brBCFromSurfaceRadical", <|"c" -> c|>]]]];

brBCLogFree::usage = "brBCLogFree[e] is the part of the sum e free of Log, ArcTan, ArcTanh, ArcCot and RootSum terms (the algebraic part of an elementary integral).";
brBCLogFree[e_] := Module[{t = Expand[e]},
  If[Head[t] === Plus, Select[t, FreeQ[#, Log | ArcTan | ArcTanh | ArcCot | RootSum] &],
    If[FreeQ[t, Log | ArcTan | ArcTanh | ArcCot | RootSum], t, 0]]];

brBCLimited::usage = "brBCLimited[C, r] solves the limited integration problem r = D w + c eta in k0 (eta = D theta, theta = Log[h]): {w, c} with w in x and Y, c rational. Throws NonElementary when there is none (Trager proves r non-integrable, or r - D v is not a constant multiple of eta: see the header), Decline when Trager is undecided.";
brBCLimited[C_Association, r_] := Module[{F0, v, L, c, eta = C["Dth"]},
  If[Together[r] === 0, Return[{0, 0}]];
  F0 = brBCIntegrateK0[C, r];
  v = brBCFromSurface[C, brBCLogFree[F0]];
  L = brBCNorm[C, r - brBCD[C, v]];
  If[L === 0, Return[{v, 0}]];
  c = Together[brBCNorm[C, L/eta]];
  If[!FreeQ[c, C["x"] | C["Y"] | C["th"]] || !NumericQ[c],
    throwNonElementary["no algebraic part: the limited integration problem has no solution over the algebraic base"]];
  {v, c}];

brBCPolyExp::usage = "brBCPolyExp[C, L] integrates the Laurent polynomial Sum e_j theta^j (exponential theta = E^h, L = {{j, e_j}}) over k0: b_j' + j h' b_j = e_j by the Risch DE on the curve (RDEAlgebraic.m) for j != 0, the theta^0 coefficient in k0. It returns {symbolic part (in x, Y, theta), surface part}.";
brBCPolyExp[C_Association, L_List] := Module[{th = C["th"], B = C["B"], F0 = 0, alg = 0, r, e0 = 0},
  Do[If[lj[[1]] === 0, e0 = lj[[2]],
      r = brAlgRDE[B, afFromY[fieldQ[], B, brBCNorm[C, lj[[1]] C["Dh"]], C["Y"]],
                      afFromY[fieldQ[], B, lj[[2]], C["Y"]]];
      Switch[First[r],
        "Solution", alg = alg + brBCNorm[C, afToY[B, r[[2]], C["Y"]]] th^lj[[1]],
        "NoSolution", throwNonElementary["no algebraic part: the Risch DE D b + j h' b = c has no solution over the algebraic base ([B] Part 3)"],
        _, throwUnimplemented["RDEAlgebraic:undecided", "Risch DE on the curve undecided: " <> ToString[Last[r]]]]],
    {lj, L}];
  {alg, brBCIntegrateK0[C, e0]}];

brBCPolyLog::usage = "brBCPolyLog[C, L] integrates the polynomial Sum e_j theta^j (logarithmic theta = Log[h], L = {{j, e_j}}, j >= 0) over k0: top down, the limited integration problems e_j - (j+1) b_(j+1) eta = D w + c eta (brBCLimited), then the theta^0 coefficient in k0. It returns {symbolic part, surface part}.";
brBCPolyLog[C_Association, L_List] := Module[{th = C["th"], es = Association[Rule @@@ L], p, beta = 0, alg = 0, vc, eta = C["Dth"]},
  If[L === {}, Return[{0, 0}]];
  If[Min[First /@ L] < 0, throwInternal["brBCLogNegativePower"]];
  p = Max[First /@ L];
  Do[
    vc = brBCLimited[C, brBCNorm[C, Lookup[es, j, 0] - (j + 1) beta eta]];
    (* b_(j+1) = beta + c/(j+1) is now fixed *)
    alg = alg + brBCNorm[C, beta + vc[[2]]/(j + 1)] th^(j + 1);
    beta = vc[[1]],
    {j, p, 1, -1}];
  alg = alg + beta th;
  {alg, brBCIntegrateK0[C, brBCNorm[C, Lookup[es, 0, 0] - beta eta]]}];

brShapeBMono::usage = "brShapeBMono[T, f, x] integrates the curve-monomial integrand (tower record T, shape \"BCurve\") over k0 = Q(x)(Y): Hermite, Rothstein-Trager and the polynomial part over k0 (see the header). It returns <|\"F\" -> antiderivative, \"verified\" -> False|> or throws.";
brShapeBMono[T_Association, f_, x_Symbol] := Module[{C, her, sp, lp, L, pe, symb, surf},
  C = brBCContext[T];
  her = brBCHermite[C, C["f"]];
  sp = brBCSplit[C, her["R"]];
  lp = brBCLogPart[C, sp["a"], sp["d"]];
  L = Normal[Merge[{Association[Rule @@@ sp["L"]], Association[Rule @@@ brBCLaurent[C, lp["Rem"]]]},
    brBCNorm[C, Total[#]] &]];
  L = Select[List @@@ L, #[[2]] =!= 0 &];
  pe = If[C["kind"] === "Exp", brBCPolyExp[C, L], brBCPolyLog[C, L]];
  symb = her["H"] + pe[[1]] +
    Total[If[Length[#] === 3, #[[1]] #[[2]], #[[1]] Log[#[[2]]]] & /@ lp["Logs"]];
  <|"F" -> brBCSurface[C, symb] + pe[[2]], "verified" -> False|>];
