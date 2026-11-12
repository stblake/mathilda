(* ::Package:: *)

(* ========================================================================== *)
(* BaseAlgebraic.m -- Shape B, the component method: a radical y^n = g(x) at *)
(* the BOTTOM of a Log/Exp tower whose kernel arguments are free of y         *)
(* (BRONSTEIN_PLAN.md "Shape B").                                             *)
(*                                                                            *)
(* In Bronstein's terms (1990, Parts 2-3) this is the transcendental Risch    *)
(* algorithm over the algebraic base k0 = Q(x)(y) with n = 1: k = Q(x)(t_i)   *)
(* is not algebraically closed in K, so y joins the base.  With the kernel    *)
(* arguments free of y every D t_i lies in k, and the integrand splits into   *)
(* components a_i y^i ([T] App A Sec 1; the Galois projection y -> zeta y     *)
(* fixes every t_i and commutes with D).  Each lives on its reduced curve     *)
(* Y^m = G(x), Y = y^i/s (componentCurve), as c Y with c in k(theta), and the *)
(* space k(theta) Y is closed under D: D(b Y) = (D b + s b) Y, s = G'/(m G)   *)
(* in Q(x).  So, for one component, with theta the top monomial:             *)
(*                                                                            *)
(*   1. Hermite at the normal places: brHermiteReduce (Shape A's, with G free *)
(*      of theta its Prop 2.2 invertibility needs only V normal).             *)
(*   2. The proper normal part a/d: Rothstein-Trager over the FIELD k0.  The  *)
(*      residues are the roots of brResultantR's R(z) (whose roots are the    *)
(*      residues times m-th roots of unity, so constancy is equivalent: the   *)
(*      [B] Prop 2.9 certificate); each root z0 gives the logand               *)
(*      v = gcd_{k0[theta]}(d, a Y - z0 D d), trivial unless z0 is a residue  *)
(*      (K0.m).  The remainder h - Sum z0 Dv/v is computed exactly in k0: a    *)
(*      Laurent polynomial in theta, whose weight-0 part (for an exponential  *)
(*      theta, -Sum z0 deg(v) eta') is integrated in k.                       *)
(*   3. The polynomial part.  Weight >= 1 has no nonzero constants, so the    *)
(*      coefficients are unique and need no limited integration:              *)
(*        log theta:  D(beta_j theta^j Y) = (D_s beta_j) theta^j Y            *)
(*                                          + j eta beta_j theta^(j-1) Y,     *)
(*                    top down beta_j' + s beta_j = e_j - (j+1) eta beta_(j+1) *)
(*                    for j >= 1 (a Risch DE over k, brRischDE);              *)
(*        exp theta:  beta_j' + (s + j eta') beta_j = e_j for j != 0;         *)
(*      the theta^0 coefficient is a full integral over the lower tower: the  *)
(*      same method one level down, ending in Trager's TragerIntegrateCurve.  *)
(*                                                                            *)
(* Certificates come only from theorems over k0 (constants Q: the reduced     *)
(* curve is absolutely irreducible): the residue criterion, a certified       *)
(* "no solution" of a Risch DE over k (brRischDE), and the lower level's or   *)
(* Trager's own certificates.  Every answer is verified by the caller.        *)
(*                                                                            *)
(* A component result is <|"Alg" -> e, "Logs" -> {{c, u}, ...,                *)
(* {c, rootsum, "RootSum"}, ...}, "Extra" -> e0, "Y" -> Y|>: the integral of  *)
(* c Y is Alg + Sum c Log[u] + Sum c rootsum + Extra, in x, the tower symbols *)
(* and the symbol Y (Extra is free of Y).                                     *)
(* ========================================================================== *)

brShapeBCurve::usage = "brShapeBCurve[m, G, x] is Shape B's curve record of Y^m = G(x): <|\"m\", \"G\", \"x\", \"s\" -> G'/(m G) (D Y = s Y), \"B\" -> curveBasis[m, G, x]|>.";
brShapeBCurve[m_Integer, G_, x_Symbol] := <|"m" -> m, "G" -> Expand[G], "x" -> x,
  "s" -> Together[D[G, x]/(m G)], "B" -> curveBasis[m, Expand[G], x]|>;

brShapeBTop::usage = "brShapeBTop[T] is the tower record T with its LAST monomial as theta (\"theta\", \"kind\", \"lower\" set), the form the theta-level routines (brHermiteReduce, brSplit, brResultantR) take.";
brShapeBTop[T_Association] := Join[T, <|"theta" -> Last[T["syms"]], "kind" -> Last[T["kinds"]],
  "lower" -> Most[T["syms"]]|>];

brShapeBRest::usage = "brShapeBRest[T] is the tower record of the monomials of T below the last one (shape \"B\", theta None).";
brShapeBRest[T_Association] := Module[{monos = Most[T["monos"]]},
  <|"x" -> T["x"], "monos" -> monos, "syms" -> (#["sym"] & /@ monos), "kinds" -> (#["kind"] & /@ monos),
    "args" -> (#["arg"] & /@ monos), "rules" -> brRules[T["x"], monos], "back" -> brBack[monos],
    "kernels" -> (#["kernel"] & /@ monos), "f" -> None, "radical" -> T["radical"], "shape" -> "B",
    "theta" -> None, "kind" -> None, "lower" -> (#["sym"] & /@ monos)|>];

brShapeBResult::usage = "brShapeBResult[alg, logs, extra, Y] is a component result Association (see the header).";
brShapeBResult[alg_, logs_List, extra_, Y_Symbol] := <|"Alg" -> alg, "Logs" -> logs, "Extra" -> extra, "Y" -> Y|>;

brShapeBJoin::usage = "brShapeBJoin[r1, r2] is the sum of two component results with the same symbol Y.";
brShapeBJoin[r1_Association, r2_Association] := (
  If[r1["Y"] =!= r2["Y"], throwInternal["brShapeBJoinY"]];
  brShapeBResult[r1["Alg"] + r2["Alg"], Join[r1["Logs"], r2["Logs"]], r1["Extra"] + r2["Extra"], r1["Y"]]);

brShapeBExpression::usage = "brShapeBExpression[r] is the integral a component result r stands for, as an expression in x, the tower symbols and r[\"Y\"].";
brShapeBExpression[r_Association] := r["Alg"] + r["Extra"] +
  Total[If[Length[#] === 3, #[[1]] #[[2]], #[[1]] Log[#[[2]]]] & /@ r["Logs"]];

(* -------------------------------------------------------------------------- *)
(* The base of the recursion: no monomial left, c in Q(x).                    *)
(* -------------------------------------------------------------------------- *)

brShapeBBase::usage = "brShapeBBase[cv, c, Y] integrates c(x) Y on the curve cv (Y^m = G) by Trager's algorithm (TragerIntegrateCurve), as a component result in the symbol Y; Trager's certificate becomes NonElementary, a decline a Decline.";
brShapeBBase[cv_Association, c_, Y_Symbol] := Module[{x = cv["x"], r, Yt},
  If[Together[c] === 0, Return[brShapeBResult[0, {}, 0, Y]]];
  r = Trager`TragerIntegrateCurve[Together[c], x, cv["m"], cv["G"], "Trager"];
  Which[
    MatchQ[r, {"not elementary", _}], throwNonElementary["Trager (algebraic base): " <> r[[2]]],
    MatchQ[r, {"failed", _}], throwDecline["Trager (algebraic base): " <> r[[2]]],
    !AssociationQ[r], throwInternal["brShapeBBaseResult", <|"r" -> r|>]];
  Yt = r["Y"];
  brShapeBResult[r["Algebraic"] /. Yt -> Y, r["Logs"] /. Yt -> Y, 0, Y]];

(* -------------------------------------------------------------------------- *)
(* The theta level.                                                           *)
(* -------------------------------------------------------------------------- *)

brShapeBSplitRemainder::usage = "brShapeBSplitRemainder[TB, r] splits the Hermite remainder r (an element of k(theta) whose normal denominator is squarefree) as r = L + a/d: L a Laurent polynomial in theta (a polynomial when theta is a logarithm), d the normal part of the denominator (monic-free, in k[theta]) and deg a < deg d. It returns <|\"L\" -> {{j, e_j}, ...} (nonzero coefficients, e_j in k), \"a\" -> a, \"d\" -> d|>.";
brShapeBSplitRemainder[TB_Association, r0_] := Module[{th = TB["theta"], r = Together[r0], A, den, q, d, a, L, cl},
  If[r === 0, Return[<|"L" -> {}, "a" -> 0, "d" -> 1|>]];
  A = Expand[Numerator[r]]; den = Expand[Denominator[r]];
  (* the special part of the denominator: the power of theta (exponential case) *)
  q = 0;
  If[TB["kind"] === "Exp",
    While[!FreeQ[den, th] && Expand[den /. th -> 0] === 0, den = Expand[Together[den/th]]; q++]];
  If[FreeQ[den, th],
    a = 0; d = 1,
    d = den;
    a = brHermiteSolve[TB, A, th^q, d]];                  (* a th^q == A (mod d), deg a < deg d *)
  L = Together[(Together[A/(th^q den)] - a/d) th^q];
  If[!FreeQ[Denominator[L], th], throwInternal["brShapeBSplitNotLaurent", <|"r" -> r0, "L" -> L|>]];
  cl = CoefficientList[Expand[Numerator[L]], th];
  <|"L" -> Select[Table[{j - 1 - q, Together[cl[[j]]/Denominator[L]]}, {j, Length[cl]}], #[[2]] =!= 0 &],
    "a" -> a, "d" -> d|>];

brShapeBResidueRoots::usage = "brShapeBResidueRoots[P, z] is the list of the irreducible factors of the residue polynomial P(z) (constant coefficients, nonzero roots only) as {mu(z), degree}.";
brShapeBResidueRoots[P_, z_Symbol] := Module[{p = brNonzeroResiduePolynomial[P, z]},
  If[FreeQ[p, z], Return[{}]];
  {#[[1]], Exponent[#[[1]], z]} & /@ Select[FactorList[p], !FreeQ[#[[1]], z] &]];

brShapeBLogPart::usage = "brShapeBLogPart[TB, cv, a, d, Y] is the logarithmic part of (a/d) Y over k0 = k(Y) (Rothstein-Trager in theta over the field k0, K0.m): it returns <|\"Logs\" -> {{c, u}, ...}, \"Rem\" -> rem|> with rem = (a/d) Y - Sum c Du/u, computed exactly as an element of k0 (coordinate list on cv[\"B\"]) that is a Laurent polynomial in theta. Throws NonElementary when the residues are not constant ([B] Prop 2.9 over k0).";
brShapeBLogPart[TB_Association, cv_Association, a_, d_, Y_Symbol] := Module[
  {th = TB["theta"], B = cv["B"], F = brK0Field[TB], z = Unique["Bronstein`Private`brZb"], R, crit, P, facs,
   h, dEl, Dd, logs = {}, sumDl, t, Ft, v, aY, w, mu, rootsQ, val},
  h = brK0YMultiple[B, Together[a/d]];
  If[Together[a] === 0, Return[<|"Logs" -> {}, "Rem" -> afZero[B]|>]];
  R = brResultantR[TB, cv, a/d, z];
  crit = brResultantCriterion[TB, R, z];
  If[First[crit] === "NonConstant",
    throwNonElementary["non-constant residues (resultant criterion over the algebraic base, [B] Prop 2.9)"]];
  P = crit[[2]];
  facs = brShapeBResidueRoots[P, z];
  dEl = brK0Scalar[B, d];
  Dd = brK0Scalar[B, Together[brD[TB, d]]];
  sumDl = afZero[B];
  Do[
    mu = fac[[1]];
    If[fac[[2]] === 1,
      (* a rational residue *)
      val = Together[-Coefficient[mu, z, 0]/Coefficient[mu, z, 1]];
      v = brK0PolyGCD[F, B, dEl, afMinus[F, brK0YMultiple[B, a], afScale[F, val, Dd]], th];
      If[brK0Degree[v, th] > 0,
        AppendTo[logs, {val, brK0ToY[B, v, Y]}];
        sumDl = afPlus[F, sumDl, afScale[F, val, brK0Times[F, B, brK0D[TB, F, B, v], brK0Inverse[F, B, v]]]]],
      (* an algebraic residue t, mu(t) = 0: the logand over k0(t), the sum over the conjugates as a RootSum *)
      t = Unique["Bronstein`Private`brTb"];
      Ft = brK0Field[TB, {t, mu /. z -> t, Root[Function[Evaluate[mu /. z -> #]], 1]}];
      v = brK0PolyGCD[Ft, B, dEl, afMinus[Ft, brK0YMultiple[B, a], afScale[Ft, t, Dd]], th];
      If[brK0Degree[v, th] > 0,
        w = afScale[Ft, t, brK0Times[Ft, B, brK0D[TB, Ft, B, v], brK0Inverse[Ft, B, v]]];
        (* the sum over the conjugates of t: the trace from Q(t) down to Q *)
        sumDl = afPlus[F, sumDl, brK0Trace[Ft, #] & /@ w];
        With[{vt = brK0ToY[B, v, Y], mf = Function[Evaluate[mu /. z -> #]]},
          AppendTo[logs, {1, RootSum[mf, Function[Evaluate[(# Log[vt /. t -> #])]]], "RootSum"}]]]],
    {fac, facs}];
  <|"Logs" -> logs, "Rem" -> afMinus[F, h, sumDl]|>];

brShapeBLaurent::usage = "brShapeBLaurent[TB, cv, rem] reads the Laurent polynomial rem (an element of k0, its coordinates Laurent polynomials in theta) as {{j, e_j}, ...} for its weight-1 coordinate (rem = Sum e_j theta^j Y + ...) and the weight-0 part e0 (in k(theta)); any other coordinate is an internal error (the log part of a weight-1 integrand has weights 0 and 1 only).";
brShapeBLaurent[TB_Association, cv_Association, rem_List] := Module[{th = TB["theta"], B = cv["B"], m = cv["m"], c1, q, cl, t},
  If[m > 2 && AnyTrue[rem[[3 ;;]], Together[#] =!= 0 &], throwInternal["brShapeBLogWeight", <|"rem" -> rem|>]];
  c1 = Together[If[m >= 2, rem[[2]]/B["d"][[2]], 0]];
  If[c1 === 0, Return[{{}, Together[rem[[1]]]}]];
  t = Together[c1];
  q = 0;
  While[!FreeQ[Denominator[t], th], t = Together[t th]; q++;
    If[q > 200, throwInternal["brShapeBNotLaurent", <|"rem" -> rem|>]]];
  cl = CoefficientList[Expand[Numerator[t]], th];
  {Select[Table[{j - 1 - q, Together[cl[[j]]/Denominator[t]]}, {j, Length[cl]}], #[[2]] =!= 0 &], Together[rem[[1]]]}];

brShapeBRDE::usage = "brShapeBRDE[TL, f, g] solves D y + f y = g for y in the field of the tower record TL (Q(x) and its kernels; f, g in its symbols) by the certified Risch DE solver: the solution y in the symbols, or a throw: NonElementary for a certified \"no solution\" (with the given reason), Decline when undecided.";
brShapeBRDE[TL_Association, f_, g_, why_String] := Module[{r},
  If[Together[g] === 0, Return[0]];
  r = brRischDE[brToSurface[TL, f], brToSurface[TL, g], TL["x"], TL["kernels"]];
  Switch[First[r],
    "Solution", Together[brFromSurface[TL, r[[2]]]],
    "NoSolution", throwNonElementary[why],
    _, throwDecline["Risch DE over the base undecided: " <> ToString[Last[r]]]]];

brShapeBPolyPart::usage = "brShapeBPolyPart[TB, cv, L, Y] integrates the weight-1 Laurent polynomial Sum_j e_j theta^j Y (L = {{j, e_j}, ...}) over the tower TB: the coefficients theta^j, j != 0, by Risch DEs over k (the weight forbids constants, so each is unique), the theta^0 coefficient by the component method one level down. It returns a component result.";
brShapeBPolyPart[TB_Association, cv_Association, L_List, Y_Symbol] := Module[
  {th = TB["theta"], s = cv["s"], TL = brShapeBRest[TB], eta, es, jmax, jmin, beta, e0, alg = 0, low},
  es = Association[Rule @@@ L];
  eta = brEta[TB];
  If[L === {}, Return[brShapeBResult[0, {}, 0, Y]]];
  jmax = Max[First /@ L]; jmin = Min[First /@ L];
  If[TB["kind"] === "Log",
    If[jmin < 0, throwInternal["brShapeBLogNegativePower"]];
    beta = 0;                                                  (* beta_(j+1) *)
    Do[
      beta = brShapeBRDE[TL, s, Together[Lookup[es, j, 0] - (j + 1) eta beta],
        "no algebraic part: the coefficient Risch DE D b + s b = c has no solution over the base"];
      alg = alg + beta th^j,
      {j, jmax, 1, -1}];
    e0 = Together[Lookup[es, 0, 0] - eta beta],
    (* Exp *)
    Do[If[j =!= 0,
        beta = brShapeBRDE[TL, Together[s + j eta], Lookup[es, j, 0],
          "no algebraic part: the coefficient Risch DE D b + (s + j eta') b = c has no solution over the base"];
        alg = alg + beta th^j],
      {j, jmin, jmax}];
    e0 = Lookup[es, 0, 0]];
  low = brShapeBComponent[TL, cv, e0, Y];
  brShapeBJoin[brShapeBResult[Together[alg] Y, {}, 0, Y], low]];

brShapeBComponent::usage = "brShapeBComponent[T, cv, c, Y] integrates the component c Y (c in the field of the tower record T, Y^m = G(x) the curve cv, D Y = s Y) by the component method (see the header): a component result in the symbol Y, or a throw (NonElementary with its theorem, Decline, Internal).";
brShapeBComponent[T_Association, cv_Association, c0_, Y_Symbol] := Module[{c = Together[c0], TB, hr, sp, lp, lau, poly, extra},
  If[c === 0, Return[brShapeBResult[0, {}, 0, Y]]];
  If[T["syms"] === {}, Return[brShapeBBase[cv, c, Y]]];
  TB = brShapeBTop[T];
  hr = brHermiteReduce[TB, cv, c];
  sp = brShapeBSplitRemainder[TB, hr["Remainder"]];
  lp = brShapeBLogPart[TB, cv, sp["a"], sp["d"], Y];
  lau = brShapeBLaurent[TB, cv, lp["Rem"]];
  (* weight 0 of the log part's remainder (exponential theta: -Sum z0 deg v eta'),
     an element of k(theta), integrated by the transcendental algorithm *)
  extra = If[Together[lau[[2]]] === 0, 0,
    brFromSurface[TB, brIntegrateTranscendental[brToSurface[TB, -lau[[2]]], TB["x"]]]];
  (* the polynomial part: the split's Laurent part plus the log remainder's weight-1 part *)
  poly = Normal[Merge[{Association[Rule @@@ sp["L"]], Association[Rule @@@ lau[[1]]]}, Together[Total[#]] &]];
  poly = Select[List @@@ poly, #[[2]] =!= 0 &];
  brShapeBJoin[
    brShapeBResult[Together[hr["Algebraic"]] Y, lp["Logs"], extra, Y],
    brShapeBPolyPart[TB, cv, poly, Y]]];

(* -------------------------------------------------------------------------- *)
(* The whole integrand.                                                       *)
(* -------------------------------------------------------------------------- *)

brShapeB::usage = "brShapeB[T, f, x] integrates the Shape B integrand f (tower record T, shape \"B\": radical y^n = b(x), kernel arguments free of y) by the component method: coordinate 0 by the transcendental Risch algorithm, every component a_i y^i on its reduced curve over Q(x) by brShapeBComponent. It returns <|\"F\" -> antiderivative, \"verified\" -> False|> (the caller verifies).";
brShapeB[T_Association, f_, x_Symbol] := Module[{cf, n, g, total, cc, cv, res, Yexpr, rho, Y, comp},
  cf = brCanonicalForm[T];
  n = cf["n"]; g = cf["g"];
  total = brFromSurface[T, brIntegrateTranscendental[brToSurface[T, cf["coeffs"][[1]]], x]];
  Do[
    If[Together[cf["coeffs"][[i + 1]]] =!= 0,
      cc = componentCurve[n, g, i, x];
      If[cc["m"] === 1,
        (* the component is rational in the tower times the locally constant rho *)
        rho = (cf["yexpr"]^i/cc["s"])/cc["G"]^(1/cc["m"]);
        total = total + rho brFromSurface[T, brIntegrateTranscendental[
          brToSurface[T, Together[cf["coeffs"][[i + 1]] cc["s"] cc["G"]]], x]],
        If[!reducedCurveQ[cc["m"], cc["G"], x], throwDecline["component curve not reduced over Q(x)"]];
        cv = brShapeBCurve[cc["m"], cc["G"], x];
        Y = Unique["Bronstein`Private`brYB"];
        comp = brShapeBComponent[T, cv, Together[cf["coeffs"][[i + 1]] cc["s"]], Y];
        If[cc["exact"],
          rho = 1; Yexpr = Together[cf["yexpr"]^i/cc["s"]],
          Yexpr = cc["G"]^(1/cc["m"]);
          rho = (cf["yexpr"]^i/cc["s"])/Yexpr];
        total = total + rho (brShapeBExpression[comp] /. Y -> Yexpr)]],
    {i, 1, n - 1}];
  <|"F" -> brToSurface[T, total], "verified" -> False|>];
