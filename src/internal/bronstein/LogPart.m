(* ::Package:: *)

(* ========================================================================== *)
(* LogPart.m -- the logarithmic part of a component and its weight-1 form.    *)
(*                                                                            *)
(* L = Sum_j (beta_j/N_j) Log[u_j] has the residues of the component c Y      *)
(* ([B] Prop 2.13).  Its derivative has coordinates on every basis element;   *)
(* the component lives in coordinate 1 (weight 1 under sigma: Y -> zeta Y).   *)
(* The weight-1 projection L1 = (1/m) Sum_k zeta^-k sigma^k(L) has the same   *)
(* residues and derivative = coordinate 1 of L'.  The coordinate-0 part of    *)
(* u'/u is (1/m) Norm(u)'/Norm(u) (the trace of a logarithmic derivative), so *)
(* when coordinates 2..m-1 of L' vanish (always for m = 2),                   *)
(*     L1 = Sum_j (beta_j/N_j) (Log[u_j] - (1/m) Log[Norm(u_j)]),            *)
(* for m = 2 the familiar (beta/(2N)) Log[u/sigma(u)].                        *)
(* The remainder c - [L']_1 is then the algebraic part's right-hand side.     *)
(* ========================================================================== *)

brAfD::usage = "brAfD[T, cv, F, B, u] is the TOWER derivative of the element u (coordinates on the basis w_i = Y^i/d_i of B over F): (D u_i + u_i (i s - D d_i/d_i)), s = cv[\"s\"].";
brAfD[T_Association, cv_Association, F_Association, B_Association, u_List] := Module[{s = cv["s"]},
  Table[rc[F, brD[T, u[[i + 1]]] + u[[i + 1]] (i s - brD[T, B["d"][[i + 1]]]/B["d"][[i + 1]])], {i, 0, B["m"] - 1}]];

brLogDerivative::usage = "brLogDerivative[T, cv, F, B, u] is u'/u (tower derivative) as a coordinate list.";
brLogDerivative[T_Association, cv_Association, F_Association, B_Association, u_List] :=
  afTimes[F, B, brAfD[T, cv, F, B, u], rationalElementInverse[F, B, u]];

brLogPartDerivative::usage = "brLogPartDerivative[T, cv, F, B, terms] is Sum c (u'/u) over the log terms {{c, u}, ...} (coordinate list).";
brLogPartDerivative[T_Association, cv_Association, F_Association, B_Association, terms_List] :=
  Fold[afPlus[F, #1, afScale[F, #2[[1]], brLogDerivative[T, cv, F, B, #2[[2]]]]] &, afZero[B], terms];

brComponentLogPart::usage = "brComponentLogPart[T, cv, c] is the log part of the normal (reduced) component c Y: <|\"terms\" -> {{coef, u, normU}, ...} (L1 = Sum coef (Log[u] - Log[normU]/m), u a coordinate list, normU its norm), \"remainder\" -> c - [L']_1 (in k(theta)), \"F\" -> field, \"B\" -> basis|>. Throws NonElementary from the resultant criterion or a residue divisor of infinite order, Decline when special-place residues are nonzero (not yet supported) or a residue/specialization step is undecided.";
brComponentLogPart[T_Association, cv_Association, c_] := Module[
  {nr, spec, rs, F, B, terms = {}, Nj, u, Lp, rem, vars},
  If[Together[c] === 0, Return[<|"terms" -> {}, "remainder" -> 0, "F" -> None, "B" -> cv["B"]|>]];
  nr = brNormalResidues[T, cv, c];
  spec = If[T["kind"] === "Log", brLogInfinityResidues[T, cv, c], brExpSpecialResidues[T, cv, c]];
  If[!brSpecialResiduesZeroQ[T, spec],
    Return[brComponentLogPartSpecial[T, cv, c, nr, spec]]];
  If[nr["P"] === 1, Return[<|"terms" -> {}, "remainder" -> Together[c], "F" -> None, "B" -> cv["B"]|>]];
  rs = brResidueSetup[T, cv, c, nr];
  F = rs["Fk"]; B = rs["B"];
  Do[
    Nj = brTorsionOrder[T, cv, rs, rs["coeffs"][[All, j]]];
    u = brPrincipalGenerator[rs, rs["coeffs"][[All, j]], Nj];
    If[u === None,
      throwNonElementary["residue divisor is not torsion (N delta not principal at the only possible order N)",
        <|"N" -> Nj|>]];
    AppendTo[terms, {kc[F, rs["basis"][[j]]/Nj], u, afNorm[F, B, u]}],
    {j, Length[rs["basis"]]}];
  Lp = brLogPartDerivative[T, cv, F, B, terms[[All, {1, 2}]]];
  If[B["m"] > 2 && !AllTrue[Lp[[3 ;;]], rcZeroQ[F, #] &],
    throwDecline["log part with components of other weights (m >= 3)"]];
  brLogPartFinish[T, cv, c, F, B, terms]];

brLogPartFinish::usage = "brLogPartFinish[T, cv, c, F, B, terms] checks the log terms {{coef, u, normU}, ...} (coordinates on the theta-basis B over F) of the component c Y and returns the brComponentLogPart Association: coordinates 2..m-1 of L' must vanish (else Decline), and the remainder c - [L']_1 must be free of the residue-field generator (else Internal).";
brLogPartFinish[T_Association, cv_Association, c_, F_Association, B_Association, terms_List] := Module[{Lp, rem},
  Lp = brLogPartDerivative[T, cv, F, B, terms[[All, {1, 2}]]];
  If[B["m"] > 2 && !AllTrue[Lp[[3 ;;]], rcZeroQ[F, #] &],
    throwDecline["log part with components of other weights (m >= 3)"]];
  rem = rc[F, c - Lp[[2]]];
  If[fieldExtQ[F] && !FreeQ[rem, F["t"]], throwInternal["brLogPartRemainderNotRational", <|"rem" -> rem|>]];
  <|"terms" -> terms, "remainder" -> Together[rem], "F" -> F, "B" -> B|>];

brComponentLogPartSpecial::usage = "brComponentLogPartSpecial[T, cv, c, nr, spec] is brComponentLogPart when some residue at a special place is nonzero (genus 0 only, see Divisors.m): residues at the finite and special places, a Q-basis of all of them, divisors on the shifted frame, principal generators (genus 0: always principal) mapped back to the theta-basis.";
brComponentLogPartSpecial[T_Association, cv_Association, c_, nr_Association, spec_Association] := Module[
  {th = T["theta"], m = cv["m"], G = cv["G"], locs, z = nr["z"], polys, Pall, t = Unique["Bronstein`Private`brTau"],
   rf, F0, Fk, vars, finRoots, rhoLists, taus, allRes, qb, tc, N, E, Ep, En, Bp, rs, a, u = Unique["Bronstein`Private`brU"],
   fr, F, Bth, ideals, I, gen, Yh = Unique["Bronstein`Private`brYh"], Y = Unique["Bronstein`Private`brYs"], expr, terms = {}, nfin},
  If[curveGenus[m, G, th] =!= 0,
    throwDecline["residues at a special place on a curve of positive genus (not yet supported)"]];
  locs = Select[brSpecialLocations[T, cv, spec], #["ramified"] || #["coef"] =!= 0 || T["kind"] === "Exp" &];
  (* the field: finite residues and the m-th roots of every kappa *)
  polys = Join[If[nr["P"] === 1, {}, {nr["P"]}],
    (z^m - #["kappa"]) & /@ Select[locs, !#["ramified"] &]];
  Pall = Expand[Times @@ (First /@ Select[FactorList[Times @@ polys], !FreeQ[#[[1]], z] &])];
  rf = residueField[Pall, z, t];
  F0 = rf["F"];
  vars = Prepend[T["lower"], T["x"]];
  Fk = If[F0["Type"] === "NF", fieldFn[vars, {F0["t"], F0["mu"], F0["theta"]}], fieldFn[vars]];
  finRoots = If[nr["P"] === 1, {}, Select[rf["roots"], kc[F0, nr["P"] /. z -> #] === 0 &]];
  rhoLists = If[#["ramified"], {}, Select[rf["roots"], kc[F0, #^m] === kc[F0, #["kappa"]] &]] & /@ locs;
  rhoLists = MapThread[If[#1["ramified"], {}, Select[rf["roots"], Function[r, kc[F0, r^m - #1["kappa"]] === 0]]] &, {locs, rhoLists}];
  taus = kc[F0, #] & /@ brSpecialTaus[T, cv, locs, rhoLists];
  allRes = Join[finRoots, taus];
  If[AllTrue[allRes, # === 0 &], Return[<|"terms" -> {}, "remainder" -> Together[c], "F" -> None, "B" -> cv["B"]|>]];
  qb = qBasis[F0, allRes];
  (* residue data in theta and the frame *)
  tc = Together[c]; N = Expand[Numerator[tc]]; E = Expand[Denominator[tc]];
  Ep = Together[brD[T, E]];
  En = Numerator[Together[First[brSplit[T, E]]]];
  Bp = If[FreeQ[En, th], 1, Expand[Numerator[Together[En/PolynomialGCD[En, G]]]]];
  rs = <|"N" -> N, "Ep" -> Ep, "Bp" -> Bp|>;
  a = firstOrdinaryInteger[Join[{Expand[G], E}, If[T["kind"] === "Exp", {th}, {}]], th];
  If[a === 0, a = firstOrdinaryInteger[{Expand[G] , E, th}, th]];
  fr = brShiftFrame[T, cv, rs, a, u];
  F = Fk; Bth = cv["B"];
  ideals = brSpecialPlaceIdeals[F, fr, #[[1]], #[[2]]] & /@ Transpose[{locs, rhoLists}];
  ideals = Flatten[ideals, 1];
  nfin = Length[finRoots];
  Do[
    I = idealUnit[fr["Bu"]];
    Do[If[qb["coeffs"][[i, j]] =!= 0,
         I = idealTimes[F, fr["Bu"], I, idealPower[F, fr["Bu"],
           residueIdeal[F, fr["Bu"], finRoots[[i]], fr["Nh"], fr["Eh"], fr["Bph"]], qb["coeffs"][[i, j]]]]],
      {i, nfin}];
    Do[If[qb["coeffs"][[nfin + l, j]] =!= 0,
         I = idealTimes[F, fr["Bu"], I, idealPower[F, fr["Bu"], ideals[[l, 1]], qb["coeffs"][[nfin + l, j]]]]],
      {l, Length[ideals]}];
    gen = principalGenerator[F, fr["Bu"], I];
    If[gen === None, throwInternal["brGenus0DivisorNotPrincipal"]];
    (* back to theta: u -> 1/(theta - a), Yh -> Y (theta - a)^(-k) *)
    expr = Together[(afToY[fr["Bu"], gen, Yh] /. {u -> 1/(th - a)}) /. Yh -> Y (th - a)^(-fr["k"])];
    gen = afFromY[F, Bth, expr, Y];
    AppendTo[terms, {kc[F, qb["basis"][[j]]], gen, afNorm[F, Bth, gen]}],
    {j, Length[qb["basis"]]}];
  brLogPartFinish[T, cv, c, F, Bth, terms]];
