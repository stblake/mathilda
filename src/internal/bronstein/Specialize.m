(* ::Package:: *)

(* ========================================================================== *)
(* Specialize.m -- the order of a residue divisor class on a curve over the   *)
(*                 function field k = Q(x, t_1..t_r), by specialization.      *)
(*                                                                            *)
(* Trager's torsion bound uses reduction modulo primes of a curve over a      *)
(* number field ([T] Ch 6 Sec 2-3); over k the Weil bound does not apply      *)
(* directly ([B] p.130).  Specialize every transcendental of k (x and the     *)
(* lower monomials, independent indeterminates) to rationals sigma such that  *)
(* the curve keeps its shape -- deg G, the squarefree profile of G (hence the *)
(* genus), the degree and squarefreeness of the branch-free denominator Bp    *)
(* and its coprimality with G, every denominator nonzero, and the degree of   *)
(* every residue place A_alpha.  Specialization of a smooth model is a        *)
(* homomorphism of divisor class groups, injective on torsion in              *)
(* characteristic 0, so                                                       *)
(*   - delta_sigma of infinite order  =>  delta of infinite order;            *)
(*   - delta torsion of order N       =>  delta_sigma of order exactly N.     *)
(* Hence N_sigma (Trager's torsionOrderSpec on the specialized Q-curve) is the *)
(* only possible order, and principality of N_sigma delta, tested GENERICALLY *)
(* over k(t), decides.                                                        *)
(*                                                                            *)
(* Trager's frame arithmetic recomputes residue data with d/dtheta, while the *)
(* divisor here is defined with the tower derivative D.  The specialized      *)
(* remainder handed to Trager is therefore Nbar/Bp with                       *)
(*     Nbar == N Bp'/(D E)   (mod Bp),     ' = d/dtheta,                      *)
(* which has, at every root of Bp, the d/dtheta-residue N Y/(D E): the same    *)
(* places carry the same residues, so Delta(alpha) is unchanged.              *)
(* ========================================================================== *)

brSpecializationPoint::usage = "brSpecializationPoint[vars, try] is the deterministic try-th candidate specialization: a list of rules v -> rational for the symbols vars. The height grows with try: candidates of round r = Quotient[try - 1, brOption[\"SpecializationTries\"]] have numerators below brOption[\"SpecializationHeight\"] 2^r and denominators up to 7 2^r, so the search never cycles through a finite set.";
brSpecializationPoint[vars_List, try_Integer] := Module[{r = Quotient[try - 1, brOption["SpecializationTries"]], h, d},
  h = brOption["SpecializationHeight"] 2^r; d = 7 2^r;
  MapIndexed[
    #1 -> (Mod[7919 First[#2]^2 + 104729 try + 31 First[#2] try^2, h] + 2)/
          (Mod[13 try + 5 First[#2], d] + 1) &, vars]];

brSquarefreeProfile::usage = "brSquarefreeProfile[p, th] is the sorted list {{deg_th, multiplicity}, ...} of the squarefree factors of p involving th.";
brSquarefreeProfile[p_, th_Symbol] := Sort[{Exponent[#[[1]], th], #[[2]]} & /@
  Select[tragerSquareFreeList[Expand[p]], !FreeQ[#[[1]], th] &]];

brSpecializeValue::usage = "brSpecializeValue[e, sigma] is e with the rules sigma applied, or $Failed when a denominator vanishes.";
brSpecializeValue[e_, sigma_List] := Module[{t = Together[e], d},
  d = Denominator[t] /. sigma;
  If[Together[d] === 0, $Failed, Together[(Numerator[t] /. sigma)/d]]];

brGoodSpecializationQ::usage = "brGoodSpecializationQ[T, cv, rs, sigma] is True when sigma preserves the shape of the curve and of the residue data rs (see the header).";
brGoodSpecializationQ[T_Association, cv_Association, rs_Association, sigma_List, refc_ : None] := Module[
  {th = T["theta"], G = cv["G"], m = cv["m"], Gs, Ns, Es, Eps, Bps, ok = True, F0 = rs["F0"], Fk = rs["Fk"], B = rs["B"]},
  Gs = brSpecializeValue[G, sigma]; Ns = brSpecializeValue[rs["N"], sigma];
  Es = brSpecializeValue[rs["E"], sigma]; Eps = brSpecializeValue[rs["Ep"], sigma];
  Bps = brSpecializeValue[rs["Bp"], sigma];
  If[MemberQ[{Gs, Ns, Es, Eps, Bps}, $Failed], Return[False]];
  If[Exponent[Gs, th] =!= Exponent[G, th] || brSquarefreeProfile[Gs, th] =!= brSquarefreeProfile[G, th], Return[False]];
  If[Exponent[Bps, th] =!= Exponent[rs["Bp"], th] || Exponent[Es, th] =!= Exponent[rs["E"], th], Return[False]];
  If[Exponent[Bps, th] > 0 && (Exponent[tragerGCD[Bps, D[Bps, th]], th] > 0 ||
       Exponent[tragerGCD[Bps, Gs], th] > 0 || Resultant[Bps, Eps, th] === 0), Return[False]];
  (* a reference point theta = refc (orbit divisors) stays ordinary: unramified
     and off every residue place *)
  If[refc =!= None && (Together[Gs /. th -> refc] === 0 || Together[Es /. th -> refc] === 0), Return[False]];
  (* every residue place keeps its degree *)
  Do[
    If[brPlaceDegree[rs, a] =!= kDeg[F0, kGCD[F0, Bps, kc[F0, Ns^m Gs - a^m Eps^m], th], th], ok = False],
    {a, rs["roots"]}];
  ok];

brGoodSpecialization::usage = "brGoodSpecialization[T, cv, rs] is the first good specialization (brGoodSpecializationQ) among the deterministic candidates of growing height (brSpecializationPoint), at most brOption[\"SpecializationRounds\"] rounds of brOption[\"SpecializationTries\"]. Every goodness condition fails only on a thin set of points (Hilbert irreducibility for the place degrees, a proper hypersurface for the rest), so the search succeeds; exhausting it means a condition is false identically, an Internal error (never a certificate).";
brGoodSpecialization[T_Association, cv_Association, rs_Association, refc_ : None] := Module[{s},
  Do[
    s = brSpecializationPoint[rs["vars"], try];
    If[brGoodSpecializationQ[T, cv, rs, s, refc], Return[s, Module]],
    {try, 1, brOption["SpecializationTries"] brOption["SpecializationRounds"]}];
  throwInternal["NoGoodSpecialization", <|"vars" -> rs["vars"]|>]];

brSpecializedContext::usage = "brSpecializedContext[T, cv, rs, sigma] is the Trager component context of the specialized curve Y^m = G_sigma over Q with the residue data rewritten for d/dtheta (Nbar/Bp, see the header): keys \"B\", \"genus\", \"F\", \"residues\", \"N\", \"E\", \"Ep\", \"Bp\".";
brSpecializedContext[T_Association, cv_Association, rs_Association, sigma_List] := Module[
  {th = T["theta"], m = cv["m"], Gs, Ns, Eps, Bps, inv, Nbar, Nt},
  Gs = Expand[brSpecializeValue[cv["G"], sigma]];
  Ns = brSpecializeValue[rs["N"], sigma]; Eps = brSpecializeValue[rs["Ep"], sigma];
  Bps = Expand[brSpecializeValue[rs["Bp"], sigma]];
  inv = PolynomialExtendedGCD[Numerator[Together[Eps]], Bps, th];
  Nbar = Expand[PolynomialRemainder[Together[Ns D[Bps, th] Denominator[Together[Eps]] inv[[2, 1]]/inv[[1]]], Bps, th]];
  (* The stand-in remainder is Nt/(Bps Gs), Nt = Nbar Gs mod Bps: the same
     residues at the roots of Bps, its other poles on branch points (which
     Trager's frame excludes from its residue places), and it decays like
     theta^(-1 - deg G) at infinity.  The bare Nbar/Bps (decay theta^-1) was
     wrong: Trager's mod-p frame (jacFrame) moves theta = infinity to t = 0,
     and a pole there beyond the t-multiplicity of the frame's G becomes a
     spurious residue place, a divisor of nonzero degree whose "order" runs
     past the Weil bound -- Int 1/(x Log[x] (2 + Log[x]^4)^(2/3)) threw
     TorsionBeyondWeilBound (benchmark 257, cM-1-74).  With the factor Gs the
     pole order at t = 0 is at most k + 1 - deg G <= k m - deg G (k = Ceiling[
     deg G/m], m >= 2), so the frame keeps it out of its residue places. *)
  Nt = Expand[PolynomialRemainder[Nbar Gs, Bps, th]];
  <|"B" -> curveBasis[m, Gs, th], "genus" -> curveGenus[m, Gs, th], "F" -> rs["F0"],
    "residues" -> rs["roots"], "N" -> Nt, "E" -> Expand[Bps Gs], "Ep" -> Expand[D[Bps Gs, th]], "Bp" -> Bps|>];

brTorsionOrder::usage = "brTorsionOrder[T, cv, rs, b] is the only possible order of the class of delta = Sum_i b_i Delta(alpha_i): Trager's torsionOrderSpec on a good specialization (1 in genus 0). It throws NonElementary when the specialized divisor has infinite order (then so has delta).";
brTorsionOrder[T_Association, cv_Association, rs_Association, b_List] := Module[{sigma, ctx},
  sigma = brGoodSpecialization[T, cv, rs];
  ctx = brSpecializedContext[T, cv, rs, sigma];
  If[ctx["genus"] === 0, Return[1]];
  torsionOrderSpec[ctx, <|"res" -> rs["roots"], "b" -> b, "ref" -> None|>]];

brOrbitOrder::usage = "brOrbitOrder[T, cv, rs, spec, refc, nr] is the only possible order of the orbit divisor spec (m Delta(gamma) - kappa (theta - refc)_0, residue data rs over Q(gamma)) by specialization, as brTorsionOrder (1 in genus 0), or $Failed when the specialized divisor is shown of infinite order (orders modulo two primes disagree, or the order is refuted at a further good prime, torsionScreen). $Failed is NOT a certificate: an individual orbit divisor may have infinite order while every Q-basis combination is torsion; the caller falls back to the Q-basis construction.";
brOrbitOrder[T_Association, cv_Association, rs_Association, spec_Association, refc_, nr_Association] := Module[
  {sigma, ctx, r},
  sigma = brGoodSpecialization[T, cv, rs, refc];
  ctx = brSpecializedContext[T, cv, rs, sigma];
  If[ctx["genus"] === 0, Return[1]];
  ctx = Join[ctx, <|"refc" -> refc, "alpha" -> None, "R" -> nr["P"], "Z" -> nr["z"]|>];
  r = brCatch[Module[{n = torsionOrderSpec[ctx, spec]}, If[torsionScreen[ctx, spec, n], n, $Failed]]];
  Which[
    brFailQ[r] && First[brFailTag[r]] === "NonElementary", $Failed,
    brFailQ[r], Throw[r, "TragerFail"],
    True, r]];

(* -------------------------------------------------------------------------- *)
(* Torsion with the special places included (gap #9, LogPart:special-genus;  *)
(* docs/notes/bronstein_extensions.tex, E3).  The divisor lives on the       *)
(* shifted frame theta = a + 1/u of the genus-0 construction, where the      *)
(* special places (theta = infinity: u = 0; theta = 0: u = -1/a) are finite.  *)
(* A good specialization keeps the shape of the curve, of the finite residue  *)
(* places (brGoodSpecializationQ) AND of the special places: the frame values *)
(* of their Y-coordinates stay finite and nonzero, the multiplicity of G at   *)
(* theta = 0 is kept, and a branch-free denominator coprime to theta stays    *)
(* so.  The order of the specialized class is then Trager's torsionOrderWith  *)
(* with a builder that constructs the same divisor over GF(p).               *)
(* -------------------------------------------------------------------------- *)

brLowOrder::usage = "brLowOrder[p, th] is the multiplicity of th = 0 as a root of the polynomial p (0 when p(0) != 0).";
brLowOrder[p_, th_Symbol] := Module[{cl = CoefficientList[Expand[p], th]},
  If[AllTrue[cl, # === 0 &], Infinity, LengthWhile[cl, # === 0 &]]];

brSpecialLocationsGoodQ::usage = "brSpecialLocationsGoodQ[T, cv, rs, locs, sigma] is True when the specialization sigma keeps the special locations locs (brSpecialLocations) of the component: every location's frame coefficient u stays finite and nonzero, and at theta = 0 (exponential) the multiplicity of G and the non-vanishing of Bp are kept.";
brSpecialLocationsGoodQ[T_Association, cv_Association, rs_Association, locs_List, sigma_List] := Module[{th = T["theta"]},
  AllTrue[locs, With[{v = brSpecializeValue[#["u"], sigma]}, v =!= $Failed && v =!= 0] &] &&
  If[AnyTrue[locs, #["where"] === "zero" &],
    brLowOrder[brSpecializeValue[cv["G"], sigma], th] === brLowOrder[cv["G"], th] &&
    (brLowOrder[rs["Bp"], th] =!= 0 || brLowOrder[brSpecializeValue[rs["Bp"], sigma], th] === 0),
    True]];

brFrameDegreesKeptQ::usage = "brFrameDegreesKeptQ[fr, sigma] is True when the specialization sigma keeps the u-degree of every polynomial of the shifted frame fr (the curve Gu, the branch-free denominator Bph, the residue data Nh, Eh) and keeps Bph and Eh nonzero: no place moves to u = infinity (theta = a) and the divisor built from the specialized frame is the specialization of the divisor.";
brFrameDegreesKeptQ[fr_Association, sigma_List] := Module[{u = fr["u"]},
  AllTrue[{fr["Bu"]["G"], fr["Bph"], fr["Nh"], fr["Eh"]},
    Function[p, With[{ps = brSpecializeValue[p, sigma]},
      ps =!= $Failed && Together[ps] =!= 0 && Exponent[Expand[Numerator[Together[ps]]], u] === Exponent[Expand[Numerator[Together[p]]], u]]]]];

brFrameIdealDegree::usage = "brFrameIdealDegree[I, u, m] is the degree of the divisor of the fractional ideal I = M/den (Hermite form on an integral basis of rank m in the curve variable u): deg_u det M - m deg_u den. On the shifted frame a divisor avoiding u = infinity has this degree; a residue divisor's is 0.";
brFrameIdealDegree[I_Association, u_Symbol, m_Integer] :=
  Exponent[Expand[Numerator[Together[Det[I["M"]]]]], u] - Exponent[Expand[Denominator[Together[Det[I["M"]]]]], u] -
    m (Exponent[Expand[Numerator[Together[I["den"]]]], u] - Exponent[Expand[Denominator[Together[I["den"]]]], u]);

brSpecialDivisorIdeal::usage = "brSpecialDivisorIdeal[F, fr, finRoots, locs, rhoLists, cj] is the ideal over F, on the frame fr (brShiftFrame: \"Bu\", \"Nh\", \"Eh\", \"Bph\", \"u\", \"a\", \"k\"), of the divisor Sum_i cj_i Delta(alpha_i) + Sum_l cj_(nfin + l) P_l: the finite residue places of the residues finRoots and the special places P_l of the locations locs (in the order of brSpecialPlaceIdeals, one per root in rhoLists for an unramified location). Field-generic: the same construction over GF(p) builds the reduced divisor.";
brSpecialDivisorIdeal[F_Association, fr_Association, finRoots_List, locs_List, rhoLists_List, cj_List] := Module[
  {B = fr["Bu"], I, ideals, nfin = Length[finRoots]},
  ideals = Flatten[brSpecialPlaceIdeals[F, fr, #[[1]], #[[2]]] & /@ Transpose[{locs, rhoLists}], 1];
  I = idealUnit[B];
  Do[If[cj[[i]] =!= 0,
      I = idealTimes[F, B, I, idealPower[F, B, residueIdeal[F, B, finRoots[[i]], fr["Nh"], fr["Eh"], fr["Bph"]], cj[[i]]]]],
    {i, nfin}];
  Do[If[cj[[nfin + l]] =!= 0,
      I = idealTimes[F, B, I, idealPower[F, B, ideals[[l, 1]], cj[[nfin + l]]]]],
    {l, Length[ideals]}];
  I];

brSpecialTorsionOrder::usage = "brSpecialTorsionOrder[T, cv, rs, fr, locs, rhoLists, cj] is the only possible order of the class of the divisor of brSpecialDivisorIdeal (finite residues rs[\"roots\"] over rs[\"F0\"], special places of locs, coefficients cj) on the positive-genus curve of the component: specialize x and the lower monomials (brGoodSpecializationQ and brSpecialLocationsGoodQ), then Trager's torsionOrderWith on the specialized frame, with the special places' good-reduction conditions (frame values p-adic units, distinct roots rho, the point off the branch-free denominator, the multiplicity of the frame curve at the point kept) as its goodQ. It throws NonElementary when the orders modulo two primes disagree: the specialized divisor, and so the divisor, has infinite order.";
brSpecialTorsionOrder[T_Association, cv_Association, rs_Association, fr_Association, locs_List, rhoLists_List, cj_List] := Module[
  {sigma = None, s, m = cv["m"], u = fr["u"], F0 = rs["F0"], frs, locsS, g, pts, goodQ, ctx, mk, Gu0},
  Do[
    s = brSpecializationPoint[rs["vars"], try];
    (* the frame theta = a + 1/u must keep every u-degree: a residue place
       moved onto theta = a would leave the divisor (review, 2026-11-14:
       Bp_sigma(a) = 0 gave a divisor of nonzero degree and a false certificate) *)
    If[brGoodSpecializationQ[T, cv, rs, s, fr["a"]] && brSpecialLocationsGoodQ[T, cv, rs, locs, s] &&
       brFrameDegreesKeptQ[fr, s], sigma = s; Break[]],
    {try, 1, brOption["SpecializationTries"] brOption["SpecializationRounds"]}];
  If[sigma === None, throwInternal["NoGoodSpecialization", <|"vars" -> rs["vars"], "special" -> True|>]];
  Gu0 = Expand[brSpecializeValue[fr["Bu"]["G"], sigma]];
  frs = <|"u" -> u, "a" -> fr["a"], "k" -> fr["k"], "Bu" -> curveBasis[m, Gu0, u],
    "Nh" -> Expand[brSpecializeValue[fr["Nh"], sigma]], "Eh" -> Expand[brSpecializeValue[fr["Eh"], sigma]],
    "Bph" -> Expand[brSpecializeValue[fr["Bph"], sigma]]|>;
  locsS = Join[#, <|"u" -> brSpecializeValue[#["u"], sigma]|>] & /@ locs;
  g = curveGenus[m, Gu0, u];
  If[g === 0, Return[1]];
  (* Trager's reduced class arithmetic needs a root of Gu of multiplicity
     prime to m; with none in characteristic 0, no prime ever has one and the
     prime search would not end *)
  If[!AnyTrue[Select[tragerSquareFreeList[Gu0], !FreeQ[#[[1]], u] &], GCD[#[[2]], m] === 1 &],
    throwUnimplemented["TorsionModP:frame-root", "no root of the frame curve of multiplicity prime to m"]];
  pts = If[#["where"] === "inf", 0, -1/fr["a"]] & /@ locsS;
  (* the goodQ conditions must hold in characteristic 0, else no prime meets
     them and the (unbounded, Chebotarev) prime search would never end *)
  Do[
    If[Together[frs["Bph"] /. u -> pts[[l]]] === 0 ||
       (!locsS[[l]]["ramified"] && (Together[Gu0 /. u -> pts[[l]]] === 0 ||
          !DuplicateFreeQ[kc[F0, #] & /@ rhoLists[[l]]] || MemberQ[kc[F0, #] & /@ rhoLists[[l]], 0])),
      throwInternal["brSpecialPlaceDegenerate", <|"location" -> locsS[[l]]|>]],
    {l, Length[locsS]}];
  goodQ = Function[{p, rho}, Module[{red = Function[e, reduceResidue[F0, e, p, rho]]},
    AllTrue[Range[Length[locsS]], Function[l,
      denominatorsPrimeToQ[{locsS[[l]]["u"], pts[[l]], locsS[[l]]["kappa"], rhoLists[[l]]}, p] && red[locsS[[l]]["u"]] =!= 0 &&
      red[frs["Bph"] /. u -> pts[[l]]] =!= 0 &&
      If[locsS[[l]]["ramified"],
        rootMultiplicity[fieldGF[p], gfMod[Gu0, p], gfMod[pts[[l]], p], u] === rootMultiplicity[F0, Gu0, pts[[l]], u],
        red[Gu0 /. u -> pts[[l]]] =!= 0 &&
          With[{rr = red /@ rhoLists[[l]]}, !MemberQ[rr, 0] && DuplicateFreeQ[rr]]]]]]];
  (* a degree-0 divisor stays one: any mismatch between the specialized data
     and the divisor (a place lost or moved) is an Internal error here, never
     a certificate from the order test below *)
  With[{d = brFrameIdealDegree[brSpecialDivisorIdeal[F0, frs, rs["roots"], locsS, rhoLists, cj], u, m]},
    If[d =!= 0, throwInternal["brSpecialDivisorDegree", <|"degree" -> d, "sigma" -> sigma|>]]];
  ctx = <|"B" -> frs["Bu"], "genus" -> g, "F" -> F0, "residues" -> rs["roots"],
    "N" -> frs["Nh"], "Ep" -> frs["Eh"], "Bp" -> frs["Bph"], "goodQ" -> goodQ|>;
  mk = Function[{Fx, Bx, red},
    brSpecialDivisorIdeal[Fx,
      Join[frs, <|"Bu" -> Bx, "Nh" -> red[frs["Nh"]], "Eh" -> red[frs["Eh"]], "Bph" -> red[frs["Bph"]]|>],
      red /@ rs["roots"],
      Map[Join[#, <|"u" -> red[#["u"]], "kappa" -> If[#["kappa"] === None, None, red[#["kappa"]]]|>] &, locsS],
      Map[red, rhoLists, {2}], cj]];
  torsionOrderWith[ctx, mk]];
