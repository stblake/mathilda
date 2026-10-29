(* ::Package:: *)

(* ========================================================================== *)
(* TragerLogTerms.m -- the logarithmic part by Trager's construction,          *)
(*                     [T] Ch 5 Sec 1-3 and Ch 6; [B] Sec 2.5.                 *)
(*                                                                            *)
(* Remainder r = (N/E) Y dx, residues alpha_i in K with alpha_i =              *)
(* Sum_j b_ij beta_j (b_ij integers).  For each j the divisor                  *)
(*     delta_j = Sum_i b_ij Delta(alpha_i)                                     *)
(* has degree 0 (the residues sum to zero) and no place at infinity.  If      *)
(* N_j delta_j = div(u_j), then Sum_j beta_j/N_j du_j/u_j has residue          *)
(* Sum_j beta_j b_ij = alpha_i at every place with residue alpha_i, so        *)
(*     omega = r - Sum_j (beta_j/N_j) du_j/u_j                                 *)
(* has no residues and at most simple poles: it is holomorphic.  It is zero   *)
(* iff the integral is elementary ([T] p.71, [B] p.19), and then              *)
(*     Int r = Sum_j (beta_j/N_j) Log[u_j].                                    *)
(* ========================================================================== *)

tragerLogTerms::usage = "tragerLogTerms[ctx] constructs the logarithmic part of the component ctx (keys \"B\", \"F\", \"genus\", \"residues\", \"basis\", \"coeffs\", \"N\", \"E\", \"Ep\", \"Bp\"). It returns {{c, u}, ...} with c in K and u an element of K(x)(Y) (coordinate list) such that Sum c du/u has the residues of the remainder. For each basis divisor delta_j the order N_j comes from torsionOrderSpec and the generator of N_j delta_j from generatorTerms -- possibly as a short product Prod h^e, each factor its own term. Throws a NonElementary certificate when a residue divisor has infinite order.";
tragerLogTerms[ctx0_Association] := Module[
  {ctx = withFrameAlpha[ctx0], F = ctx0["F"], b = ctx0["coeffs"], out = {}, spec, Nj, terms},
  Do[
    spec = <|"res" -> ctx["residues"], "b" -> b[[All, j]], "ref" -> None|>;
    Nj = torsionOrderSpec[ctx, spec];
    terms = generatorTerms[ctx, spec, Nj];
    If[terms === $Failed,
      If[ctx["genus"] === 0, throwInternal["Genus0DivisorNotPrincipal"]];
      throwNonElementary["residue divisor is not a torsion point of the Jacobian",
        <|"CandidateOrder" -> Nj|>]];
    out = Join[out, {kc[F, ctx["basis"][[j]] #[[1]]/Nj], #[[2]]} & /@ terms],
    {j, Length[ctx["basis"]]}];
  out];

withFrameAlpha::usage = "withFrameAlpha[ctx] adds ctx[\"alpha\"]: a rational root of G of multiplicity prime to m, the centre of the characteristic-0 reduced frame (Jacobian.m), or None when G has none.";
withFrameAlpha[ctx_Association] := If[KeyExistsQ[ctx, "alpha"], ctx,
  Join[ctx, <|"alpha" -> First[Join[frameRoots[fieldQ[], ctx["B"]["G"], ctx["B"]["m"], ctx["B"]["x"]], {None}]]|>]];

generatorTerms::usage = "generatorTerms[ctx, spec, Nn] is a generator of Nn delta (delta the divisor spec, Nn >= 1) as a list {{e, u}, ...}, u elements of the working frame (coordinate lists), with div(Prod u^e) = Nn delta; $Failed when Nn delta is not principal. For Nn up to tragerOption[\"DirectPowerMaxOrder\"], or without a characteristic-0 frame root, the ideal of Nn delta is formed and tested directly ([T] Ch 6 Sec 1); for larger Nn with a frame root (ctx[\"alpha\"]) Nn delta is formed by doubling and adding in reduced form (Jacobian.m) and the generator is the product of the principal factors met.";
generatorTerms[ctx_Association, spec_Association, Nn_Integer] := Module[
  {F = ctx["F"], B = ctx["B"], x, fr, Jn, c1, cN, I, u, Yt = Unique["Yt"], Yw = Unique["Yw"]},
  x = B["x"];
  (* In characteristic 0 the direct ideal power (Hermite form modulo the     *)
  (* determinant) is cheap for moderate Nn, while reduced arithmetic pays    *)
  (* norms and adjugates over the residue field at every step: use the       *)
  (* reduced classes (and the compact generator) only for large orders.      *)
  If[ctx["alpha"] =!= None && Nn > tragerOption["DirectPowerMaxOrder"],
    fr = jacFrame[F, B["m"], B["G"], rc[F, ctx["N"]/ctx["E"]], x, ctx["alpha"]];
    Jn = specFrameIdeal[F, fr, spec, Identity];
    c1 = jacClass[F, fr["B"], Jn[[1]], Jn[[2]], True];
    cN = If[Nn === 1, c1, jacMultiple[F, fr["B"], c1, Nn, True]];
    If[!jacZeroQ[cN], Return[$Failed]];
    (* Nn delta + div(Prod h^e) = 0, so Prod h^(-e) generates Nn delta *)
    Return[combineGenerator[F, B,
      {-#[[2]], afFromY[F, B, Together[frameToX[fr, afToY[fr["B"], #[[1]], Yt], Yt, Yw]], Yw]} & /@ cN["H"]]]];
  (* Nn delta built directly from its residue coefficients (Nn b_i): one     *)
  (* Hensel-lifted two-generator product per residue, no repeated squaring  *)
  I = specIdealX[ctx, scaleSpec[spec, Nn]];
  u = principalGenerator[F, B, I];
  If[u === None, $Failed, {{1, u}}]];

logSumResidual::usage = "logSumResidual[F, B, r, logs] is r - Sum c_j u_j'/u_j for log terms {{c_j, u_j}} (u_j coordinate lists, possibly rational), each logarithmic derivative formed with the inverse of u_j (rationalElementInverse); zero exactly when r = Sum c_j du_j/u_j. Unlike logPartResidual it never multiplies the u_j together, which matters for compact generators with many factors.";
logSumResidual[F_Association, B_Association, r_List, logs_List] := Module[{acc = r},
  Do[
    acc = afMinus[F, acc, afScale[F, lg[[1]],
      afTimes[F, B, afD[F, B, lg[[2]]], rationalElementInverse[F, B, lg[[2]]]]]],
    {lg, logs}];
  acc];

logPartResidual::usage = "logPartResidual[F, B, r, logs] is the element P (r - Sum c_j u_j'/u_j) with P = Prod_j u_j, computed without division; it is zero exactly when r = Sum c_j du_j/u_j.";
logPartResidual[F_Association, B_Association, r_List, logs_List] := Module[{us, P, acc},
  us = Last /@ logs;
  P = Fold[afTimes[F, B, #1, #2] &, afOne[B], us];
  acc = afTimes[F, B, r, P];
  Do[
    acc = afMinus[F, acc,
      afScale[F, logs[[j, 1]],
        afTimes[F, B, afD[F, B, us[[j]]],
          Fold[afTimes[F, B, #1, #2] &, afOne[B], Delete[us, j]]]]],
    {j, Length[us]}];
  acc];

(* -------------------------------------------------------------------------- *)
(* The logarithmic part one Galois orbit at a time.                           *)
(*                                                                            *)
(* The Q-basis construction above needs the field generated by ALL residues,  *)
(* whose degree can be the product of the factorials of the factor degrees of *)
(* R.  For a root gamma of an irreducible factor q of R, the divisor          *)
(*     delta_gamma = m Delta(gamma) - kappa (x - c)_0,                         *)
(* (kappa = deg Delta(gamma), (x - c)_0 the m places over an ordinary point   *)
(* x = c, a divisor over Q) has degree 0 and is defined over K = Q[Z]/(q), of  *)
(* degree deg q only.  If N delta_gamma = div(u_gamma), then                  *)
(*     Sum_q Sum_{q(gamma) = 0} gamma/(m N_q) Log[u_gamma]                     *)
(*       = Sum_q 1/(m N_q) RootSum[q, # Log[u(#)] &]                           *)
(* has residue gamma at the places with residue gamma, and at the places over *)
(* x = c the residue -(1/m) Sum gamma kappa_gamma = 0 (all residues together   *)
(* sum to zero).  On a curve of genus 0 every delta_gamma is principal, so    *)
(* this succeeds; in higher genus it may fail because an individual           *)
(* delta_gamma has infinite order while the Q-basis combinations are torsion. *)
(* That failure is NOT a certificate: the caller falls back to the complete   *)
(* Q-basis construction.                                                      *)
(* -------------------------------------------------------------------------- *)

referencePoint::usage = "referencePoint[ctx] is a rational c with G(c) != 0 and E(c) != 0 in the working variable: the places over x = c are ordinary (unramified, outside the support of every residue divisor).";
referencePoint[ctx_Association] := Module[{B = ctx["B"], x, cands, c},
  x = B["x"];
  cands = Join[{0}, Flatten[Table[{i, -i}, {i, tragerOption["ShiftBound"]}]]];
  (* (SelectFirst's default argument is evaluated eagerly: test afterwards) *)
  c = SelectFirst[cands, (B["G"] /. x -> #) =!= 0 && (ctx["E"] /. x -> #) =!= 0 &, None];
  If[c === None, throwInternal["NoReferencePoint"], c]];

idealDegree::usage = "idealDegree[F, B, I] is the degree of the divisor of the integral ideal I (den = 1): the degree in x of the determinant of its Hermite form.";
idealDegree[F_Association, B_Association, I_Association] :=
  Exponent[kc[F, Times @@ Table[I["M"][[i, i]], {i, B["m"]}]], B["x"]];

scaleSpec::usage = "scaleSpec[spec, k] is the divisor spec of k times the divisor spec.";
scaleSpec[spec_Association, k_Integer] := <|"res" -> spec["res"], "b" -> k spec["b"],
  "ref" -> Which[spec["ref"] === None, None,
                 spec["ref"][[1]] === "Pinf", {"Pinf", k spec["ref"][[2]]},
                 True, {"R0", spec["ref"][[2]], k spec["ref"][[3]]}]|>;

combineGenerator::usage = "combineGenerator[F, B, terms] multiplies a compact generator {{e, u}, ...} (div(Prod u^e) = N delta) out into a single element {{1, U}} when Total[Abs[e]] <= 6, which gives one logarithm in the answer; a larger product (a large torsion order) is left as the list -- a sum of logarithms of small elements instead of one huge one.";
combineGenerator[F_Association, B_Association, terms_List] := Module[{U = afOne[B], inv},
  If[Total[Abs[First /@ terms]] > 6, Return[terms]];
  Do[
    If[e > 0,
      U = afTimes[F, B, U, afPower[F, B, u, e]],
      inv = rationalElementInverse[F, B, u];
      U = afTimes[F, B, U, afPower[F, B, inv, -e]]],
    {t0, terms}, {e, {t0[[1]]}}, {u, {t0[[2]]}}];
  {{1, U}}];

rationalElementInverse::usage = "rationalElementInverse[F, B, u] is 1/u for an element with rational coordinates over F: u = P/d with P polynomial and d over the prime field, so 1/u = d adj(P)/Norm(P) (elementInverse).";
rationalElementInverse[F_Association, B_Association, u_List] := Module[{d, P, inv},
  d = Fold[PolynomialLCM, 1, Denominator[rc[F, #]] & /@ u];
  P = kc[F, Together[# d]] & /@ u;
  inv = elementInverse[F, B, P];
  rc[F, # d/inv[[2]]] & /@ inv[[1]]];

tragerOrbitLogTerms::usage = "tragerOrbitLogTerms[ctx] is the logarithmic part one Galois orbit of residues at a time (see the comment above): {{q, k, p}, ...} in the format of kauersLogTerms (contribution (1/k) Sum_{q(gamma)=0} gamma Log[p(gamma)]), or $Failed when some orbit divisor has no principal multiple (the caller then uses the Q-basis construction). With a characteristic-0 frame root the orbit divisor is Delta(gamma) - kappa P (P the place over the root, of degree 1); otherwise m Delta(gamma) - kappa (x - c)_0.";
tragerOrbitLogTerms[ctx0_Association] := Module[
  {ctx = withFrameAlpha[ctx0], B = ctx0["B"], Z = ctx0["Z"], c, facs, out = {}, failed = False},
  c = referencePoint[ctx];
  facs = First /@ Select[FactorList[ctx["R"]], !FreeQ[#[[1]], Z] &];
  Do[
    If[!failed,
      Module[{qm = Expand[q/Coefficient[q, Z, Exponent[q, Z]]], K, gamma, octx, I0, kappa, spec, scale, Nq, r},
        If[Exponent[qm, Z] === 1,
          K = fieldQ[]; gamma = -Coefficient[qm, Z, 0],
          K = fieldNF[Z, qm, None]; gamma = Z];
        I0 = residueIdeal[K, B, gamma, ctx["N"], ctx["Ep"], ctx["Bp"]];
        kappa = idealDegree[K, B, I0];
        If[ctx["alpha"] =!= None,
          spec = <|"res" -> {gamma}, "b" -> {1}, "ref" -> {"Pinf", -kappa}|>; scale = 1,
          spec = <|"res" -> {gamma}, "b" -> {B["m"]}, "ref" -> {"R0", c, -kappa}|>; scale = B["m"]];
        octx = Join[ctx, <|"F" -> K, "residues" -> {gamma}, "refc" -> If[ctx["alpha"] === None, c, None]|>];
        r = catchTrager[
          Nq = torsionOrderSpec[octx, spec];
          generatorTerms[octx, spec, Nq]];
        Which[
          tragerFailQ[r] && r[[1]] === "NonElementary", failed = True,
          tragerFailQ[r], Throw[r, "TragerFail"],
          r === $Failed, failed = True,
          True,
            (* gamma/(scale Nq) Log[Prod u^e]: one term per factor *)
            out = Join[out, {qm, scale Nq/#[[1]], Together[afToY[B, #[[2]], ctx["Ysym"]]]} & /@ r]]]],
    {q, facs}];
  If[failed, $Failed, out]];
