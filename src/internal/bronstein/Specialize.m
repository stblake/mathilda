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

brSpecializationPoint::usage = "brSpecializationPoint[vars, try] is the deterministic try-th candidate specialization: a list of rules v -> rational for the symbols vars.";
brSpecializationPoint[vars_List, try_Integer] := MapIndexed[
  #1 -> (Mod[7919 First[#2]^2 + 104729 try + 31 First[#2] try^2, brOption["SpecializationHeight"]] + 2)/
        (Mod[13 try + 5 First[#2], 7] + 1) &, vars];

brSquarefreeProfile::usage = "brSquarefreeProfile[p, th] is the sorted list {{deg_th, multiplicity}, ...} of the squarefree factors of p involving th.";
brSquarefreeProfile[p_, th_Symbol] := Sort[{Exponent[#[[1]], th], #[[2]]} & /@
  Select[FactorSquareFreeList[Expand[p]], !FreeQ[#[[1]], th] &]];

brSpecializeValue::usage = "brSpecializeValue[e, sigma] is e with the rules sigma applied, or $Failed when a denominator vanishes.";
brSpecializeValue[e_, sigma_List] := Module[{t = Together[e], d},
  d = Denominator[t] /. sigma;
  If[Together[d] === 0, $Failed, Together[(Numerator[t] /. sigma)/d]]];

brGoodSpecializationQ::usage = "brGoodSpecializationQ[T, cv, rs, sigma] is True when sigma preserves the shape of the curve and of the residue data rs (see the header).";
brGoodSpecializationQ[T_Association, cv_Association, rs_Association, sigma_List] := Module[
  {th = T["theta"], G = cv["G"], m = cv["m"], Gs, Ns, Es, Eps, Bps, ok = True, F0 = rs["F0"], Fk = rs["Fk"], B = rs["B"]},
  Gs = brSpecializeValue[G, sigma]; Ns = brSpecializeValue[rs["N"], sigma];
  Es = brSpecializeValue[rs["E"], sigma]; Eps = brSpecializeValue[rs["Ep"], sigma];
  Bps = brSpecializeValue[rs["Bp"], sigma];
  If[MemberQ[{Gs, Ns, Es, Eps, Bps}, $Failed], Return[False]];
  If[Exponent[Gs, th] =!= Exponent[G, th] || brSquarefreeProfile[Gs, th] =!= brSquarefreeProfile[G, th], Return[False]];
  If[Exponent[Bps, th] =!= Exponent[rs["Bp"], th] || Exponent[Es, th] =!= Exponent[rs["E"], th], Return[False]];
  If[Exponent[Bps, th] > 0 && (Exponent[PolynomialGCD[Bps, D[Bps, th]], th] > 0 ||
       Exponent[PolynomialGCD[Bps, Gs], th] > 0 || Resultant[Bps, Eps, th] === 0), Return[False]];
  (* every residue place keeps its degree *)
  Do[
    If[brPlaceDegree[rs, a] =!= kDeg[F0, kGCD[F0, Bps, kc[F0, Ns^m Gs - a^m Eps^m], th], th], ok = False],
    {a, rs["roots"]}];
  ok];

brGoodSpecialization::usage = "brGoodSpecialization[T, cv, rs] is the first good specialization (brGoodSpecializationQ) among brOption[\"SpecializationTries\"] deterministic candidates; it throws a Decline when none is found (never a certificate).";
brGoodSpecialization[T_Association, cv_Association, rs_Association] := Module[{s},
  Do[
    s = brSpecializationPoint[rs["vars"], try];
    If[brGoodSpecializationQ[T, cv, rs, s], Return[s, Module]],
    {try, 1, brOption["SpecializationTries"]}];
  throwDecline["no good specialization found"]];

brSpecializedContext::usage = "brSpecializedContext[T, cv, rs, sigma] is the Trager component context of the specialized curve Y^m = G_sigma over Q with the residue data rewritten for d/dtheta (Nbar/Bp, see the header): keys \"B\", \"genus\", \"F\", \"residues\", \"N\", \"E\", \"Ep\", \"Bp\".";
brSpecializedContext[T_Association, cv_Association, rs_Association, sigma_List] := Module[
  {th = T["theta"], m = cv["m"], Gs, Ns, Eps, Bps, inv, Nbar},
  Gs = Expand[brSpecializeValue[cv["G"], sigma]];
  Ns = brSpecializeValue[rs["N"], sigma]; Eps = brSpecializeValue[rs["Ep"], sigma];
  Bps = Expand[brSpecializeValue[rs["Bp"], sigma]];
  inv = PolynomialExtendedGCD[Numerator[Together[Eps]], Bps, th];
  Nbar = Expand[PolynomialRemainder[Together[Ns D[Bps, th] Denominator[Together[Eps]] inv[[2, 1]]/inv[[1]]], Bps, th]];
  <|"B" -> curveBasis[m, Gs, th], "genus" -> curveGenus[m, Gs, th], "F" -> rs["F0"],
    "residues" -> rs["roots"], "N" -> Nbar, "E" -> Bps, "Ep" -> Expand[D[Bps, th]], "Bp" -> Bps|>];

brTorsionOrder::usage = "brTorsionOrder[T, cv, rs, b] is the only possible order of the class of delta = Sum_i b_i Delta(alpha_i): Trager's torsionOrderSpec on a good specialization (1 in genus 0). It throws NonElementary when the specialized divisor has infinite order (then so has delta).";
brTorsionOrder[T_Association, cv_Association, rs_Association, b_List] := Module[{sigma, ctx},
  sigma = brGoodSpecialization[T, cv, rs];
  ctx = brSpecializedContext[T, cv, rs, sigma];
  If[ctx["genus"] === 0, Return[1]];
  torsionOrderSpec[ctx, <|"res" -> rs["roots"], "b" -> b, "ref" -> None|>]];
