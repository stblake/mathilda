(* ::Package:: *)

(* ========================================================================== *)
(* TorsionModP.m -- the order of a divisor class by reduction modulo two      *)
(*                  primes, [T] Ch 6 Sec 2-3 pp.66-70; [B] p.19.               *)
(*                                                                            *)
(* On a curve of positive genus a degree-0 divisor delta need not be          *)
(* principal; the integral is elementary only if N delta is principal for    *)
(* some N >= 1 (delta is a torsion point of the Jacobian).  Trager bounds N   *)
(* by reduction modulo a prime p of good reduction:                           *)
(*                                                                            *)
(*   * over GF(p) the divisor class group of degree 0 is finite, of order at  *)
(*     most (1 + Sqrt[p])^(2 genus) (Weil), so the order o_p of the reduced    *)
(*     class is found by trying M = 1, 2, ...;                                *)
(*   * reduction is injective on classes of order prime to p ([T] Theorem     *)
(*     p.67): if delta has order N = n0 p^a q^b (n0 prime to p q) then         *)
(*     o_p = n0 q^b p^a' and o_q = n0 p^a q^b' with a' <= a, b' <= b.          *)
(*                                                                            *)
(* Hence, with two good primes p and q: if the prime-to-pq parts of o_p and   *)
(* o_q differ, delta has infinite order and the integral is not elementary;   *)
(* otherwise the only possible order is N = n0 p^v_p(o_q) q^v_q(o_p), and ONE *)
(* principality test of N delta over the residue field decides.               *)
(*                                                                            *)
(* Good reduction ([T] pp.66-69), checked here for the data of one component  *)
(* (curve Y^m = G, divisors built from N, E', Bp and the residues):           *)
(*   p does not divide m or any denominator; the squarefree part of G stays   *)
(*   squarefree of the same degree (so the multiplicity pattern, the integral *)
(*   basis and the genus are preserved); Bp stays squarefree and prime to G   *)
(*   and E'; the minimal polynomial of the residue-field generator has a       *)
(*   simple root rho modulo p (so K maps to GF(p)); the residues map to        *)
(*   distinct nonzero elements of GF(p).                                      *)
(* ========================================================================== *)

pPart::usage = "pPart[n, p] is the largest power of the prime p dividing the positive integer n.";
pPart[n_Integer, p_Integer] := p^IntegerExponent[n, p];

denominatorsPrimeToQ::usage = "denominatorsPrimeToQ[expr, p] is True when no rational number occurring in expr has a denominator divisible by p.";
denominatorsPrimeToQ[expr_, p_Integer] :=
  AllTrue[Cases[{expr}, r_Rational :> Denominator[r], {0, Infinity}], Mod[#, p] =!= 0 &];

nonzeroModQ::usage = "nonzeroModQ[c, p] is True when the rational number c is nonzero modulo p (p-integral by assumption).";
nonzeroModQ[c_, p_Integer] := gfMod[c, p] =!= 0;

goodPrime::usage = "goodPrime[ctx, p] tests the prime p for good reduction of the component data ctx (keys \"B\", \"genus\", \"F\", \"residues\", \"N\", \"Ep\", \"Bp\"). It returns <|\"p\" -> p, \"rho\" -> rho, \"a\" -> a|> (rho the image of the residue-field generator, None over Q; a the centre of the reduced frame mod p: the image of ctx[\"alpha\"] when given, else a root of G mod p of multiplicity prime to m) or None.";
goodPrime[ctx_Association, p_Integer] := Module[
  {B = ctx["B"], x, m, G, Gs, Bp = ctx["Bp"], Ep = ctx["Ep"], F = ctx["F"], rho = None, lin, alphas, a},
  x = B["x"]; m = B["m"]; G = B["G"];
  If[Mod[m, p] === 0, Return[None]];
  If[!denominatorsPrimeToQ[{G, ctx["N"], Ep, Bp, ctx["residues"], If[F["Type"] === "NF", F["mu"], 0]}, p], Return[None]];
  Gs = Together[G/PolynomialGCD[G, D[G, x]]];
  If[!nonzeroModQ[Coefficient[G, x, Exponent[G, x]], p], Return[None]];
  If[!nonzeroModQ[Coefficient[Bp, x, Exponent[Bp, x]], p], Return[None]];
  If[!nonzeroModQ[Discriminant[Expand[Gs Bp], x], p], Return[None]];
  If[!nonzeroModQ[Resultant[Bp, Ep, x], p], Return[None]];
  If[curveGenus[m, G, x, p] =!= ctx["genus"], Return[None]];
  If[F["Type"] === "NF",
    If[!nonzeroModQ[Discriminant[F["mu"], F["t"]], p], Return[None]];
    lin = Select[FactorList[F["mu"], Modulus -> p], Exponent[#[[1]], F["t"]] === 1 &];
    If[lin === {}, Return[None]];
    rho = gfMod[-Coefficient[lin[[1, 1]], F["t"], 0] PowerMod[Coefficient[lin[[1, 1]], F["t"], 1], -1, p], p]];
  alphas = reduceResidue[F, #, p, rho] & /@ ctx["residues"];
  If[MemberQ[alphas, 0] || !DuplicateFreeQ[alphas], Return[None]];
  (* With the residue polynomial at hand, require it to stay squarefree with  *)
  (* no zero root mod p: then NO two residues (of any factor) collide, which  *)
  (* the per-orbit construction needs (its ctx lists only its own residue).   *)
  If[KeyExistsQ[ctx, "R"] && ctx["R"] =!= 1,
    If[!denominatorsPrimeToQ[ctx["R"], p] ||
       !nonzeroModQ[Discriminant[ctx["R"], ctx["Z"]], p] ||
       !nonzeroModQ[ctx["R"] /. ctx["Z"] -> 0, p], Return[None]]];
  (* The reduced class arithmetic (Jacobian.m) needs a root a of G mod p of   *)
  (* multiplicity prime to m -- the image of the characteristic-0 frame root  *)
  (* when there is one, so that P reduces to P.  The ordinary point c of an   *)
  (* R0 reference must stay ordinary.                                         *)
  a = If[KeyExistsQ[ctx, "alpha"] && ctx["alpha"] =!= None,
    With[{ab = gfMod[ctx["alpha"], p]},
      If[GCD[rootMultiplicity[fieldGF[p], G, ab, x], m] === 1, ab, None]],
    First[Join[frameRoots[fieldGF[p], G, m, x], {None}]]];
  If[a === None, Return[None]];
  If[KeyExistsQ[ctx, "refc"] && ctx["refc"] =!= None,
    If[!nonzeroModQ[G /. x -> ctx["refc"], p] || !nonzeroModQ[ctx["E"] /. x -> ctx["refc"], p] ||
       gfMod[a - ctx["refc"], p] === 0, Return[None]]];
  <|"p" -> p, "rho" -> rho, "a" -> a|>];

reduceResidue::usage = "reduceResidue[F, a, p, rho] maps the element a of the residue field F to GF(p), sending the generator t to rho.";
reduceResidue[F_Association, a_, p_Integer, rho_] :=
  gfMod[If[F["Type"] === "NF", kc[F, a] /. F["t"] -> rho, a], p];

orderModP::usage = "orderModP[ctx, b, gp] is the order of the class of the divisor Sum_i b_i Delta(alpha_i) of the component ctx after reduction at the good prime gp = <|\"p\" -> p, \"rho\" -> rho|>: the least M >= 1 with M delta principal over GF(p). An internal error past the Weil bound (1 + Sqrt[p])^(2 genus), which bounds the order (no other limit).";
orderModP[ctx_Association, b_List, gp_Association] :=
  orderModPSpec[ctx, <|"res" -> ctx["residues"], "b" -> b, "ref" -> None|>, withFrameRoot[ctx, gp]];

withFrameRoot::usage = "withFrameRoot[ctx, gp] adds a frame centre to a good-prime record that lacks one.";
withFrameRoot[ctx_Association, gp_Association] := Module[{rts},
  If[KeyExistsQ[gp, "a"] && gp["a"] =!= None, Return[gp]];
  rts = frameRoots[fieldGF[gp["p"]], ctx["B"]["G"], ctx["B"]["m"], ctx["B"]["x"]];
  If[rts === {}, throwDecline["no frame root modulo " <> ToString[gp["p"]]]];
  Join[gp, <|"a" -> First[rts]|>]];

residueCombinationBuilder::usage = "residueCombinationBuilder[ctx, b] is the ideal builder (see orderModPWith) of the divisor Sum_i b_i Delta(alpha_i): a function of (F, B, red) returning the product of the residue ideals over F, every residue-field element and polynomial passed through red first.";
residueCombinationBuilder[ctx_Association, b_List] := Function[{Fx, Bx, red},
  divisorIdeal[Fx, Bx, red /@ ctx["residues"], b, red[ctx["N"]], red[ctx["Ep"]], red[ctx["Bp"]]]];

orderModPWith::usage = "orderModPWith[ctx, mk, gp] is the order over GF(p) of the divisor class whose ideal mk[F, B, red] builds: mk is called with the field GF(p), the reduced basis, and the reduction map red (residue-field generator -> rho, coefficients mod p). Same (theorem-only) bound as orderModP.";
orderModPWith[ctx_Association, mk_, gp_Association] := Module[
  {p = gp["p"], Fp, Bm, red, I, J, weil, M = 1},
  Fp = fieldGF[p]; Bm = basisMod[ctx["B"], p];
  red = Function[e, reduceResidue[ctx["F"], e, p, gp["rho"]]];
  I = mk[Fp, Bm, red];
  weil = Floor[(1 + Sqrt[N[p]])^(2 ctx["genus"])] + 1;
  J = I;
  While[!principalQ[Fp, Bm, J],
    M++;
    If[M > weil, throwInternal["TorsionBeyondWeilBound", <|"p" -> p, "Bound" -> weil|>]];
    J = idealTimes[Fp, Bm, J, I]];
  M];

goodPrimes::usage = "goodPrimes[ctx, k] is a list of k distinct good primes (as goodPrime records) for the component ctx, smallest first. The search is unbounded: a positive density of primes is good (Chebotarev).";
goodPrimes[ctx_Association, k_Integer] := Module[{found = {}, p = 2, gp},
  (* All but finitely many primes have good reduction except for the root   *)
  (* condition on the residue field, which a positive density of primes     *)
  (* meets (Chebotarev): the search ends, and needs no bound.               *)
  While[Length[found] < k,
    p = NextPrime[p];
    gp = goodPrime[ctx, p];
    If[gp =!= None, AppendTo[found, gp]]];
  found];

torsionOrder::usage = "torsionOrder[ctx, b] is the only possible order N of the divisor class delta = Sum_i b_i Delta(alpha_i) of the component ctx, from its orders modulo two good primes ([T] Ch 6 Sec 3). It throws a NonElementary certificate when the prime-to-pq parts of the two orders differ (delta has infinite order). The caller must still test N delta for principality over the residue field: if that fails, delta has infinite order too.";
torsionOrder[ctx_Association, b_List] := torsionOrderSpec[ctx, <|"res" -> ctx["residues"], "b" -> b, "ref" -> None|>];

torsionOrderWith::usage = "torsionOrderWith[ctx, mk] is torsionOrder for the divisor whose ideal the builder mk constructs (see orderModPWith).";
torsionOrderWith[ctx_Association, mk_] := Module[{gps, p, q, op, oq, n0p, n0q},
  If[ctx["genus"] === 0, Return[1]];
  gps = goodPrimes[ctx, 2];
  {p, q} = #["p"] & /@ gps;
  op = orderModPWith[ctx, mk, gps[[1]]];
  oq = orderModPWith[ctx, mk, gps[[2]]];
  n0p = op/(pPart[op, p] pPart[op, q]);
  n0q = oq/(pPart[oq, p] pPart[oq, q]);
  If[n0p =!= n0q,
    throwNonElementary["residue divisor of infinite order: its orders modulo " <>
      ToString[p] <> " and " <> ToString[q] <> " disagree",
      <|"Orders" -> {op, oq}, "Primes" -> {p, q}|>]];
  n0p pPart[oq, p] pPart[op, q]];

(* -------------------------------------------------------------------------- *)
(* Divisor specs and the reduced-class order search.                          *)
(*                                                                            *)
(* A divisor of the component is described independently of the frame by a  *)
(* spec <|"res" -> {alpha_i}, "b" -> {b_i}, "ref" -> ref|>:                   *)
(*     Sum_i b_i Delta(alpha_i) + (reference term),                           *)
(* with ref None, {"Pinf", n} (n times the place over the characteristic-0    *)
(* frame root ctx["alpha"]) or {"R0", c, n} (n times the m places over the    *)
(* ordinary point x = c).  specFrameIdeal builds J = I(-D_fin) and the        *)
(* coefficient of P in a reduced frame; specIdealX builds I(D) in the working *)
(* frame (no Pinf term).                                                      *)
(* -------------------------------------------------------------------------- *)

specFrameIdeal::usage = "specFrameIdeal[F, fr, spec, red] is {J, n} for the divisor spec in the reduced frame fr over F: J = I(-D_fin) (the residues and the reference point passed through red) and n the coefficient of the frame's place P (the Pinf term). An R0 term at x = c is the divisor of the polynomial (a - c) t + 1 in the frame.";
specFrameIdeal[F_Association, fr_Association, spec_Association, red_] := Module[
  {Bt = fr["B"], t = fr["t"], J, n = 0, ref = spec["ref"], poly, e},
  J = divisorIdeal[F, Bt, red /@ spec["res"], -spec["b"], fr["N"], fr["Ep"], fr["Bp"]];
  Which[
    ref === None, Null,
    ref[[1]] === "Pinf", n = ref[[2]],
    ref[[1]] === "R0",
      poly = kc[F, (fr["a"] - red[ref[[2]]]) t + 1]; e = -ref[[3]];
      J = If[e >= 0,
        idealTimesGenerators[F, Bt, J, {{PadRight[{kc[F, poly^e]}, Bt["m"]]}, 1}],
        idealReduce[F, Bt, <|"M" -> J["M"], "den" -> kMonic[F, kc[F, J["den"] poly^(-e)], t]|>]]];
  {J, n}];

specIdealX::usage = "specIdealX[ctx, spec] is the ideal I(D) of the divisor spec in the working frame over the residue field ctx[\"F\"]; a Pinf term is the place over x = ctx[\"alpha\"], a finite place here (placeIdealAt).";

placeIdealAt::usage = "placeIdealAt[F, B, a] is the prime ideal of the single (totally ramified) place over x = a, a root of G of multiplicity e prime to m: generated by x - a and the uniformiser Y^s/(x - a)^r with s e - r m = 1 (0 < s < m), which is integral since Floor[s e/m] = r makes d_s/(x - a)^r a polynomial.";
placeIdealAt[F_Association, B_Association, a_] := Module[{x = B["x"], m = B["m"], e, s0, r0, pi},
  e = rootMultiplicity[F, B["G"], a, x];
  If[e === 0 || GCD[e, m] =!= 1, throwInternal["placeIdealAtNotTotallyRamified"]];
  s0 = SelectFirst[Range[1, m - 1], Mod[# e - 1, m] === 0 &];
  r0 = (s0 e - 1)/m;
  pi = ReplacePart[Table[0, {m}], s0 + 1 -> kExactQuo[F, B["d"][[s0 + 1]], kc[F, (x - a)^r0], x]];
  idealFromGenerators[F, B, {PadRight[{kc[F, x - a]}, m], pi}]];
specIdealX[ctx_Association, spec_Association] := Module[
  {F = ctx["F"], B = ctx["B"], x, I, ref = spec["ref"], e},
  x = B["x"];
  I = divisorIdeal[F, B, spec["res"], spec["b"], ctx["N"], ctx["Ep"], ctx["Bp"]];
  If[ref =!= None && ref[[1]] === "Pinf" && ref[[2]] =!= 0,
    I = idealTimes[F, B, I, idealPower[F, B, placeIdealAt[F, B, ctx["alpha"]], ref[[2]]]]];
  If[ref =!= None && ref[[1]] === "R0",
    e = ref[[3]];
    I = If[e >= 0,
      idealTimesGenerators[F, B, I, {{PadRight[{kc[F, (x - ref[[2]])^e]}, B["m"]]}, 1}],
      idealReduce[F, B, <|"M" -> I["M"], "den" -> kMonic[F, kc[F, I["den"] (x - ref[[2]])^(-e)], x]|>]]];
  I];

modPFrame::usage = "modPFrame[ctx, gp] is {frame, red}: the reduced frame over GF(p) of the component ctx at the good prime gp (centre gp[\"a\"]) and the reduction map red.";
modPFrame[ctx_Association, gp_Association] := Module[{p = gp["p"], Fp, red, B = ctx["B"], x, c},
  Fp = fieldGF[p]; x = B["x"];
  red = Function[e, reduceResidue[ctx["F"], e, p, gp["rho"]]];
  c = rc[Fp, red[ctx["N"]]/red[ctx["E"]]];
  {jacFrame[Fp, B["m"], gfMod[B["G"], p], c, x, gp["a"]], red}];

orderModPSpec::usage = "orderModPSpec[ctx, spec, gp] is the order over GF(p) of the class of the divisor spec, by reduced class arithmetic (Jacobian.m) in the reduced frame at gp (constant cost per class addition): baby-step giant-step with a bound grown geometrically until the order is found, so the cost is O(Sqrt[order]); exceeding the Weil bound (1 + Sqrt[p])^(2 genus) is an internal error. There is no other limit.";
orderModPSpec[ctx_Association, spec_Association, gp_Association] := Module[{fr, red, Fp, Jn, c, weil, bound, r},
  {fr, red} = modPFrame[ctx, gp];
  Fp = fieldGF[gp["p"]];
  Jn = specFrameIdeal[Fp, fr, spec, red];
  c = jacWithGenerators[Fp, fr["B"], jacClass[Fp, fr["B"], Jn[[1]], Jn[[2]]]];
  weil = Floor[(1 + Sqrt[N[gp["p"]]])^(2 ctx["genus"])] + 1;
  (* Growing search: baby-step giant-step up to bound, the bound multiplied  *)
  (* by 16 until the order is found.  The work is O(Sqrt[order]) -- the     *)
  (* order, not the Weil bound, which is the theorem's ceiling only: past   *)
  (* it the search is an internal error, never a decline.                   *)
  bound = Min[4096, weil];
  While[True,
    r = catchTrager[jacOrderGens[Fp, fr["B"], c, bound]];
    If[!tragerFailQ[r], Return[r, Module]];
    If[r[[1]] =!= "Decline", Throw[r, "TragerFail"]];
    If[bound >= weil, throwInternal["TorsionBeyondWeilBound", <|"p" -> gp["p"], "Bound" -> weil|>]];
    bound = Min[16 bound, weil]]];

torsionOrderSpec::usage = "torsionOrderSpec[ctx, spec] is torsionOrder for a divisor spec ([T] Ch 6 Sec 3). The order op modulo the smaller good prime p is searched (orderModPSpec); modulo q the order must then be n0 p^a q^b with n0 the prime-to-pq part of op, so the least a, b with p^a q^b n0 c = 0 are found by repeated multiplication (O(log) class operations, no second search). Mismatches are NonElementary certificates.";
torsionOrderSpec[ctx_Association, spec_Association] := Module[
  {gps, p, q, op, oq = None, n0p, n0q, fr, red, Fq, Jn, cq, weil, A, Bq, cprime},
  If[ctx["genus"] === 0, Return[1]];
  gps = goodPrimes[ctx, 2];
  {p, q} = #["p"] & /@ gps;
  op = orderModPSpec[ctx, spec, gps[[1]]];
  n0p = op/(pPart[op, p] pPart[op, q]);
  {fr, red} = modPFrame[ctx, gps[[2]]];
  Fq = fieldGF[q];
  Jn = specFrameIdeal[Fq, fr, spec, red];
  cq = jacClass[Fq, fr["B"], Jn[[1]], Jn[[2]]];
  weil = Floor[(1 + Sqrt[N[q]])^(2 ctx["genus"])] + 1;
  A = Floor[Log[p, weil]]; Bq = Floor[Log[q, weil]];
  cprime = If[n0p === 1, cq, jacMultiple[Fq, fr["B"], cq, n0p]];
  (* c' has order p^a' q^b' (if it is torsion of the predicted shape), and    *)
  (* p^a q^b c' = 0 iff a >= a' and b >= b'.  So a' is the least a with        *)
  (* p^a (q^Bq c') = 0 and then b' the least b with q^b (p^a' c') = 0: A + 2 Bq *)
  (* multiplications instead of the (A + 1)(Bq + 1) of a joint search.         *)
  Module[{yq = cprime, xa, a = 0, b = 0, ok = True},
    Do[yq = jacMultiple[Fq, fr["B"], yq, q], {Bq}];
    xa = yq;
    While[!jacZeroQ[xa] && a < A, a++; xa = jacMultiple[Fq, fr["B"], xa, p]];
    If[!jacZeroQ[xa], ok = False];
    If[ok,
      xa = cprime;
      Do[xa = jacMultiple[Fq, fr["B"], xa, p], {a}];
      While[!jacZeroQ[xa] && b < Bq, b++; xa = jacMultiple[Fq, fr["B"], xa, q]];
      If[jacZeroQ[xa], oq = n0p p^a q^b]]];
  If[oq === None, throwNonElementary["residue divisor of infinite order: its orders modulo " <>
      ToString[p] <> " and " <> ToString[q] <> " disagree", <|"OrderModP" -> op, "Primes" -> {p, q}|>]];
  (* the prime-to-pq part modulo q must be all of n0p, not a proper divisor *)
  Do[
    If[jacZeroQ[jacMultiple[Fq, fr["B"], cq, oq/r]],
      throwNonElementary["residue divisor of infinite order: its orders modulo " <>
        ToString[p] <> " and " <> ToString[q] <> " disagree", <|"OrderModP" -> op, "Primes" -> {p, q}|>]],
    {r, If[n0p === 1, {}, First /@ FactorInteger[n0p]]}];
  n0q = oq/(pPart[oq, p] pPart[oq, q]);
  If[n0p =!= n0q, throwInternal["TorsionPrimeToPQMismatch"]];
  n0p pPart[oq, p] pPart[op, q]];

torsionScreen::usage = "torsionScreen[ctx, spec, Nn] is False when Nn delta (delta the divisor spec) is shown NOT principal by reduction: Nn times its class is nonzero over GF(r) at one of tragerOption[\"TorsionScreenPrimes\"] good primes r beyond the two that determined Nn. Reduction at a good prime is a homomorphism of divisor class groups, so Nn delta principal forces Nn delta-bar principal; a False is therefore a proof, and with Nn the only possible order (torsionOrderSpec) a proof that delta has infinite order. True means no prime refuted it (the characteristic-0 test still decides). Genus 0 is always True.";
torsionScreen[ctx_Association, spec_Association, Nn_Integer] := Module[
  {found, want, fr, red, Fr, Jn, c},
  If[ctx["genus"] === 0, Return[True]];
  want = 2 + tragerOption["TorsionScreenPrimes"];
  (* the good primes in order, as goodPrimes picks them (the first two fixed Nn) *)
  found = goodPrimes[ctx, want];
  Catch[
    Do[
      {fr, red} = modPFrame[ctx, g];
      Fr = fieldGF[g["p"]];
      Jn = specFrameIdeal[Fr, fr, spec, red];
      c = jacClass[Fr, fr["B"], Jn[[1]], Jn[[2]]];
      If[!jacZeroQ[If[Nn === 1, c, jacMultiple[Fr, fr["B"], c, Nn]]], Throw[False, "torsionScreen"]],
      {g, Drop[found, Min[2, Length[found]]]}];
    True, "torsionScreen"]];
