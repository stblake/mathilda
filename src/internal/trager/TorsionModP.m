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

goodPrime::usage = "goodPrime[ctx, p] tests the prime p for good reduction of the component data ctx (keys \"B\", \"genus\", \"F\", \"residues\", \"N\", \"Ep\", \"Bp\"). It returns <|\"p\" -> p, \"rho\" -> rho|> (rho the image of the residue-field generator, None over Q) or None.";
goodPrime[ctx_Association, p_Integer] := Module[
  {B = ctx["B"], x, m, G, Gs, Bp = ctx["Bp"], Ep = ctx["Ep"], F = ctx["F"], rho = None, lin, alphas},
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
  <|"p" -> p, "rho" -> rho|>];

reduceResidue::usage = "reduceResidue[F, a, p, rho] maps the element a of the residue field F to GF(p), sending the generator t to rho.";
reduceResidue[F_Association, a_, p_Integer, rho_] :=
  gfMod[If[F["Type"] === "NF", kc[F, a] /. F["t"] -> rho, a], p];

orderModP::usage = "orderModP[ctx, b, gp] is the order of the class of the divisor Sum_i b_i Delta(alpha_i) of the component ctx after reduction at the good prime gp = <|\"p\" -> p, \"rho\" -> rho|>: the least M >= 1 with M delta principal over GF(p). Throws a Decline past tragerOption[\"MaxTorsionSearch\"], and an internal error past the Weil bound (1 + Sqrt[p])^(2 genus).";
orderModP[ctx_Association, b_List, gp_Association] := orderModPWith[ctx, residueCombinationBuilder[ctx, b], gp];

residueCombinationBuilder::usage = "residueCombinationBuilder[ctx, b] is the ideal builder (see orderModPWith) of the divisor Sum_i b_i Delta(alpha_i): a function of (F, B, red) returning the product of the residue ideals over F, every residue-field element and polynomial passed through red first.";
residueCombinationBuilder[ctx_Association, b_List] := Function[{Fx, Bx, red},
  divisorIdeal[Fx, Bx, red /@ ctx["residues"], b, red[ctx["N"]], red[ctx["Ep"]], red[ctx["Bp"]]]];

orderModPWith::usage = "orderModPWith[ctx, mk, gp] is the order over GF(p) of the divisor class whose ideal mk[F, B, red] builds: mk is called with the field GF(p), the reduced basis, and the reduction map red (residue-field generator -> rho, coefficients mod p). Same bounds as orderModP.";
orderModPWith[ctx_Association, mk_, gp_Association] := Module[
  {p = gp["p"], Fp, Bm, red, I, J, weil, cap, M = 1},
  Fp = fieldGF[p]; Bm = basisMod[ctx["B"], p];
  red = Function[e, reduceResidue[ctx["F"], e, p, gp["rho"]]];
  I = mk[Fp, Bm, red];
  weil = Floor[(1 + Sqrt[N[p]])^(2 ctx["genus"])] + 1;
  cap = tragerOption["MaxTorsionSearch"];
  J = I;
  While[!principalQ[Fp, Bm, J],
    M++;
    If[M > weil, throwInternal["TorsionBeyondWeilBound", <|"p" -> p, "Bound" -> weil|>]];
    If[M > cap, throwDecline["torsion order search bound reached", <|"p" -> p, "Bound" -> cap|>]];
    J = idealTimes[Fp, Bm, J, I]];
  M];

goodPrimes::usage = "goodPrimes[ctx, k] is a list of k distinct good primes (as goodPrime records) for the component ctx, smallest first. Throws a Decline when fewer than k are found among the first tragerOption[\"MaxPrimes\"] primes.";
goodPrimes[ctx_Association, k_Integer] := Module[{found = {}, p = 2, tried = 0, gp},
  While[Length[found] < k && tried < tragerOption["MaxPrimes"],
    p = NextPrime[p]; tried++;
    gp = goodPrime[ctx, p];
    If[gp =!= None, AppendTo[found, gp]]];
  If[Length[found] < k, throwDecline["no primes of good reduction found"]];
  found];

torsionOrder::usage = "torsionOrder[ctx, b] is the only possible order N of the divisor class delta = Sum_i b_i Delta(alpha_i) of the component ctx, from its orders modulo two good primes ([T] Ch 6 Sec 3). It throws a NonElementary certificate when the prime-to-pq parts of the two orders differ (delta has infinite order). The caller must still test N delta for principality over the residue field: if that fails, delta has infinite order too.";
torsionOrder[ctx_Association, b_List] := torsionOrderWith[ctx, residueCombinationBuilder[ctx, b]];

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
