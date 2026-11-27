(* ========================================================================== *)
(* goursat/RatInt.m -- rational integration for the GeneralisedGoursat package *)
(*                                                                            *)
(* Every symmetry reduction ends in a rational integrand rho(m) in a uniformiser *)
(* m that is a rational function of t and y.  goursat.py: rat_int, _rat_int_alg. *)
(* Every candidate antiderivative is verified (PossibleZeroQ of P' - rho) before *)
(* it is used; there is no numerical fallback.                                *)
(* ========================================================================== *)

ratIntBudget[] := 4;      (* seconds for the Integrate candidate (sympy: 4) *)

(* ---- deferred rational integration ---------------------------------------- *)
(* The verdict of a reduction depends only on the residuals, never on the        *)
(* antiderivatives of the pieces.  The Integrate stage (GoursatIntegrate)         *)
(* therefore runs the reduction with $ggDefer = True: ratInt records rho and      *)
(* returns the placeholder ggDeferred[id][m] (a function of the uniformiser, so   *)
(* the substitutions m -> m(t, y) and the second-level back-substitutions apply  *)
(* to it unchanged), and the integrals are computed only once the verdict is      *)
(* elementary or certified non-elementary.  An undecided reduction -- a third-    *)
(* kind residual at the last level -- then costs no rational integration at all: *)
(* (2 + x)/((1 + x + x^2) (2 + x^3)^(1/3)) spent 25 s integrating the pieces of  *)
(* three second-level reductions whose answer was then discarded.                *)
$ggDefer = False;
$ggDeferred = {};

(* A table entry is {thunk, vars}: thunk[] computes the piece in the symbols   *)
(* vars, and ggDeferred[id][args] stands for it at vars -> args.              *)
(* info: {"rat", rho, m} for a rational integral (the deferred verifier reads  *)
(* its derivative from it), {} otherwise.                                       *)
deferPiece[thunk_, vars_List, info_List : {}] := (
  AppendTo[$ggDeferred, {thunk, vars, info}];
  ggDeferred[Length[$ggDeferred]] @@ vars);

ratInt[rho0_, m_] /; TrueQ[$ggDefer] := Module[{rho = gK[rho0]},
  (* With: the thunk must hold the values, not the Module's local symbols *)
  If[rho === 0, 0, With[{r = rho}, deferPiece[ratInt[r, m] &, {m}, {"rat", r, m}]]]];

(* Every placeholder replaced by its piece, at the placeholder's arguments. *)
(* The pieces are computed once each, in table order; a rational integrand    *)
(* that is the Galois conjugate of one already integrated (the two Hessian-    *)
(* root components of the cube-root reduction, conjugate over Q(Sqrt[D])) takes *)
(* the conjugate antiderivative, which is checked (ratOK) like any other        *)
(* candidate, instead of a second run of the rational integrator.              *)
resolveDeferred[e_] := Module[{ids, val = <||>, done = {}},
  If[FreeQ[e, ggDeferred], Return[e, Module]];
  ids = Union[Cases[e, ggDeferred[id_Integer][___] :> id, {0, Infinity}, Heads -> True]];
  Block[{$ggDefer = False},
    Do[With[{tv = $ggDeferred[[id]]},
         val[id] = If[MatchQ[tv[[3]], {"rat", _, _}],
           With[{P = conjugatePiece[tv[[3]][[2]], tv[[3]][[3]], done]},
             If[P =!= None, P,
               With[{Q = tv[[1]][]}, AppendTo[done, {tv[[3]][[2]], tv[[3]][[3]], Q}]; Q]]],
           tv[[1]][]]],
       {id, ids}]];
  e /. ggDeferred[id_Integer][args___] :>
    (val[id] /. Thread[$ggDeferred[[id]][[2]] -> {args}])];

(* The automorphisms tried: complex conjugation, and the sign change of one     *)
(* square root of a positive integer.  Complex conjugation fixes a real radical *)
(* and sends I to -I and a radical of a negative number, (-c)^e = c^e E^(I Pi e), *)
(* to c^e E^(-I Pi e) = (-c)^e (-1)^(-2 e).                                      *)
conjugationRules[rho_] := Join[
  If[FreeQ[rho, Complex], {}, {{Complex[a_, b_] :> Complex[a, -b],
                               Power[b_?Negative, e_Rational] :> Power[b, e] Power[-1, -2 e]}}],
  {Power[#, 1/2] -> -Power[#, 1/2], Power[#, -1/2] -> -Power[#, -1/2]} & /@
    Union[Cases[rho, Power[b_Integer?Positive, 1/2 | -1/2] :> b, {0, Infinity}]]];

(* The automorphism applied structurally, every new node evaluated: inside the  *)
(* held body of a RootSum's Function a plain ReplaceAll leaves Complex[0, -1*-1] *)
(* unevaluated (an answer the printer writes in an unparsable form).             *)
applyAuto[e_, rule_] := Which[
  MatchQ[e, Alternatives @@ (First /@ rule)], e /. rule,
  AtomQ[e], e,
  Head[e] === RootSum && Length[e] == 2, RootSum @@ {applyAuto[e[[1]], rule], applyAuto[e[[2]], rule]},
  Head[e] === Function, Function @@ Append[Most[List @@ e], applyAuto[Last[e], rule]],
  True, applyAuto[#, rule] & /@ e];

(* An antiderivative of rho in m from an integrated conjugate, verified; None. *)
conjugatePiece[rho_, m_, done_List] := Module[{c},
  Do[With[{rho1 = d[[1]] /. d[[2]] -> m},
       Do[If[gK[applyAuto[rho1, rule]] === rho,
            c = applyAuto[d[[3]] /. d[[2]] -> m, rule];
            (* a short budget: a candidate that is not one is rejected, not ground *)
            If[stage["ratInt: conjugate check", pzqTrue[D[c, m] - rho, m, 2]], Return[c, Module]]],
          {rule, conjugationRules[rho1]}]],
     {d, done}];
  None];

(* Antiderivative of a rational function rho of m (coefficients: Q(algebraics)). *)
ratInt[rho0_, m_] := Module[{rho, P},
  rho = gK[rho0];
  If[rho === 0, Return[0, Module]];
  (* 1. the system rational integrator (sympy: ratint), over Q and over       *)
  (*    Q(algebraics) alike: its log part works natively in the number field.  *)
  P = stage["ratInt: Integrate", Quiet[TimeConstrained[
        Integrate[rho, m, Method -> "BronsteinRational"], ratIntBudget[], $Failed]]];
  If[P =!= $Failed && FreeQ[P, Integrate] && stage["ratInt: check", ratOK[P, rho, m]],
    Return[P, Module]];
  (* 2. Hermite's rational part plus logarithms summed over the roots of the *)
  (*    rational squarefree norm of the denominator                           *)
  P = stage["ratInt: RootSum", Catch[ratIntAlg[rho, m], "GoursatFail"]];
  If[!MatchQ[P, _ggFail] && stage["ratInt: check", ratOK[P, rho, m]], Return[P, Module]];
  ggDecline["rational integration"]];

(* P' == rho, exactly: the RootSum-aware derivative of Verify.m over the     *)
(* trivial extension (n = 1, R = 1), PossibleZeroQ for anything else.        *)
(* A RootSum-free answer is put to PossibleZeroQ first: the system integrator's  *)
(* closed forms carry unrelated-looking radicals (2^(1/3), 3^(1/6), 6^(1/3), ...) *)
(* that the vector arithmetic, which does not canonicalise, lets swell.          *)
(* PossibleZeroQ goes first for a RootSum answer too, under a short budget:    *)
(* the vector check of a closed form carrying (1/k)^(1/3) over Q(k)(I Sqrt[3])  *)
(* runs past a minute (PossibleZeroQ: 30 ms), and a 471-leaf answer over Q(alg) *)
(* (stress_7 sqrt-nonelementary #6) past 30 s (PossibleZeroQ: 17 ms).           *)
ratOK[P_, rho_, m_] := Or[
  pzqTrue[D[P, m] - rho, m,
    Which[FreeQ[P, _RootSum], 20, parametricQ[{P, rho}, m, ggY], 5, True, 2]],
  TrueQ[Quiet[Catch[derivZeroQ[P, -rho, m, ggY, 1, 1], "GoursatFail"]]]];

(* Horowitz-Ostrogradsky: r/den = (P/D1)' + Q/D2 with D1 = gcd(den, den'),      *)
(* D2 = den/D1 squarefree, deg P < deg D1, deg Q < deg D2.  Returns {P/D1, Q/D2}. *)
hermiteRat[r_, den_, m_] := Module[{D1, D2, dP, dQ, cp, cq, P, Q, sol},
  D1 = PolynomialGCD[den, D[den, m], Extension -> Automatic];
  If[Exponent[D1, m] == 0, Return[{0, gK[r/den]}, Module]];
  D2 = gK[den/D1];
  dP = Exponent[D1, m]; dQ = Exponent[D2, m];
  P = Sum[cp[i] m^i, {i, 0, dP - 1}];
  Q = Sum[cq[i] m^i, {i, 0, dQ - 1}];
  (* r = P' D2 - P (D1' D2 / D1) + Q D1 *)
  sol = solveIdentity[
    D[P, m] D2 - P gK[D[D1, m] D2/D1] + Q D1 - r,
    Join[Table[cp[i], {i, 0, dP - 1}], Table[cq[i], {i, 0, dQ - 1}]], m];
  {gK[(P /. sol)/D1], gK[(Q /. sol)/D2]}];

(* Norm over Q of a polynomial B in m whose coefficients are algebraic numbers: *)
(* the generators (I, radicals b^(p/q), Root objects) are eliminated one at a  *)
(* time by resultants.  An unneeded conjugate factor only adds roots at which  *)
(* the residue below vanishes, so a multiple of the true norm is harmless.     *)
algNorm[B0_, m_] := Module[{B = Expand[B0], z, gens, g, base, L, mu, guard = 0},
  While[guard++ < 50,
    (* 1. a radical of any numeric base, the outermost first: its base may hold  *)
    (*    I and further radicals (the conic step's                              *)
    (*    Sqrt[(a + b I Sqrt[46])/(c + d I Sqrt[46])]), which later passes        *)
    (*    eliminate -- eliminating those first would rewrite them inside the     *)
    (*    radicand, which then is no longer a number                            *)
    gens = Union[Cases[B, Power[b_?NumericQ, e_Rational] /; Denominator[e] > 1 :> b, {0, Infinity}]];
    If[gens =!= {},
      base = Last[SortBy[gens, LeafCount]];
      L = LCM @@ Cases[B, Power[base, e_Rational] :> Denominator[e], {0, Infinity}];
      B = Expand[B /. Power[base, e_Rational] :> z^(e L)];
      B = Expand[B z^(-Min[0, Exponent[B, z, Min]])];
      B = Expand[Resultant[Together[z^L - base] Denominator[Together[base]], B, z]];
      Continue[]];
    (* 2. I *)
    If[!FreeQ[B, Complex],
      B = Expand[B /. Complex[a_, b_] :> a + b ggI];
      B = Expand[Resultant[z^2 + 1, Expand[B /. ggI -> z], z]];
      Continue[]];
    (* 3. Root objects *)
    gens = Union[Cases[B, _Root, {0, Infinity}]];
    If[gens =!= {},
      g = First[gens]; mu = MinimalPolynomial[g, z];
      B = Expand[Resultant[mu, Expand[B /. g -> z], z]];
      Continue[]];
    Break[]];
  If[!FreeQ[B, ggI | z | Complex | Power[_?NumericQ, _Rational]], ggDecline["algNorm"]];
  B = Together[B];
  Expand[Numerator[B]]];

ratIntAlg[rho_, m_] := Module[{num, den, q, r, out, C, rest, A, B, nrm, nsf, cof, inv, c},
  {num, den} = numDen[rho];
  {q, r} = PolynomialQuotientRemainder[num, den, m];
  out = Integrate[Expand[q], m];
  {C, rest} = hermiteRat[r, den, m];
  out += C;
  rest = gK[rest];
  If[rest === 0, Return[out, Module]];
  {A, B} = numDen[rest];
  nrm = algNorm[B, m];
  (* the squarefree part, restricted to the irreducible factors over Q that    *)
  (* share a root with B: the other conjugates carry zero residues (sympy keeps *)
  (* them; dropping them halves the RootSum degree and changes nothing else)    *)
  (* ... but only for fixed numbers: coefficients in a cubic field (Root       *)
  (* objects) make the answer a formula that CubicField.m conjugates, and a      *)
  (* conjugate's roots can lie in a factor the first embedding does not use --   *)
  (* so there the whole squarefree norm is kept, as sympy does                   *)
  nsf = If[!FreeQ[{A, B}, _Root],
    Times @@ Cases[FactorList[nrm], {f_, _} /; Exponent[f, m] > 0 :> f],
    Times @@ Cases[FactorList[nrm],
      {f_, _} /; Exponent[f, m] > 0 &&
                 Exponent[PolynomialGCD[f, B, Extension -> Automatic], m] > 0 :> f]];
  If[Exponent[nsf, m] < Exponent[B, m], ggInternal["ratIntAlg: norm lost a root"]];
  nsf = Expand[nsf/Coefficient[Expand[nsf], m, Exponent[nsf, m]]];
  cof = gK[nsf/B];
  If[!PolynomialQ[cof, m], ggInternal["ratIntAlg: B does not divide its norm"]];
  inv = PolynomialExtendedGCD[D[nsf, m], nsf, m][[2, 1]];      (* over Q *)
  c = With[{cs = polyrem[Expand[A cof inv], nsf, m]},
    cs.(m^Range[Length[cs] - 1, 0, -1])];
  With[{p = nsf /. m -> Slot[1], s = (c /. m -> Slot[1]) Log[m - Slot[1]]},
    out + RootSum[Function[p], Function[s]]]];
