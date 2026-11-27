(* ========================================================================== *)
(* goursat/Verify.m -- exact verification of a GeneralisedGoursat answer.      *)
(*                                                                            *)
(* goursat.py checks d/dt(antiderivative) + residual/y^k = F/y^k numerically   *)
(* at random points.  Here the check is exact: the derivative is taken          *)
(* formally (y' = R'/(n y^(n-1))) and reduced to the normal form                *)
(*   c_0(t) + c_1(t) y + ... + c_(n-1)(t) y^(n-1)  in Q(t)[y]/(y^n - R),         *)
(* term by term, as coefficient VECTORS.  A RootSum term                        *)
(*   RootSum[p, c(#) Log[U - #] &]  (the shape every rational integration here   *)
(*   produces, ratIntAlg)                                                       *)
(* is differentiated without ever forming the RootSum's symbolic collapse:      *)
(*   d/dt = U' Sum_a c(a)/(U - a),                                              *)
(*   Sum_a a^k/(U - a) = U^k p'(U)/p(U) - Sum_(j<k) U^(k-1-j) s_j               *)
(* (s_j the power sums of the roots of p), every quantity evaluated by Horner's *)
(* rule in the vector arithmetic of Q(t)[y]/(y^n - R) -- where U = 2y/(2t + 1   *)
(* + I Sqrt[3]) costs nothing, though the same p(U) put through Together has a  *)
(* degree-8 numerator over Q(I Sqrt[3], t, y) and takes a minute.  A nested     *)
(* RootSum (CubicField.m) is summed over its cubic by a companion-matrix trace. *)
(* Anything else falls back to PossibleZeroQ on the radical form.               *)
(* ========================================================================== *)

(* PossibleZeroQ under a time budget (the verifier of record; a defect in it is *)
(* fixed there, not patched around here).                                       *)
pzqTrue[e_, t_, budget_] := TrueQ[Quiet[TimeConstrained[PossibleZeroQ[e], budget, False]]];

(* Formal derivative of an expression in (t, y), y^n = R. *)
totalD[e_, t_, y_, R_, n_] := D[e, t] + D[e, y] D[R, t]/(n y^(n - 1));

(* ---- vector arithmetic in Q(t)[y]/(y^n - R) -------------------------------- *)
(* The vector arithmetic needs no canonical forms, only an exact zero test at  *)
(* the end, so it uses plain Together (vK): no algebraic-extension gcds, which  *)
(* over a field like Q(I, Sqrt[3], 2^(1/3)) cost seconds per operation.         *)
(* Below degree 4 the extension gcds are cheap and stop expression swell, so the *)
(* canonical gK is used there.                                                  *)
vK[e_] := If[algDegree[e] <= 4, gK[e], Together[e]];
algDegree[e_] := Times @@ Join[
  If[FreeQ[e, Complex], {}, {2}],
  Union[Cases[e, Power[b_?NumericQ, q_Rational] :> {b, Denominator[q]}, {0, Infinity}]][[All, 2]],
  Cases[Union[Cases[e, _Root, {0, Infinity}]], r_ :> Exponent[r[[1]][ggZ], ggZ]]];
(* The vector of e: read off directly when the denominator is free of y (every  *)
(* argument tidyArgs produced), the exact normal form otherwise.                *)
vecOf[e_, R_, t_, y_, n_] := Module[{r = Together[e], num, den, cl, v},
  num = Numerator[r]; den = Denominator[r];
  If[!FreeQ[den, y] || !PolynomialQ[num, y], Return[normalizeCoeffs[e, R, t, y, n], Module]];
  cl = CoefficientList[Expand[num], y];
  v = ConstantArray[0, n];
  Do[v[[Mod[j - 1, n] + 1]] += cl[[j]] R^Quotient[j - 1, n], {j, Length[cl]}];
  vK /@ (v/den)];
constV[c_, n_] := Prepend[ConstantArray[0, n - 1], c];
addV[a_, b_] := vK /@ (a + b);
scaleV[c_, a_] := vK /@ (c a);
mulV[a_, b_, R_, n_] := Module[{c = ConstantArray[0, n]},
  Do[If[a[[i + 1]] =!= 0 && b[[j + 1]] =!= 0,
       With[{k = i + j},
         If[k < n, c[[k + 1]] += a[[i + 1]] b[[j + 1]],
                   c[[k - n + 1]] += R a[[i + 1]] b[[j + 1]]]]],
     {i, 0, n - 1}, {j, 0, n - 1}];
  vK /@ c];
(* Inverse: the first column of adj(M)/det(M), M the multiplication matrix. *)
invV[a_, R_, n_] := Module[{M = ConstantArray[0, {n, n}], dt, ad},
  If[n == 1 || AllTrue[Rest[a], # === 0 &], Return[constV[vK[1/a[[1]]], n], Module]];
  Do[With[{k = i + j},
       M[[Mod[k, n] + 1, j + 1]] += a[[i + 1]] R^Quotient[k, n]],
     {i, 0, n - 1}, {j, 0, n - 1}];
  If[n > 3, Return[gK /@ linSolve[Map[gK, M, {2}], constV[1, n]], Module]];
  dt = vK[Det[M]];
  ad = adjN[M];
  vK /@ (ad[[All, 1]]/dt)];
(* Horner: Sum_i cs[[i]] U^(len - i), cs descending. *)
(* A monomial U = alpha y^e (every uniformiser y/(t - gamma), H y/R, v/(X - x1)) *)
(* needs no Horner steps: Sum_j a_j alpha^j y^(e j) splits by e j mod n into n   *)
(* sums, each put through vK once -- where Horner canonicalises 18 growing        *)
(* intermediate results (a degree-18 norm over Q(I Sqrt[3]): minutes).           *)
hornerV[cs_List, U_, R_, n_] := Module[{acc = constV[0, n], nz, e, al, d = Length[cs] - 1, out},
  nz = Flatten[Position[U, x_ /; x =!= 0, {1}, Heads -> False]];
  (* only for a small alpha: the powers of a large one (a back-substituted       *)
  (* second-level uniformiser) are cheaper canonicalised step by step           *)
  If[Length[nz] == 1 && LeafCount[U[[nz[[1]]]]] <= 40,
    e = nz[[1]] - 1; al = U[[e + 1]];
    out = ConstantArray[0, n];
    Do[With[{a = cs[[d - j + 1]]},
         If[a =!= 0, out[[Mod[e j, n] + 1]] += a al^j R^Quotient[e j, n]]],
       {j, 0, d}];
    Return[vK /@ out, Module]];
  Do[acc = addV[mulV[acc, U, R, n], constV[c, n]], {c, cs}];
  acc];

(* Derivative of an element given as a vector: d/dt Sum c_j y^j               *)
(*   = Sum (c_j' + j c_j R'/(n R)) y^j,   since y' = R' y/(n R).               *)
dVec[v_, t_, R_, n_] := Table[vK[D[v[[j + 1]], t] + j v[[j + 1]] D[R, t]/(n R)], {j, 0, n - 1}];

transHeads[] := Log | ArcTan | ArcTanh | ArcCot | ArcCoth;

(* The derivative of h[w] for a transcendental head h, all in vector arithmetic:  *)
(* w'/w, w'/(1 + w^2), w'/(1 - w^2), -w'/(1 + w^2), w'/(1 - w^2).                *)
transDerivVec[h_, w_, t_, y_, R_, n_] := Module[{wv, dw, one = constV[1, n], w2},
  wv = vecOf[w, R, t, y, n];
  dw = dVec[wv, t, R, n];
  w2 = mulV[wv, wv, R, n];
  Switch[h,
    Log, mulV[dw, invV[wv, R, n], R, n],
    ArcTan, mulV[dw, invV[addV[one, w2], R, n], R, n],
    ArcCot, scaleV[-1, mulV[dw, invV[addV[one, w2], R, n], R, n]],
    ArcTanh | ArcCoth, mulV[dw, invV[addV[one, scaleV[-1, w2]], R, n], R, n]]];

(* ---- derivative of a term, as a vector --------------------------------------- *)
derivVec[e_, t_, y_, R_, n_] := Module[{rs, c, d, tr},
  Which[
    FreeQ[e, transHeads[] | _RootSum | _Integrate],      (* algebraic *)
      dVec[vecOf[e, R, t, y, n], t, R, n],
    MatchQ[Head[e], transHeads[]] && Length[e] == 1 && FreeQ[e[[1]], transHeads[] | _RootSum],
      transDerivVec[Head[e], e[[1]], t, y, R, n],
    Head[e] === Times && Length[Select[List @@ e, !FreeQ[#, t | y] &]] == 1 &&
      MatchQ[tr = First[Select[List @@ e, !FreeQ[#, t | y] &]], _?(MatchQ[Head[#], transHeads[]] &)],
      scaleV[Times @@ Select[List @@ e, FreeQ[#, t | y] &], derivVec[tr, t, y, R, n]],
    FreeQ[e, _RootSum], vecOf[Together[totalD[e, t, y, R, n]], R, t, y, n],
    Head[e] === Plus, Module[{parts = derivVec[#, t, y, R, n] & /@ List @@ e},
      If[MemberQ[parts, $Failed], $Failed, Fold[addV, constV[0, n], parts]]],
    Head[e] === Times,
      rs = Select[List @@ e, !FreeQ[#, _RootSum] &];
      c = Times @@ Select[List @@ e, FreeQ[#, _RootSum] &];
      If[Length[rs] != 1 || !FreeQ[c, t | y], Return[$Failed, Module]];
      d = derivVec[First[rs], t, y, R, n];
      If[d === $Failed, $Failed, scaleV[c, d]],
    Head[e] === RootSum && Length[e] == 2, rootSumDerivVec[e, t, y, R, n],
    True, $Failed]];

rootSumDerivVec[e_, t_, y_, R_, n_] := Module[
  {z, p, d, body, c, w, w0, w1, U, Uv, s, cs, Sm, Q, sum, bv},
  p = Expand[e[[1]][z]];
  p = Expand[p/Coefficient[p, z, Exponent[p, z]]];
  d = Exponent[p, z];
  body = e[[2]][z];
  (* the log shape c(z) Log[w], w = w0 + w1 z, c and w1 free of t and y *)
  {c, w} = Replace[body, {cc_. Log[ww_] /; FreeQ[cc, t | y] :> {cc, ww}, _ -> {None, None}}];
  If[c =!= None && ratInQ[c, z] && PolynomialQ[w, z] && Exponent[w, z] == 1 &&
     FreeQ[w1 = Coefficient[w, z, 1], t | y] &&
     !FreeQ[w0 = Coefficient[w, z, 0], t | y],
    U = gK[-w0/w1];                                 (* w = -w1 (U - z) *)
    Uv = vecOf[U, R, t, y, n];
    s = powerSums[p, z, d];
    cs = polyrem[c, p, z];                          (* descending, degree < d *)
    (* S(m) = Sum_a c(a)/(m - a) = c(m) p'(m)/p(m) - Q(m),                       *)
    (*   Q(m) = Sum_k c_k Sum_(j<k) s_j m^(k-1-j),                              *)
    (* formed and canonicalised in ONE variable m first: the residues of a norm  *)
    (* RootSum are hundreds of digits long and cancel down to the (small)        *)
    (* rational integrand, which is then evaluated at U.                         *)
    Q = Sum[Coefficient[cs.(z^Range[Length[cs] - 1, 0, -1]), z, k] *
            Sum[s[[j + 1]] z^(k - 1 - j), {j, 0, k - 1}], {k, 1, d - 1}];
    Sm = gK[(cs.(z^Range[Length[cs] - 1, 0, -1])) D[p, z]/p - Q];
    sum = mulV[hornerV[allCoeffs[Numerator[Sm], z], Uv, R, n],
               invV[hornerV[allCoeffs[Denominator[Sm], z], Uv, R, n], R, n], R, n];
    Return[mulV[sum, vecOf[Together[totalD[U, t, y, R, n]], R, t, y, n], R, n], Module]];
  (* any other body: differentiate it as a vector with z symbolic, then sum each *)
  (* entry over the roots of p                                                    *)
  bv = derivVec[body, t, y, R, n];
  If[bv === $Failed, Return[$Failed, Module]];
  sum = rootSumOf[p, z, #] & /@ bv;
  If[MemberQ[sum, $Failed], $Failed, gK /@ sum]];

(* Sum_{p(z)=0} g(z) for g rational in z: power sums when the denominator is free *)
(* of z, a companion-matrix trace (deg p <= 3) otherwise.                        *)
rootSumOf[p_, z_, g_] := Module[{d = Exponent[p, z], A, B, cs, s, Th, pows, mat},
  If[FreeQ[g, z], Return[d g, Module]];
  (* soundness: the trace and the power sums are valid only for g RATIONAL in z; *)
  (* a radical of z (Sqrt[A(z)] in a conjugated cubic-field piece) is not a       *)
  (* coefficient                                                                  *)
  If[!ratInQ[g, z], Return[$Failed, Module]];
  A = Numerator[Together[g]]; B = Denominator[Together[g]];
  Which[
    FreeQ[B, z],
      s = powerSums[p, z, d];
      cs = Reverse[polyrem[A, p, z]];
      gK[cs.s[[1 ;; Length[cs]]]/B],
    d <= 3,
      Th = companion[p, z];
      pows = NestList[Th.# &, IdentityMatrix[d], d - 1];
      mat[f_] := With[{c = Reverse[polyrem[f, p, z]]}, Sum[c[[i]] pows[[i]], {i, Length[c]}]];
      With[{Ma = mat[A], Mb = mat[B]}, gK[Tr[Ma.adjN[Mb]]/Det[Mb]]],
    True, $Failed]];

(* Is e a rational function of z: z occurs only in sums, products and integer  *)
(* powers (never under a radical, a logarithm or any other head)?               *)
ratInQ[e_, z_] := FreeQ[e /. {Power[b_, _Integer] :> b},
  (h_ /; !MemberQ[{Plus, Times, Power, List, Complex, Rational}, h])[___, a_ /; !FreeQ[a, z], ___] |
  Power[a_ /; !FreeQ[a, z], _]];

(* Power sums s_0 .. s_(K-1) of the roots of the monic p (Newton's identities):  *)
(* with p = z^d + a1 z^(d-1) + ... + ad,                                        *)
(*   s_k = -(a1 s_(k-1) + ... + a_(k-1) s_1 + k a_k)   for k <= d,               *)
(*   s_k = -(a1 s_(k-1) + ... + a_d s_(k-d))             for k > d.                *)
powerSums[p_, z_, K_] := Module[{d = Exponent[p, z], c, a, s},
  c = Reverse[allCoeffs[p, z]];                 (* c[[i+1]] = coefficient of z^i *)
  a[j_] := c[[d - j + 1]];
  s = ConstantArray[0, K];
  If[K >= 1, s[[1]] = d];
  Do[s[[k + 1]] = If[k <= d,
       -(Sum[a[j] s[[k - j + 1]], {j, 1, k - 1}] + k a[k]),
       -Sum[a[j] s[[k - j + 1]], {j, 1, d}]],
     {k, 1, K - 1}];
  s];

(* Companion matrix of the monic p (multiplication by z on 1, z, ..., z^(d-1)). *)
companion[p_, z_] := Module[{d = Exponent[p, z], c = Reverse[allCoeffs[p, z]]},
  Table[Which[j == d, -c[[i]], i == j + 1, 1, True, 0], {i, d}, {j, d}]];

(* Adjugate of a square matrix of size <= 3. *)
adjN[M_] := Switch[Length[M],
  1, {{1}},
  2, {{M[[2, 2]], -M[[1, 2]]}, {-M[[2, 1]], M[[1, 1]]}},
  3, adj3[M]];

(* ---- the checks ------------------------------------------------------------- *)
(* Is  D[anti] + rest  zero on y^n = R?  Terms whose derivative has no vector    *)
(* form are collected and decided, with the rest, by PossibleZeroQ.             *)
(* Two verifiers, cheapest first: PossibleZeroQ on the radical form (seconds   *)
(* at most for a closed form, but it false-negatives on nested RootSums), then  *)
(* the exact vector check (which, not canonicalising, can swell over many        *)
(* unrelated radicals -- hence its budget).                                      *)
pzqBudget[] := 5;
vecBudget[] := 60;
derivZeroQ[anti_, rest_, t_, y_, R_, n_] := Module[{v, terms, rs, cl, rv},
  terms = If[Head[anti] === Plus, List @@ anti, {anti}];
  rs = Select[terms, !FreeQ[#, _RootSum] &];
  cl = Complement[terms, rs];
  (* 1. no RootSum: PossibleZeroQ on the radical form *)
  If[rs === {} && pzqTrue[(totalD[anti, t, y, R, n] + rest) /. y -> R^(1/n), t, pzqBudget[]],
    Return[True, Module]];
  (* 2. hybrid: the RootSum terms differentiated exactly (vector form, so no    *)
  (*    RootSum is left), the closed-form terms -- whose radicals the vector     *)
  (*    arithmetic would let swell -- left to PossibleZeroQ with the rest        *)
  If[rs =!= {} && cl =!= {},
    rv = Quiet[TimeConstrained[Catch[
      Fold[addV, constV[0, n], derivVec[#, t, y, R, n] & /@ rs], "GoursatFail"], vecBudget[], $Failed]];
    If[ListQ[rv] && FreeQ[rv, $Failed] &&
       pzqTrue[(totalD[Plus @@ cl, t, y, R, n] + rv.(y^Range[0, n - 1]) + rest) /. y -> R^(1/n),
               t, pzqBudget[]],
      Return[True, Module]]];
  (* 3. everything exactly *)
  v = TimeConstrained[derivZeroVec[anti, rest, t, y, R, n], vecBudget[], $Failed];
  TrueQ[v]];

derivZeroVec[anti_, rest_, t_, y_, R_, n_] := Module[{terms, vec, left = 0, cv},
  terms = If[Head[anti] === Plus, List @@ anti, {anti}];
  vec = vecOf[rest, R, t, y, n];
  Do[cv = Quiet[Catch[derivVec[term, t, y, R, n], "GoursatFail"]];
     If[ListQ[cv], vec = addV[vec, cv], left += term],
     {term, terms}];
  If[left === 0, Return[AllTrue[vec, zeroEntryQ[#, t] &], Module]];
  pzqTrue[(Together[totalD[left, t, y, R, n]] + vec.(y^Range[0, n - 1])) /. y -> R^(1/n), t, 60]];

(* The exact zero test of a vector entry (a rational function of t over       *)
(* Q(algebraic numbers), not canonicalised): every coefficient of its expanded  *)
(* numerator is an algebraic number, decided by PossibleZeroQ (exact on         *)
(* algebraic numbers); anything else by the canonical form gK.                  *)
zeroEntryQ[e_, t_] := Module[{num = Expand[Numerator[Together[e]]], cs},
  If[num === 0, Return[True, Module]];
  cs = If[PolynomialQ[num, t], CoefficientList[num, t], {num}];
  If[AllTrue[cs, NumericQ],
    AllTrue[cs, numZeroQ],
    gK[e] === 0]];

(* An algebraic number is zero: its canonical form, else PossibleZeroQ. *)
numZeroQ[c_] := c === 0 || gK[c] === 0 || TrueQ[Quiet[PossibleZeroQ[c]]];

unevaluatedQ[anti_] := transcendentalQ[anti /. {_Log -> 0, _ArcTan -> 0, _ArcTanh -> 0,
                                                 _ArcCot -> 0, _ArcCoth -> 0, _RootSum -> 0}];

verifyResult[res_, F_, t_] := Module[{y = res["y"], n = res["n"], k = res["k"], R = res["R"]},
  If[unevaluatedQ[res["Antiderivative"]], Return[False, Module]];
  derivZeroQ[res["Antiderivative"], res["Residual"]/y^k - F/y^k, t, y, R, n]];

verifyExprResult[res_, f_, t_] := Module[{fy},
  If[unevaluatedQ[res["Antiderivative"]], Return[False, Module]];
  If[!KeyExistsQ[res, "y"],
    Return[TrueQ[Quiet[PossibleZeroQ[Together[D[res["Antiderivative"], t] + res["Residual"] - f]]]],
           Module]];
  fy = f /. Power[b_, q_Rational] /; !FreeQ[b, t] :> res["y"]^(q res["n"]);
  derivZeroQ[res["AntiderivativeY"], res["ResidualY"] - fy, t, res["y"], res["R"], res["n"]]];
