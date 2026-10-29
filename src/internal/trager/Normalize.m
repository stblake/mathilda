(* ::Package:: *)

(* ========================================================================== *)
(* Normalize.m -- Hermite normal form over K[x] and normalisation of a module *)
(*                basis at infinity.                                          *)
(*                                                                            *)
(* HNF ([T] Ch 2 p.25).  A K[x]-submodule of K(x)^m given by polynomial       *)
(* generating rows is put in upper-triangular Hermite form by unimodular row  *)
(* operations built from the extended Euclidean algorithm: for a pivot a and  *)
(* an entry b below it, with s a + t b = g, the pair of rows (R1, R2) becomes *)
(*     (s R1 + t R2,  (a/g) R2 - (b/g) R1),                                   *)
(* a transformation of determinant 1.  Pivots are made monic and the entries  *)
(* above a pivot are reduced modulo it.                                       *)
(*                                                                            *)
(* Normal at infinity ([T] Ch 2 Sec 3 pp.27-29).  The basis w_i is normal at  *)
(* infinity with local integral basis v_j = x^(-delta_j) w_j (Basis.m).  A    *)
(* row (c_1, ..., c_m)/den represents Sum c_j w_j / den = Sum (c_j x^delta_j  *)
(* / den) v_j; its order at infinity in the x-adic valuation of K(x) is       *)
(*     k = min_j ord_inf(c_j x^delta_j / den) = min_j -(deg c_j + delta_j     *)
(*                                                   - deg den),              *)
(* and its leading vector N_j = lc(c_j)/lc(den) where that minimum is         *)
(* attained, 0 elsewhere.  The basis is normal iff the leading vectors are    *)
(* independent over K.  Otherwise take c with Sum_i c_i N_i = 0 and i0 of     *)
(* least k among c_i != 0, and replace row i0 by Sum_i c_i x^(k_i - k_i0)      *)
(* row_i: the module is unchanged (the coefficient of row i0 is a nonzero     *)
(* constant, the others polynomials), and k_i0 strictly increases.  All       *)
(* exponents are integers -- the package's fractional variant could leave     *)
(* the module.                                                                *)
(* ========================================================================== *)

hnf::usage = "hnf[F, rows, x] is the Hermite normal form over F[x] of the K[x]-module spanned by the polynomial row vectors rows: an upper-triangular m x m matrix (m = number of columns) with monic diagonal and entries above each pivot of lower degree. The module must have full rank m; otherwise this is an internal error.";
hnf[F_Association, rows0_List, x_Symbol] := Module[
  {M, nr, nc, piv, a, b, eg, g, s, t, r1, r2, q, inv, nat},
  If[nativeQ[],
    nat = FLINT`PolynomialHermiteForm[rows0, x, Sequence @@ nativeOpts[F]];
    If[ListQ[nat], Return[nat]]];
  M = Map[kc[F, #] &, rows0, {2}];
  nr = Length[M]; nc = Length[First[M]];
  If[nr < nc, throwInternal["hnfRankDeficient"]];
  Do[
    (* bring the nonzero entry of least degree to the pivot position *)
    piv = SelectFirst[SortBy[Select[Range[j, nr], M[[#, j]] =!= 0 &], Exponent[M[[#, j]], x] &], True &, 0];
    If[piv === 0, throwInternal["hnfRankDeficient", <|"Column" -> j|>]];
    (* {M[[j]], M[[piv]]} = {...} is a silent no-op in Mathilda: swap by a list of positions *)
    If[piv =!= j, M[[{j, piv}]] = M[[{piv, j}]]];
    Do[
      If[M[[i, j]] =!= 0,
        a = M[[j, j]]; b = M[[i, j]];
        {g, s, t} = kExtGCD[F, a, b, x];
        r1 = M[[j]]; r2 = M[[i]];
        M[[j]] = kc[F, #] & /@ (s r1 + t r2);
        M[[i]] = kc[F, #] & /@ (kExactQuo[F, a, g, x] r2 - kExactQuo[F, b, g, x] r1)],
      {i, j + 1, nr}];
    inv = kInv[F, kLC[F, M[[j, j]], x]];
    M[[j]] = kc[F, # inv] & /@ M[[j]];
    Do[
      If[M[[i, j]] =!= 0 && Exponent[M[[i, j]], x] >= Exponent[M[[j, j]], x],
        q = kQuo[F, M[[i, j]], M[[j, j]], x];
        M[[i]] = MapThread[kc[F, #1 - q #2] &, {M[[i]], M[[j]]}]],
      {i, 1, j - 1}],
    {j, nc}];
  If[!AllTrue[Flatten[M[[nc + 1 ;;]]], # === 0 &], throwInternal["hnfResidualRows"]];
  Take[M, nc]];

hnfModDet::usage = "hnfModDet[F, rows, x, D] is hnf[F, rows, x] when D is a multiple of the determinant of the module spanned by rows (a full-rank lattice L contains D K[x]^m as soon as det L divides D, by the adjugate). The rows D e_k are adjoined and every non-pivot entry is reduced modulo D after each row operation, so entries stay of degree < deg D instead of growing through the elimination. The result is the same Hermite normal form.";
hnfModDet[F_Association, rows0_List, x_Symbol, D_] := Module[
  {m = Length[First[rows0]], Dm = kMonic[F, D, x], M, nr, piv, a, b, g, s, t, r1, r2, q, inv, red, nat},
  If[kDeg[F, Dm, x] <= 0, Return[hnf[F, rows0, x]]];
  If[nativeQ[],
    nat = FLINT`PolynomialHermiteForm[rows0, x, "DeterminantMultiple" -> Dm, Sequence @@ nativeOpts[F]];
    If[ListQ[nat], Return[nat]]];
  red[row_, j_] := MapIndexed[If[First[#2] > j, kRem[F, #1, Dm, x], kc[F, #1]] &, row];
  M = Join[Map[kRem[F, #, Dm, x] &, rows0, {2}], Dm IdentityMatrix[m]];
  M = Map[kc[F, #] &, M, {2}];
  nr = Length[M];
  Do[
    piv = SelectFirst[SortBy[Select[Range[j, nr], M[[#, j]] =!= 0 &], Exponent[M[[#, j]], x] &], True &, 0];
    If[piv === 0, throwInternal["hnfRankDeficient", <|"Column" -> j|>]];
    If[piv =!= j, M[[{j, piv}]] = M[[{piv, j}]]];
    Do[
      If[M[[i, j]] =!= 0,
        a = M[[j, j]]; b = M[[i, j]];
        {g, s, t} = kExtGCD[F, a, b, x];
        r1 = M[[j]]; r2 = M[[i]];
        M[[j]] = red[s r1 + t r2, j];
        M[[i]] = red[kExactQuo[F, a, g, x] r2 - kExactQuo[F, b, g, x] r1, j]],
      {i, j + 1, nr}];
    inv = kInv[F, kLC[F, M[[j, j]], x]];
    M[[j]] = kc[F, # inv] & /@ M[[j]],
    {j, m}];
  M = Take[M, m];
  (* reduce above the pivots, as hnf does *)
  Do[
    Do[
      If[M[[i, j]] =!= 0 && Exponent[M[[i, j]], x] >= Exponent[M[[j, j]], x],
        q = kQuo[F, M[[i, j]], M[[j, j]], x];
        M[[i]] = MapThread[kc[F, #1 - q #2] &, {M[[i]], M[[j]]}]],
      {i, 1, j - 1}],
    {j, m}];
  M];

rowInfData::usage = "rowInfData[F, B, row, den] is {k, N}: the x-adic order k at infinity of the element (row/den) on the local integral basis v_j = x^(-delta_j) w_j, and its leading vector N (the leading coefficients of the coordinates of minimal order, 0 elsewhere). k is Infinity for a zero row.";
rowInfData[F_Association, B_Association, row_List, den_] := Module[
  {x = B["x"], dd, lcd, ords, k},
  dd = Exponent[den, x];
  lcd = kLC[F, den, x];
  ords = Table[
    If[row[[j]] === 0, Infinity, -(Exponent[row[[j]], x] + B["delta"][[j]] - dd)],
    {j, Length[row]}];
  k = Min[ords];
  {k, Table[If[ords[[j]] === k && k =!= Infinity,
      kc[F, kLC[F, row[[j]], x] kInv[F, lcd]], 0], {j, Length[row]}]}];

normalizeAtInfinity::usage = "normalizeAtInfinity[F, B, M, den, x] normalises at infinity ([T] Ch 2 Sec 3) the basis rows M (polynomial over F[x]) of the module M/den of elements on the integral basis B. It returns <|\"rows\" -> M', \"k\" -> {k_i}|>: a basis of the same module whose leading vectors at infinity are independent, so that an element Sum p_i row_i/den (p_i in F[x]) has no pole at infinity iff every p_i row_i/den has none, and the orders k_i of the rows at infinity.";
normalizeAtInfinity[F_Association, B_Association, M0_List, den_, x_Symbol] := Module[
  {M = M0, n = Length[M0], data, ks, Nm, c, nz, i0, newRow, it = 0, cap},
  data = rowInfData[F, B, #, den] & /@ M;
  (* Every step raises one k_i by at least 1, and Sum k_i is bounded by the   *)
  (* order at infinity of the determinant ([T] p.28), so this cap is never   *)
  (* reached by a correct run: it guards against a non-terminating bug.      *)
  cap = 4 (Total[Abs[Exponent[#, x]] & /@ Select[Flatten[M0], # =!= 0 &]] +
           n Abs[Exponent[den, x]] + Total[Abs[B["delta"]]] + n) + 200;
  While[True,
    If[++it > cap, throwInternal["normalizeAtInfinityCap"]];
    ks = First /@ data; Nm = Last /@ data;
    If[MemberQ[ks, Infinity], throwInternal["normalizeAtInfinityZeroRow"]];
    c = kLeftNullVector[F, Nm];
    If[c === None, Break[]];
    nz = Select[Range[n], c[[#]] =!= 0 &];
    i0 = First[SortBy[nz, ks[[#]] &]];
    newRow = Sum[c[[i]] x^(ks[[i]] - ks[[i0]]) M[[i]], {i, nz}];
    M[[i0]] = kc[F, #] & /@ newRow;
    data[[i0]] = rowInfData[F, B, M[[i0]], den];
    If[data[[i0, 1]] <= ks[[i0]], throwInternal["normalizeAtInfinityNoProgress"]]];
  <|"rows" -> M, "k" -> ks|>];
