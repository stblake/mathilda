(* ::Package:: *)

(* ========================================================================== *)
(* Jacobian.m -- divisor classes in reduced form: constant-size class         *)
(*               arithmetic for the torsion order and the principal           *)
(*               generator.                                                   *)
(*                                                                            *)
(* The order of a residue divisor delta (TorsionModP.m) is found by testing   *)
(* M delta for M = 1, 2, ..., and the logand is a generator of N delta.  With *)
(* ideals of M delta the sizes grow linearly in M, which is what made the     *)
(* large-order cases time out.  Divisor classes have bounded representatives: *)
(*                                                                            *)
(* Frame.  Put infinity at a point x = a where G has a root of multiplicity e *)
(* prime to m:  x = a + 1/t, Yt = Y t^k (k = Ceiling[M/m]), so                 *)
(*     Yt^m = Gt(t) = t^(k m - M) t^M G(a + 1/t),   deg Gt = k m - e.          *)
(* deg Gt is prime to m, so over t = Infinity there is ONE place P, totally    *)
(* ramified and of degree 1, and for c(t) w_j                                 *)
(*     ord_P(c w_j) = -m deg c + m deg d_j - j deg Gt                         *)
(* exactly (the values are distinct modulo m for distinct j, so in a sum the  *)
(* order is the least order of the terms).                                    *)
(*                                                                            *)
(* Reduced representatives.  Every class of degree 0 is [A - deg(A) P] with A *)
(* effective of degree s <= genus: for D = D_fin + n P of degree 0 take the   *)
(* least s with L(D + s P) != 0; an h in it gives                             *)
(*     D + div(h) = A - s P,   A = div_fin(h) + D_fin >= 0                     *)
(* (Riemann-Roch: s <= genus).  The class is zero iff s = 0.  L(D + s P) is   *)
(* read off a basis of I(-D_fin) in which the leading terms at P are in       *)
(* distinct columns (a Popov form, reducedRows): its element of greatest      *)
(* order at P has order -(n + s).                                             *)
(*                                                                            *)
(* A class is stored as <|"J" -> I(-A), "d" -> deg A, "H" -> {{h, e}, ...}|>  *)
(* with  k delta + div(Prod h^e) = A - d P  (the "H" list -- the principal    *)
(* factors met on the way, a compact representation of the generator -- is   *)
(* kept only when asked).  Sums multiply the J ideals (I(-(A1 + A2)) =        *)
(* I(-A1) I(-A2)) and reduce; no ideal inverse is needed.                     *)
(* ========================================================================== *)

curveBasisOver::usage = "curveBasisOver[F, m, G, t] is the integral basis descriptor of Y^m = G(t) over the field F: curveBasis over Q (and over a number field, the curve being defined over Q), its analogue from the factorisation over GF(p) otherwise.";
curveBasisOver[F_Association, m_Integer, G_, t_Symbol] := If[F["Type"] =!= "GF",
  curveBasis[m, G, t],
  Module[{p = F["p"], Gp = gfMod[G, F["p"]], pairs, d, s, M, delta, mult},
    pairs = Select[FactorList[Gp, Modulus -> p], !FreeQ[#[[1]], t] &];
    d = Table[gfMod[Times @@ (#[[1]]^Floor[i #[[2]]/m] & /@ pairs), p], {i, 0, m - 1}];
    s = Table[rc[F, i D[Gp, t]/(m Gp) - D[d[[i + 1]], t]/d[[i + 1]]], {i, 0, m - 1}];
    M = Exponent[Gp, t];
    delta = Table[Ceiling[i M/m] - Exponent[d[[i + 1]], t], {i, 0, m - 1}];
    mult = Table[
      With[{k = Mod[i + j, m], q = Quotient[i + j, m]},
        {k, gfMod[PolynomialQuotient[gfMod[Gp^q d[[k + 1]], p], gfMod[d[[i + 1]] d[[j + 1]], p], t, Modulus -> p], p]}],
      {i, 0, m - 1}, {j, 0, m - 1}];
    <|"m" -> m, "G" -> Gp, "x" -> t, "deg" -> M, "factors" -> pairs,
      "d" -> d, "s" -> s, "delta" -> delta, "mult" -> mult, "p" -> p|>]];

rootMultiplicity::usage = "rootMultiplicity[F, G, a, x] is the multiplicity of the root a of the polynomial G over F (0 when G(a) != 0).";
rootMultiplicity[F_Association, G_, a_, x_Symbol] := Module[{g = kc[F, G], e = 0, qr},
  While[g =!= 0 && kc[F, g /. x -> a] === 0,
    qr = kDivRem[F, g, x - a, x]; g = qr[[1]]; e++];
  e];

frameRoots::usage = "frameRoots[F, G, m, x] lists the roots a in F of the polynomial G (over Q: rational roots; over GF(p): roots mod p) whose multiplicity is prime to m: the points at which the reduced frame can be centred.";
frameRoots[F_Association, G_, m_Integer, x_Symbol] := Module[{lin},
  lin = Select[
    If[F["Type"] === "GF", FactorList[gfMod[G, F["p"]], Modulus -> F["p"]], FactorList[G]],
    Exponent[#[[1]], x] === 1 && GCD[#[[2]], m] === 1 &];
  kc[F, -Coefficient[#[[1]], x, 0] kInv[F, Coefficient[#[[1]], x, 1]]] & /@ lin];

jacFrame::usage = "jacFrame[F, m, G, c, x, a] is the reduced frame centred at x = a (a root of G of multiplicity prime to m) for the remainder c(x) Y dx on Y^m = G: <|\"B\" -> basis in t, \"N\", \"E\", \"Ep\", \"Bp\" (the remainder c_t(t) Yt dt = (N/E) Yt dt and its residue data in t), \"a\", \"k\", \"t\", \"x\"|>. Residues are invariant under the change of variable, so residue ideals in the frame use the same residues.";
jacFrame[F_Association, m_Integer, G_, c_, x_Symbol, a_] := Module[{t = Unique["t"], M, k, Gt, Bt, ct, nd, E, g},
  M = Exponent[kc[F, G], x];
  k = Ceiling[M/m];
  Gt = kc[F, t^(k m - M) Expand[t^M (G /. x -> a + 1/t)]];
  Bt = curveBasisOver[F, m, Gt, t];
  ct = rc[F, -(c /. x -> a + 1/t) t^(-k - 2)];
  nd = rcNumDen[F, ct];
  E = kMonic[F, nd[[2]], t];
  g = kGCD[F, E, Bt["G"], t];
  <|"B" -> Bt, "N" -> kc[F, nd[[1]] kInv[F, kLC[F, nd[[2]], t]]], "E" -> E, "Ep" -> kc[F, D[E, t]],
    "Bp" -> kExactQuo[F, E, g, t], "a" -> a, "k" -> k, "t" -> t, "x" -> x|>];

frameToX::usage = "frameToX[fr, expr, Yt, Y] maps an expression in the frame variable t and generator Yt back to x and Y: t -> 1/(x - a), Yt -> Y/(x - a)^k.";
frameToX[fr_Association, expr_, Yt_Symbol, Y_Symbol] := With[{x = fr["x"], a = fr["a"]},
  expr /. {fr["t"] -> 1/(x - a), Yt -> Y (x - a)^(-fr["k"])}];

placeOrder::usage = "placeOrder[F, Bt, row, den] is the exact order at the single place P over t = Infinity of the element row/den of the frame basis Bt (Infinity for zero).";
placeOrder[F_Association, Bt_Association, row_List, den_] := Module[{t = Bt["x"], m = Bt["m"], dd, ords},
  dd = kDeg[F, den, t];
  ords = Table[
    If[row[[j + 1]] === 0, Infinity,
      -m (Exponent[row[[j + 1]], t] - dd) + m Exponent[Bt["d"][[j + 1]], t] - j Bt["deg"]],
    {j, 0, m - 1}];
  Min[ords]];

leadingColumn::usage = "leadingColumn[F, Bt, row] is the column (0-based) where the order at P of the polynomial row is attained (unique: the orders in different columns differ modulo m).";
leadingColumn[F_Association, Bt_Association, row_List] := Module[{t = Bt["x"], m = Bt["m"], best = Infinity, col = None, o},
  Do[
    If[row[[j + 1]] =!= 0,
      o = -m Exponent[row[[j + 1]], t] + m Exponent[Bt["d"][[j + 1]], t] - j Bt["deg"];
      If[o < best, best = o; col = j]],
    {j, 0, m - 1}];
  col];

reducedRows::usage = "reducedRows[F, Bt, rows] brings the basis rows (polynomial over F[t]) of a module to a form whose leading terms at P lie in pairwise distinct columns (a Popov form), by repeatedly cancelling the leading term of the row of lower order against a t-multiple of a row with the same leading column. The module is unchanged; afterwards the order at P of Sum p_i row_i is the least order of its terms. Returns the rows.";
reducedRows[F_Association, Bt_Association, rows0_List] := Module[
  {t = Bt["x"], R = rows0, n = Length[rows0], cols, ords, clash, i, j, lo, hi, s, f, it = 0, cap},
  cap = 4 (Total[Exponent[#, t] & /@ Select[Flatten[rows0], # =!= 0 &]] + n) + 200;
  While[True,
    If[++it > cap, throwInternal["reducedRowsCap"]];
    cols = leadingColumn[F, Bt, #] & /@ R;
    clash = SelectFirst[Subsets[Range[n], {2}], cols[[#[[1]]]] === cols[[#[[2]]]] &, None];
    If[clash === None, Break[]];
    {i, j} = clash;
    ords = placeOrder[F, Bt, #, 1] & /@ R[[{i, j}]];
    (* lo: the row of lower order (larger pole); cancel its leading term *)
    {lo, hi} = If[ords[[1]] <= ords[[2]], {i, j}, {j, i}];
    s = (Max[ords] - Min[ords])/Bt["m"];
    With[{c = cols[[i]] + 1},
      f = kc[F, kLC[F, R[[lo, c]], t] kInv[F, kLC[F, R[[hi, c]], t]]];
      R[[lo]] = kc[F, #] & /@ (R[[lo]] - f t^s R[[hi]])]];
  R];

idealNormDegree::usage = "idealNormDegree[F, B, I] is the degree of the divisor of the fractional ideal I = M/den: deg det M - m deg den.";
idealNormDegree[F_Association, B_Association, I_Association] :=
  kDeg[F, Times @@ Diagonal[I["M"]], B["x"]] - B["m"] kDeg[F, I["den"], B["x"]];

elementInverse::usage = "elementInverse[F, B, h] is {num, den} with 1/h = num/den: num the first row of the adjugate of the multiplication matrix of h times the norm cofactor, den the norm of h over the prime field (h with polynomial coordinates).";
elementInverse[F_Association, B_Association, h_List] := Module[{x = B["x"], M, det, adj1, nrm, cof},
  M = Table[afTimesPoly[F, B, afUnit[B, k], h], {k, 0, B["m"] - 1}];
  det = kc[F, Det[M]];
  If[det === 0, throwInternal["elementInverseZero"]];
  adj1 = If[B["m"] === 1, {1}, Table[kc[F, (-1)^(1 + j) Det[Drop[M, {j}, {1}]]], {j, B["m"]}]];
  nrm = kMonic[F, kNorm[F, det], x];
  cof = kExactQuo[F, nrm, det, x];
  {kc[F, # cof] & /@ adj1, nrm}];

jacReduce::usage = "jacReduce[F, Bt, J, n] reduces the degree-0 divisor D = D_fin + n P, given J = I(-D_fin): it returns {Jnew, s, h} with D + div(h) = A - s P, A >= 0 of degree s, Jnew = I(-A) = h^(-1) J, and h (coordinates over F[t], divided by J's denominator) the element of greatest order at P in J.";
jacReduce[F_Association, Bt_Association, J_Association, n_Integer] := Module[
  {t = Bt["x"], rows, ords, best, h, s, inv, newRows, Jn},
  rows = reducedRows[F, Bt, J["M"]];
  ords = placeOrder[F, Bt, #, J["den"]] & /@ rows;
  best = First[Ordering[ords, -1]];
  h = rc[F, #/J["den"]] & /@ rows[[best]];
  s = -ords[[best]] - n;
  If[s < 0, throwInternal["jacReduceNegativeDegree", <|"s" -> s|>]];
  (* Jnew = h^(-1) J: rows of J times the adjugate of the polynomial part of h *)
  inv = elementInverse[F, Bt, kc[F, #] & /@ rows[[best]]];
  newRows = afTimesPoly[F, Bt, #, inv[[1]]] & /@ J["M"];
  Jn = idealReduce[F, Bt, <|"M" -> hnf[F, newRows, t], "den" -> inv[[2]]|>];
  {Jn, s, h}];

jacClass::usage = "jacClass[F, Bt, J, n, track] is the reduced class of the degree-0 divisor D_fin + n P with J = I(-D_fin): <|\"J\" -> I(-A), \"d\" -> deg A, \"H\" -> {{h, 1}} or {}|> (H kept when track is True).";
jacClass[F_Association, Bt_Association, J_Association, n_Integer, track_ : False] := Module[{r = jacReduce[F, Bt, J, n]},
  <|"J" -> r[[1]], "d" -> r[[2]], "H" -> If[TrueQ[track], {{r[[3]], 1}}, {}]|>];

jacWithGenerators::usage = "jacWithGenerators[F, Bt, c] adds to the class c a verified two-generator form of its J ideal (key \"gens\"), used by jacAdd when c is the second summand: worth it for an addend used repeatedly (the stepping in jacOrder).";
jacWithGenerators[F_Association, Bt_Association, c_Association] := Module[{g = twoGenerators[F, Bt, c["J"]]},
  If[g === None, c, Join[c, <|"gens" -> g|>]]];

jacZeroQ::usage = "jacZeroQ[c] is True for the zero class (reduced degree 0).";
jacZeroQ[c_Association] := c["d"] === 0;

jacAdd::usage = "jacAdd[F, Bt, c1, c2, track] is the reduced class of c1 + c2: the J ideals multiply (I(-(A1 + A2)) = I(-A1) I(-A2)) and the sum A1 + A2 - (d1 + d2) P is reduced. With track, H collects the principal factors: (k1 + k2) delta + div(Prod h^e) = A - d P.";
jacAdd[F_Association, Bt_Association, c1_Association, c2_Association, track_ : False] := Module[{J, r, H, g},
  (* a two-generator form of the second summand, when it carries one       *)
  (* (jacWithGenerators), makes the product 2m element products, not m^2   *)
  g = Lookup[c2, "gens", None];
  J = If[g === None, idealTimes[F, Bt, c1["J"], c2["J"]],
    idealTimesGenerators[F, Bt, c1["J"], {g[[1]], c2["J"]["den"]}, g[[2]]]];
  r = jacReduce[F, Bt, J, -(c1["d"] + c2["d"])];
  H = If[TrueQ[track], mergeFactors[Join[c1["H"], c2["H"], {{r[[3]], 1}}]], {}];
  <|"J" -> r[[1]], "d" -> r[[2]], "H" -> H|>];

mergeFactors::usage = "mergeFactors[H] combines the entries {h, e} of H with identical h, dropping exponent 0.";
mergeFactors[H_List] := Select[{#[[1, 1]], Total[#[[All, 2]]]} & /@ GatherBy[H, First], #[[2]] =!= 0 &];

jacMultiple::usage = "jacMultiple[F, Bt, c, N, track] is the class N c (N >= 1) by doubling and adding.";
jacMultiple[F_Association, Bt_Association, c_Association, N_Integer?Positive, track_ : False] := Module[
  {bits = IntegerDigits[N, 2], acc},
  acc = c;
  Do[
    acc = jacAdd[F, Bt, acc, acc, track];
    If[b === 1, acc = jacAdd[F, Bt, acc, c, track]],
    {b, Rest[bits]}];
  acc];

jacKey::usage = "jacKey[c] is a canonical key of the class c: with a single place P of degree 1 at infinity and a Popov-reduced basis, the space L(D + s P) at the least s is one-dimensional (only one basis row can attain the greatest order, the orders of different columns being distinct modulo m), so the reduced representative A is unique and the Hermite form of I(-A) identifies the class.";
jacKey[c_Association] := {c["d"], c["J"]["M"], c["J"]["den"]};

jacOrder::usage = "jacOrder[F, Bt, c, bound] is the least M >= 1 with M c = 0, M <= bound (a Decline otherwise): baby-step giant-step on the canonical keys (about 2 Sqrt[bound] class additions), then the least divisor of the multiple found that still annihilates c.";
jacOrder[F_Association, Bt_Association, c0_Association, bound_Integer] := Module[
  {c = jacWithGenerators[F, Bt, c0]},
  jacOrderGens[F, Bt, c, bound]];

jacOrderGens::usage = "jacOrderGens[F, Bt, c, bound] is jacOrder for a class already carrying its two generators.";
jacOrderGens[F_Association, Bt_Association, c_Association, bound_Integer] := Module[
  {k, baby, acc, step, M0 = None, i = 0, jj = 1, key, found = None},
  If[jacZeroQ[c], Return[1]];
  (* Small orders are the common case: step linearly first. *)
  acc = c;
  Do[
    acc = jacAdd[F, Bt, acc, c];
    If[jacZeroQ[acc], found = jj; Break[]],
    {jj, 2, Min[bound, 24]}];
  If[found =!= None, Return[found]];
  If[bound <= 24, throwDecline["torsion order search bound reached", <|"Bound" -> bound|>]];
  k = Ceiling[Sqrt[N[bound]]] + 1;
  (* baby steps: j c for j = 1..k; a zero on the way is the order *)
  baby = <||>; acc = c;
  While[found === None && jj <= k,
    If[jacZeroQ[acc], found = jj,
      key = jacKey[acc];
      If[!KeyExistsQ[baby, key], baby[key] = jj];
      If[jj < k, acc = jacAdd[F, Bt, acc, c]];
      jj++]];
  If[found =!= None, Return[found]];
  (* giant steps: i k c; a match with j c gives (i k - j) c = 0 *)
  step = jacWithGenerators[F, Bt, acc]; acc = step; i = 1;
  While[M0 === None,
    If[i k > bound + k, throwDecline["torsion order search bound reached", <|"Bound" -> bound|>]];
    key = jacKey[acc];
    Which[
      jacZeroQ[acc], M0 = i k,
      KeyExistsQ[baby, key] && i k - baby[key] > 0, M0 = i k - baby[key],
      True, i++; acc = jacAdd[F, Bt, acc, step]]];
  jacMinimalOrder[F, Bt, c, M0]];

jacMinimalOrder::usage = "jacMinimalOrder[F, Bt, c, M0] is the order of c given that M0 c = 0: M0 with every prime factor removed while the result still annihilates c.";
jacMinimalOrder[F_Association, Bt_Association, c_Association, M0_Integer] := Module[{M = M0},
  Do[
    While[Mod[M, pr] === 0 && M > 1 && jacZeroQ[jacMultiple[F, Bt, c, M/pr]], M = M/pr],
    {pr, First /@ FactorInteger[M0]}];
  M];
