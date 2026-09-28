(* logrewrite.m -- the real form of the logarithmic part (Rioboo).

   The Mathilda port of logrewrite.wl of the research directory (itself the counterpart of
   logrewrite.py and maxima/logrewrite.mac), one function per function; loaded by
   ParallelMixed.m through LoadModule["mixed/logrewrite.m"] inside its private context.
   A post-processing layer, not part of the parallel algorithm: the integrals come out of
   the residue analysis and the ansatz as a rational part plus Sum_k c_k Log[u_k] with
   algebraic constants c_k, and a real integrand gives its non-real terms in conjugate
   pairs (c, u), (conj c, conj u).  Every such pair is collapsed to a real logarithm and
   arctangents, and every pair conjugate under Sqrt[p] -> -Sqrt[p] to a real logarithm
   and a hyperbolic arctangent where the argument is a polynomial.  The one hook is
   LogToReal, called by iPIM on the terms in the tower coordinates, before the generators
   are substituted back; the RootSum terms of the Root-object primes are not touched.

   Rioboo (Bronstein, Symbolic Integration I, Section 2.8): for a conjugate pair

     (a + I b) Log[A + I B] + (a - I b) Log[A - I B]  =  a Log[A^2 + B^2] + b LogToAtan[A, B],

   LogToAtan[A, B] a sum of arctangents of POLYNOMIALS when A, B are polynomials in the
   top generator over the field of the lower ones (the Euclidean algorithm of LogToAtan,
   so the result does not jump at a zero of B); when A, B involve the radical, one
   arctangent of the quotient A/B or -B/A, the denominator chosen free of real zeros when
   the norm shows it (Sturm's theorem), else A/B, the quotient written as a polynomial in the
   radical over the generators (YQuot: the surface layer rewrites the positive powers of the
   radical only, so a radical in a denominator would be read on the principal branch).  The
   pairing is formal (I -> -I with every other symbol fixed), so the derivative identity
   holds for complex parameters too; logands sharing a factor with the conjugate of another
   are split first and equal logands are merged.  When the problem is real the formal
   conjugate of an antiderivative is an antiderivative, and a non-real term left without a
   partner is collapsed with half its coefficient, as the average of the two -- its partner
   is then an imaginary multiple of a generator in the polynomial part (ArcSin[x]/x^2), or a
   logand conjugate to it only on the curve (I x + Sqrt[1 - x^2]) -- and the polynomial part
   is taken real.  The same identity with I replaced by Sqrt[p], p a prime, gives the
   hyperbolic arctangent: a pair (c, u), (c, u)|_{Sqrt p -> -Sqrt p} with u = A + Sqrt[p] B
   is a Log[A^2 - p B^2] + 2 b Sqrt[p] ArcTanh[A/(Sqrt[p] B)] when B divides A (the argument
   a polynomial, as in ArcTanh[x/Sqrt[2]] for 1/(x^2 - 2)); otherwise the two real
   logarithms are kept.  Two safety valves leave a sum as it came: constants that are not
   tame (TameQ: a square root of a non-real number with irrational parts, a higher root of a
   non-rational base, a Root object without a radical form, a trigonometric constant) and a
   budget of $RewriteBudget seconds (TimeConstrained) on the rewrite as a whole.

   Where this file differs from logrewrite.wl (the kernel divergences are catalogued in
   MATHILDA_DIVERGENCES.md, A19-A23):
     - Can leaves nested radicals in non-canonical forms: SameQ0 settles a numerically zero
       difference coefficient by coefficient (RootReduce), and a square root of a non-real
       number with real algebraic parts and a rational squared modulus counts as tame
       (RationalModulusQ), since the places of this kernel arrive in that form;
     - FreeQ[e, Complex] is True on a complex atom here, so every test for I is written
       with the pattern _Complex (and _Root, _Re, ... for the heads);
     - ComplexExpand writes a real nested radical as Cos[Arg[...]] terms and leaves
       Arg[1 - Sqrt[5]] unevaluated, so RectConst is rule-based (RectPow: the sign of a real
       base decides, the principal square root of a non-real one is written by hand,
       RootReduce denests where it can) and ComplexExpand is applied only to (-1)^k;
     - CountRoots is not implemented: SturmCount (Sturm's theorem) counts the real roots;
     - NumericQ is False on a Root object and ToRadicals writes every quartic root in the
       Ferrari form: RootRadicals finds the quadratic-in-disguise forms (biquadratic,
       palindromic) that Mathematica's ToRadicals returns itself, and the Root objects are
       replaced before any NumericQ test;
     - a Do iterator captures a caller's symbol of the same name: the iterators are
       lr-prefixed;
     - an outer TimeConstrained cannot interrupt an inner one, so the package budget
       cannot cut the rewrite short; the rewrite's own budget can.                         *)

(* CanL[e]: the canonical form of the package -- CanRaw, Cancel over the extension of the
   constants, which is what Can is in ParallelMixed.wl; the number-field detour of this
   kernel's Can (FieldData) is bypassed here: on the nested-radical constants of the rewrite it
   is slow and it has returned a Dot[{}, Inverse[{}], {}] coefficient (A40 after P4) -- with
   any Root object written back in radicals, so that the formal split I -> -I can read the
   constants *)
CanL[e_] := With[{r = CanRaw[e]}, If[FreeQ[r, _Root], r, RootFree[r]]];

(* RootFree[e]: e with every Root object in radicals where it has a usable radical form *)
RootFree[e_] := e /. rr_Root :> RootRadicals[rr];

(* QRoots[a, b, c]: the two roots of a z^2 + b z + c, the discriminant expanded *)
QRoots[a_, b_, c_] := With[{dd = Sqrt[Expand[b^2 - 4 a c]]}, {(-b + dd)/(2 a), (-b - dd)/(2 a)}];

(* RootRadicals[r]: a Root object in radicals: the kernel's ToRadicals when its radicals are
   tame, else for a quartic the quadratic-in-disguise forms -- biquadratic z^4 + a z^2 + b,
   palindromic (w = z + 1/z) and antipalindromic (w = z - 1/z) -- with the root picked by its
   30-digit value (Mathematica's ToRadicals finds these forms itself; the kernel's gives the
   Ferrari form, which is not tame); the Root object itself when nothing applies *)
RootRadicals[r_Root] := Module[{tr, z, pol, cl, cands = {}, num},
  tr = ToRadicals[r];
  If[FreeQ[tr, _Root] && RadicalsTameQ[tr], Return[tr]];
  pol = Expand[r[[1]][z]];
  If[Exponent[pol, z] != 4, Return[r]];
  cl = CoefficientList[pol, z]/Coefficient[pol, z, 4];                (* {c0, c1, c2, c3, 1} *)
  Which[
    cl[[2]] == 0 && cl[[4]] == 0,
      cands = Flatten[{Sqrt[#], -Sqrt[#]} & /@ QRoots[1, cl[[3]], cl[[1]]]],
    cl[[1]] == 1 && cl[[2]] == cl[[4]],
      cands = Flatten[QRoots[1, -#, 1] & /@ QRoots[1, cl[[2]], cl[[3]] - 2]],
    cl[[1]] == 1 && cl[[2]] == -cl[[4]],                              (* z^4 + a z^3 + b z^2 - a z + 1: w^2 + a w + b + 2 = 0 *)
      cands = Flatten[QRoots[1, -#, -1] & /@ QRoots[1, cl[[4]], cl[[3]] + 2]]];
  If[cands === {}, Return[r]];
  num = N[r, 30];
  cands = Select[cands, Abs[N[#, 30] - num] < 10^-20 &];
  If[cands === {}, r, First[cands]]];

(* RadicalsTameQ[e]: the radical test of TameQ on an expression free of Root objects.  Beyond
   the rule of logrewrite.wl, a square root of a non-real number with real algebraic parts and
   a RATIONAL squared modulus is tame here (Sqrt[(1 - I Sqrt[3])/2] = (Sqrt[3] - I)/2): the
   places of this kernel arrive in that form where Mathematica's are already rectangular, and
   RectPow writes its principal root with square roots of real algebraic numbers *)
RadicalsTameQ[e_] := AllTrue[Cases[e, Power[b_, k_Rational] :> {b, k}, {0, Infinity}],
  Function[{bk}, With[{b = bk[[1]], k = bk[[2]]},
    MatchQ[b, _Integer | _Rational] ||
    (Denominator[k] == 2 && NumericQ[b] && (FreeQ[b, _Complex] || MatchQ[b, Complex[_Integer | _Rational, _Integer | _Rational]] || RationalModulusQ[b]))]]];

(* RationalModulusQ[b]: b = p + q I formally with p, q real algebraic and p^2 + q^2 rational *)
RationalModulusQ[b_] := Module[{ii, ri},
  ri = FormalReIm[b, ii];
  ri =!= None && FreeQ[ri, _Complex] && MatchQ[RootReduce[Expand[ri[[1]]^2 + ri[[2]]^2]], _Integer | _Rational]];

(* TameQ[e]: True when every constant of e can be worked with: its radicals are roots of
   rational numbers, square roots of real algebraic numbers, or square roots of Gaussian
   rationals (Sqrt[1 - 2 I]).  A square root of a non-real number with irrational parts (the
   roots of a generic quartic), a higher root of a non-rational base (a Cardano form), a Root
   object without a radical form or a trigonometric constant leaves the sum alone *)
TameQ[e0_] := Module[{e = RootFree[e0]},
  FreeQ[e, _Root | _Re | _Im | _Abs | _Arg | _Cos | _Sin | _ArcTan] && RadicalsTameQ[e]];

(* Denest[e]: a numeric radical expression denested by RootReduce when the result is again
   in radicals (Sqrt[3 + 2 Sqrt[2]] = 1 + Sqrt[2]); e itself otherwise *)
Denest[e_] := If[NumericQ[e] && ! FreeQ[e, Power[_, _Rational]],
  With[{r = RootReduce[e]}, If[FreeQ[r, _Root], r, e]], e];

(* RectConst[c]: a constant in rectangular form re + I im with real algebraic parts: every
   radical b^k rewritten from the inside out by RectPow, then Expand; an expression with
   symbols is returned as it is *)
RectConst[c0_] := Module[{c = c0},
  If[! FreeQ[c, _Root], c = RootFree[c]];               (* a Root object in radicals when it has them *)
  If[! NumericQ[c] || MatchQ[c, _Integer | _Rational], Return[c]];
  Expand[c /. Power[b_, k_Rational] :> RectPow[RectConst[b], k]]];

(* RectPow[b, k]: b^k in rectangular form for a base b already rectangular: a real base by
   its sign, (-b)^k (-1)^k when negative (Sqrt[1 - Sqrt[5]] = I Sqrt[Sqrt[5] - 1], the sign
   read from 30 digits); a non-real base with k = j/2 by the principal square root
   Sqrt[p + q I] = Sqrt[(r+p)/2] + Sign[q] I Sqrt[(r-p)/2], r = Abs[p + q I], raised to the
   odd power 2k.  ComplexExpand is applied to (-1)^k only, where the kernel is exact *)
RectPow[b_, k_] := Module[{ri, re, im, r, s, ii},
  If[FreeQ[b, _Complex],
    s = If[MatchQ[b, _Integer | _Rational], Sign[b], Sign[N[b, 30]]];
    Return[If[s < 0, Denest[Expand[-b]^k] ComplexExpand[(-1)^k], Denest[b^k]]]];
  If[Denominator[k] != 2, Return[b^k]];
  ri = FormalReIm[b, ii];
  If[ri === None, Return[b^k]];
  {re, im} = ri;
  If[IsZero[im], Return[RectPow[re, k]]];
  r = Denest[Sqrt[Expand[re^2 + im^2]]];
  s = Denest[Sqrt[Expand[(r + re)/2]]] + Sign[N[im, 30]] I Denest[Sqrt[Expand[(r - re)/2]]];
  If[k == 1/2, s, s^(2 k)]];

(* RectPoly[pol, vars]: the polynomial pol in vars with every constant coefficient in
   rectangular form *)
RectPoly[pol_, vars_] := Module[{rules = CoefficientRules[pol, vars]},
  If[! ListQ[rules], Return[pol]];
  FromCoefficientRules[MapAt[RectConst, rules, {All, 2}], vars]];

(* SplitOn[e, z, p]: {A, B} with e = A + z B formally, z a kernel with z^2 = p (the symbol
   standing for I with p = -1, or for Sqrt[p]): the polynomial in z reduced modulo z^2 - p;
   None when e is not polynomial in z *)
SplitOn[e_, z_, p_] := Module[{n, d, dz, r},
  {n, d} = NumDen[Expand[e]];
  If[! FreeQ[d, z],                                        (* rationalised by the conjugate d(-z) *)
    dz = d /. z -> -z;
    n = Expand[n dz]; d = PolynomialRemainder[Expand[d dz], z^2 - p, z];
    If[! FreeQ[d, z], Return[None]]];
  If[! PolynomialQ[Expand[n], z], Return[None]];
  r = PolynomialRemainder[Expand[n], z^2 - p, z];
  {CanL[Coefficient[r, z, 0]/d], CanL[Coefficient[r, z, 1]/d]}];

(* FormalReIm[e, ii]: {A, B} with e = A + I B formally, I replaced by the symbol ii *)
FormalReIm[e_, ii_] := If[! FreeQ[e, _Root | _Re | _Im | _Abs | _Arg | _Conjugate], None,
  SplitOn[e /. Complex[a_, b_] :> a + b ii, ii, -1]];

(* SameQ0[a, b, vars]: a - b = 0, tested structurally after Expand, then numerically at a fixed
   point of the variables (a clear nonzero settles it), and exactly only when the value is
   tiny -- IsZero's Simplify on two large Miller logands with Root-object coefficients is what
   the cheap tests avoid *)
SameQ0[a_, b_, vars_] := Module[{d = Expand[a - b], v, rules},
  If[d === 0, Return[True]];
  v = Quiet[N[d /. Thread[vars -> Table[13/10 + lrk/7, {lrk, Length[vars]}]], 30]];
  If[NumericQ[v] && Abs[v] > 10^-15, Return[False]];
  (* the exact test coefficient by coefficient: Can of this kernel leaves nested radicals in
     non-canonical forms, which IsZero's RootReduce settles one constant at a time *)
  rules = CoefficientRules[d, vars];
  If[ListQ[rules] && AllTrue[rules[[All, 2]], NumericQ], AllTrue[rules[[All, 2]], IsZero], IsZero[d]]];

(* MonicIn[u, vars]: u divided by the coefficient of its leading monomial in vars, canonical *)
MonicIn[u_, vars_] := CanL[u/First[CoefficientRules[u, vars]][[2]]];

(* TotalDegree[p, vars]: the total degree of a polynomial in vars *)
TotalDegree[p_, vars_] := Max[0, Total /@ Keys[CoefficientRules[p, vars]]];

(* KGcd[f, g]: the gcd of two polynomials over the number field of their constants *)
KGcd[f_, g_] := PolynomialGCD[f, g, Extension -> Automatic];

(* RiobooAtan[A, B, g]: LogToAtan of Bronstein (Rioboo): A, B polynomials in g over a field
   (the constants and the other generators), B != 0; a sum of arctangents of polynomials in
   g whose derivative is that of I Log[(A + I B)/(A - I B)] *)
RiobooAtan[A0_, B0_, g_] := Module[{A = CanL[A0], B = CanL[B0], gg, dd, cc},
  If[CanL[PolynomialRemainder[A, B, g]] === 0, Return[2 ArcTan[CanL[A/B]]]];
  If[Exponent[A, g] < Exponent[B, g], Return[RiobooAtan[-B, A, g]]];
  {gg, {dd, cc}} = PolynomialExtendedGCD[B, -A, g];      (* dd B - cc A = gg = gcd(A, B) *)
  {dd, cc, gg} = CanRaw /@ {dd, cc, gg};
  2 ArcTan[CanL[(A dd + B cc)/gg]] + RiobooAtan[dd, cc, g]];

(* SturmCount[P, g]: the number of distinct real roots of the polynomial P in g with real
   algebraic coefficients, by Sturm's theorem (the kernel has no CountRoots): the squarefree
   part, the remainder sequence with every coefficient canonical (RootReduce, so that a
   vanishing leading coefficient is seen as 0), the signs of the leading coefficients read
   from 30 digits when they are not rational.  None when a coefficient is not a real
   algebraic number *)
SturmCount[P0_, g_] := Module[{P, canon, seq, r, lc, sgn, changes},
  canon[p_] := Collect[Expand[p], g, RootReduce];
  P = canon[P0];
  If[! AllTrue[CoefficientList[P, g], NumericQ[#] && FreeQ[#, _Complex] &], Return[None]];
  If[Exponent[P, g] <= 0, Return[0]];
  P = canon[PolynomialQuotient[P, PolynomialGCD[P, D[P, g], Extension -> Automatic], g]];
  seq = {P, canon[D[P, g]]};
  While[Exponent[Last[seq], g] > 0,
    r = canon[-PolynomialRemainder[seq[[-2]], seq[[-1]], g]];
    If[r === 0, Break[]];
    AppendTo[seq, r]];
  lc[p_] := Coefficient[p, g, Exponent[p, g]];
  sgn[v_] := If[MatchQ[v, _Integer | _Rational], Sign[v], Sign[N[v, 30]]];
  changes[s_] := Module[{t = Select[s, # != 0 &], n = 0},
    Do[If[t[[lri]] != t[[lri + 1]], n++], {lri, Length[t] - 1}]; n];
  changes[sgn[lc[#] (-1)^Exponent[#, g]] & /@ seq] - changes[sgn[lc[#]] & /@ seq]];

(* RealZeroFreeQ[P, T, Y]: True when the polynomial P in the generators and Y is known to
   have no zero on the real locus of the curve: its norm, a polynomial in the curve variable
   alone with numeric coefficients, has no real root (SturmCount) apart from the roots of q
   for even m, which bound the real locus; None when undecided *)
RealZeroFreeQ[P_, T_, Y_] := Module[{gq, g, NN, h, m = T["m"], q = T["q"], cnt},
  gq = Select[T["gens"], SubsetQ[{#}, Variables[q]] &];
  If[Length[gq] != 1, Return[None]];
  g = First[gq];
  NN = Quiet[Check[CanL[TNorm[T, FromY[T, P, Y]]], $Failed]];
  If[NN === $Failed || ! SubsetQ[{g}, Variables[NN]] || ! AllTrue[CoefficientList[NN, g], NumericQ], Return[None]];
  If[EvenQ[m], While[True, h = KGcd[NN, q]; If[Exponent[h, g] <= 0, Break[]]; NN = CanL[NN/h]]];
  If[FreeQ[NN, g], Return[True]];
  cnt = SturmCount[NN, g];
  If[cnt === None, None, cnt == 0]];

(* YQuot[A, B, T, Y]: the quotient A/B of two polynomials in the generators and Y as a
   polynomial in Y over the generators (the division of the tower coordinates), so that the
   radical never sits in a denominator: the surface layer rewrites its positive powers only.
   For m = 2 the conjugate formula A B(-Y)/N(B) is written here with CanRaw (the package's
   Pdiv canonicalises through the field detour of Can, see CanL); TDiv for m >= 3 *)
YQuot[A_, B_, T_, Y_] := Module[{Bbar, rel = Y^2 - T["q"]},
  If[T["m"] != 2, Return[CanL[ToY[T, TDiv[T, FromY[T, A, Y], FromY[T, B, Y]], Y]]]];
  Bbar = B /. Y -> -Y;
  CanL[PolynomialRemainder[Expand[A Bbar], rel, Y]/PolynomialRemainder[Expand[B Bbar], rel, Y]]];

(* LogToAtan[A, B, T, Y]: LogToAtan for a pair of logands A +- I B (polynomials in the
   generators and Y): Rioboo's sum of arctangents of polynomials in the top generator when
   the radical is absent; over the radical, 2 ArcTan[A/B] or -2 ArcTan[B/A] with the
   denominator known free of real zeros when one is, else 2 ArcTan[A/B] *)
LogToAtan[A0_, B0_, T_, Y_] := Module[{A = CanL[A0], B = CanL[B0], tops},
  If[B === 0 || A === 0, Return[0]];
  If[T["q"] === None || (FreeQ[A, Y] && FreeQ[B, Y]),
    tops = Select[T["gens"], ! FreeQ[A, #] || ! FreeQ[B, #] &];
    If[tops === {}, Return[0]];                            (* the arctangent of a constant *)
    Return[RiobooAtan[A, B, Last[tops]]]];
  If[RealZeroFreeQ[B, T, Y] === True, Return[2 ArcTan[YQuot[A, B, T, Y]]]];
  If[RealZeroFreeQ[A, T, Y] === True, Return[-2 ArcTan[YQuot[B, A, T, Y]]]];
  2 ArcTan[YQuot[A, B, T, Y]]];

(* $lrZs: prime p -> the symbol standing for Sqrt[p] in the current call of RewriteLogs *)
$lrZs = <||>;
ZSym[p_] := If[KeyExistsQ[$lrZs, p], $lrZs[p], $lrZs[p] = Unique["r"]];

(* PrimeRadicals[e]: e with every power n^(k/2) of a positive rational written on the square
   roots of primes, each the symbol ZSym[p], so that Sqrt[6] and Sqrt[2] Sqrt[3] agree and
   Sqrt[p] -> -Sqrt[p] is a substitution *)
PrimeRadicals[e_] := e /. Power[b_?(MatchQ[#, _Integer | _Rational] && # > 0 &), k_Rational /; Denominator[k] == 2] :>
  With[{j = (Numerator[k] - 1)/2, fi = FactorInteger[Numerator[b] Denominator[b]]},
    b^j (1/Denominator[b]) Times @@ (#[[1]]^Quotient[#[[2]], 2] & /@ fi) Times @@ (If[OddQ[#[[2]]], ZSym[#[[1]]], 1] & /@ fi)];

(* LogToAtanh[A, B, z, p, T, Y]: the pair c Log[A + z B] + conj as 2 ArcTanh[A/(z B)] when B
   divides A (its argument a polynomial; z^2 = p), else None -- the logarithms are kept *)
LogToAtanh[A0_, B0_, z_, p_, T_, Y_] := Module[{A = CanL[A0], B = CanL[B0], ok, tops},
  If[B === 0 || A === 0, Return[0]];
  If[T["q"] =!= None && (! FreeQ[A, Y] || ! FreeQ[B, Y]),
    ok = FreeQ[B, Alternatives @@ Append[T["gens"], Y]],
    tops = Select[T["gens"], ! FreeQ[A, #] || ! FreeQ[B, #] &];
    ok = tops === {} || CanL[PolynomialRemainder[A, B, Last[tops]]] === 0];
  If[ok, 2 ArcTanh[CanL[A/(z B)]], None]];

(* RewriteLogs[terms, T, Y, real, verbose]: the logarithmic part Sum_k c_k Log[u_k] --
   terms = {{c_k, u_k}}, u_k expressions in the generators and the symbol Y of the radical
   -- with its conjugate pairs collapsed to real logarithms, arctangents and hyperbolic
   arctangents (LogToReal).  real: the integrand, derivations and radicand are free of I,
   so that a non-real term without a partner may be averaged with its formal conjugate.
   Returns {expression, pairs to arctan, pairs to artanh, averaged}, averaged when the sum was
   replaced by its average with the formal conjugate (the polynomial part must then be too). *)
RewriteLogs[terms_, T_, Y_, real0_, verbose_] := Module[
  {real = real0, vars, rel = None, red, work = {}, c, u, n, d, pol, hasI, i, j, cands, cj, uj, h, hit, merged, found,
   out = {}, reals = {}, natan = 0, natanh = 0, ii = Unique["i"], ric, riu, a, b, A, B, ubar, cbar,
   partner, ck, rest, half, kernels, conjs, done, sc, su, hh, L, back},
  vars = If[T["q"] === None, T["gens"], Append[T["gens"], Y]];
  If[T["q"] =!= None, rel = Y^T["m"] - T["q"]];
  If[! AllTrue[terms, TameQ[#[[1]]] && TameQ[#[[2]]] &], Return[{Total[#[[1]] Log[#[[2]]] & /@ terms], 0, 0, False}]];
  red[e_] := Module[{e2 = Expand[e]}, If[rel =!= None && ! FreeQ[e2, Y], e2 = PolynomialRemainder[e2, rel, Y]]; CanL[e2]];
  Do[c = lrcu[[1]]; u = lrcu[[2]];
    If[CanL[c] === 0, Continue[]];
    {n, d} = NumDen[CanL[u]];
    Do[pol = red[RectPoly[lrcp[[2]], vars]];
      If[! FreeQ[pol, Alternatives @@ vars], AppendTo[work, {RectConst[lrcp[[1]]], pol}]], {lrcp, {{c, n}, {-c, d}}}],
    {lrcu, terms}];
  hasI = ! FreeQ[work, _Complex];
  If[! hasI && FreeQ[work, Power[_?(MatchQ[#, _Integer | _Rational] && # > 0 &), k_Rational /; Denominator[k] == 2]],
    Return[{Total[#[[1]] Log[#[[2]]] & /@ work], 0, 0, False}]];
  If[hasI,
    (* logands sharing a factor with another logand or its conjugate are split *)
    i = 1;
    While[i <= Length[work],
      {c, u} = work[[i]];
      cands = If[FreeQ[u, _Complex], {u}, {u, CanL[u /. Complex[a1_, b1_] :> Complex[a1, -b1]]}];
      j = 1;
      While[j <= Length[work],
        If[j != i && FreeQ[work[[j, 2]], _Complex] && ! (TotalDegree[u, vars] <= 1 && TotalDegree[work[[j, 2]], vars] <= 1),
          {cj, uj} = work[[j]]; hit = False;        (* a product of conjugate factors is real; two linear logands share a factor only when equal: merged below *)
          Do[If[! hit,
              h = KGcd[uj, lrv];
              If[! FreeQ[h, Alternatives @@ vars] && ! SameQ0[MonicIn[h, vars], MonicIn[uj, vars], vars],
                work[[j]] = {cj, red[CanL[uj/h]]};
                AppendTo[work, {cj, red[h]}];
                hit = True]], {lrv, cands}]];
        j++];
      i++]];
  (* equal logands (up to a constant factor) merged *)
  merged = {};
  Do[u = MonicIn[lrcu[[2]], vars]; found = False;
    Do[If[! found && SameQ0[merged[[lrk, 2]], u, vars], merged[[lrk]] = {CanL[merged[[lrk, 1]] + lrcu[[1]]], merged[[lrk, 2]]}; found = True], {lrk, Length[merged]}];
    If[! found, AppendTo[merged, {lrcu[[1]], u}]], {lrcu, work}];
  work = Select[merged, ! IsZero[#[[1]]] &];
  (* the averaging with the formal conjugate is all or nothing: a term that does not split
     would keep an imaginary part that the others no longer compensate *)
  If[real && AnyTrue[work, FormalReIm[#[[1]], ii] === None || FormalReIm[#[[2]], ii] === None &], real = False];
  (* the conjugate pairs under I *)
  While[work =!= {},
    {c, u} = Last[work]; work = Most[work];
    ric = FormalReIm[c, ii]; riu = FormalReIm[u, ii];
    If[ric === None || riu === None, AppendTo[out, c Log[u]]; Continue[]];
    {a, b} = ric; {A, B} = riu;
    If[IsZero[b] && IsZero[B], AppendTo[reals, {c, u}]; Continue[]];
    ubar = MonicIn[A - I B, vars];
    If[SameQ0[ubar, u, vars],                              (* a real logand with a non-real coefficient *)
      If[real, AppendTo[reals, {a, u}], AppendTo[out, c Log[u]]]; Continue[]];
    cbar = CanL[a - I b];
    partner = 0; Do[If[partner == 0 && SameQ0[work[[lrk, 2]], ubar, vars], partner = lrk], {lrk, Length[work]}];
    half = None;
    If[partner > 0,
      ck = work[[partner, 1]]; work = Delete[work, partner];
      rest = CanL[ck - cbar];
      If[! IsZero[rest], AppendTo[work, {rest, ubar}]];
      half = 1,
      If[real, half = 1/2]];                               (* (c Log[u] + conj) / 2 *)
    If[half === None, AppendTo[out, c Log[u]]; Continue[]];
    If[! IsZero[a], AppendTo[reals, {half a, red[A^2 + B^2]}]];
    AppendTo[out, half b LogToAtan[A, B, T, Y]];
    natan++];
  (* the pairs under Sqrt[p] -> -Sqrt[p]: a kernel whose conjugation is a substitution (no
     nested radical contains it) and whose pair has a polynomial ArcTanh argument *)
  $lrZs = <||>;
  reals = {PrimeRadicals[#[[1]]], PrimeRadicals[#[[2]]]} & /@ reals;
  kernels = Cases[reals, Power[_, _Rational], {0, Infinity}];
  conjs = Select[Normal[$lrZs], With[{z = #[[2]]}, AllTrue[kernels, FreeQ[#, z] &]] &];
  While[reals =!= {},
    {c, u} = Last[reals]; reals = Most[reals];
    done = False;
    Do[If[! done && (! FreeQ[c, lrpz[[2]]] || ! FreeQ[u, lrpz[[2]]]),
        sc = SplitOn[c, lrpz[[2]], lrpz[[1]]]; su = SplitOn[u, lrpz[[2]], lrpz[[1]]];
        If[sc =!= None && su =!= None && ! IsZero[su[[2]]],
          {a, b} = sc; {A, B} = su;
          ubar = MonicIn[A - lrpz[[2]] B, vars];
          partner = 0; Do[If[partner == 0 && SameQ0[reals[[lrk, 2]], ubar, vars], partner = lrk], {lrk, Length[reals]}];
          If[partner > 0,
            hh = LogToAtanh[A, B, lrpz[[2]], lrpz[[1]], T, Y];
            If[hh =!= None,
              ck = reals[[partner, 1]]; reals = Delete[reals, partner];
              rest = CanL[ck - (a - lrpz[[2]] b)];
              If[! IsZero[rest], AppendTo[reals, {rest, ubar}]];
              AppendTo[out, a Log[red[A^2 - lrpz[[1]] B^2]] + b lrpz[[2]] hh];
              natanh++; done = True]]]], {lrpz, conjs}];
    If[! done, AppendTo[out, c Log[u]]]];
  back = (#[[2]] -> Sqrt[#[[1]]]) & /@ Normal[$lrZs];
  L = Total[out] /. back;
  L = L /. (fn : (ArcTan | ArcTanh | Log))[arg_] /; ! FreeQ[arg, Alternatives @@ vars] :> fn[CanL[Expand[arg]]];
  If[verbose && natan + natanh > 0,
    Print["  real form: ", natan + natanh, " pair(s) of logarithms collapsed (Rioboo): ", natan, " to arctan, ", natanh, " to artanh"]];
  {L, natan, natanh, real && hasI}];

(* LogToReal[terms, rat, f, T, Y, verbose]: the hook: terms = {{c_k, u_k}} the logarithmic
   part in the tower coordinates (u_k an expression in the generators and the symbol Y of
   the radical), rat the list of the polynomial part, f the integrand list.  Returns {the
   logarithmic part rewritten, rat}, rat with its imaginary part dropped when the
   logarithmic part was averaged with its formal conjugate. *)
$RewriteBudget = 30;         (* seconds the rewrite may take before the sum is left as it is *)
LogToReal[terms_, rat_, f_, T_, Y_, verbose_] := Module[{real, r, rat2 = rat},
  real = FreeQ[{f, T["derivs"], T["q"]}, _Complex];
  r = TimeConstrained[RewriteLogs[terms, T, Y, real, verbose], $RewriteBudget, $Aborted];
  If[r === $Aborted,
    If[verbose, Print["  real form: not attempted within ", $RewriteBudget, " s, the logarithms are kept"]];
    Return[{Total[#[[1]] Log[#[[2]]] & /@ terms], rat}]];
  If[r[[4]], rat2 = RealPart[#, T["gens"]] & /@ rat];
  {r[[1]], rat2}];

(* RealPart[r, gens]: the formal real part of a rational function of the generators with
   algebraic constants (its constants in rectangular form first: Sqrt[1 - Sqrt[5]] hides an I) *)
RealPart[r_, gens_] := Module[{ii = Unique["i"], n, d, dc, ri},
  {n, d} = NumDen[r]; n = RectPoly[n, gens]; d = RectPoly[d, gens];
  If[! FreeQ[d, _Complex], dc = d /. Complex[a_, b_] :> Complex[a, -b]; n = Expand[n dc]; d = Expand[d dc]];
  If[FreeQ[n, _Complex], Return[If[FreeQ[d, _Complex], r, CanL[n/d]]]];
  ri = FormalReIm[n, ii];
  If[ri === None, r, CanL[ri[[1]]/d]]];
