(* ::Package:: *)

(* ========================================================================== *)
(* LogToArcTan.m -- presentation of the logarithmic part: tidy logarithms,    *)
(*                  expand small RootSums, and rewrite complex logarithms in   *)
(*                  real form with ArcTan.                                    *)
(*                                                                            *)
(* Real form.  Treat x and every radical as real symbols and write each       *)
(* coefficient and logand as c = c1 + I c2, u = A + I B (A, B, c1, c2 free of  *)
(* I).  Formally                                                              *)
(*     Re(c Log[u]) = c1/2 Log[A^2 + B^2] - c2 ArcTan[B/A],                   *)
(* and d/dx of Log[A^2+B^2]/2 and ArcTan[B/A] are the formal real and         *)
(* imaginary parts of u'/u.  Because the integrand is free of I, the formal   *)
(* real part of an antiderivative is again an antiderivative, so the rewrite  *)
(* is applied to the whole answer at once; the caller re-checks the result    *)
(* and keeps the complex form if the check fails.  (This replaces the         *)
(* package's paired-term rules A/B, which were applied one after the other    *)
(* rather than jointly and missed unpaired terms.)                            *)
(* ========================================================================== *)

togetherAtoms::usage = "togetherAtoms[e] is Together[e] with every radical (Power[_, _Rational]) first replaced by a fresh symbol and restored afterwards: the same result for expressions in which the radicals are independent atoms, computed by polynomial arithmetic.";
togetherAtoms[e_] := Module[{atoms, rules},
  atoms = DeleteDuplicates[Cases[e, Power[_, _Rational], {0, Infinity}]];
  rules = Thread[atoms -> Table[Unique["ra"], {Length[atoms]}]];
  Together[e /. rules] /. (Reverse /@ rules)];

tidyLog::usage = "tidyLog[u, x] simplifies a logand u: Together, then the numeric content of numerator and denominator is dropped (Log[c u] - Log[u] is a constant).";
tidyLog[u_, x_Symbol] := Module[{t = togetherAtoms[u], num, den},
  num = Numerator[t]; den = Denominator[t];
  num = dropNumericContent[num]; den = dropNumericContent[den];
  num/den];

dropNumericContent::usage = "dropNumericContent[p] divides the expression p by its numeric content (FactorTermsList), leaving a primitive expression; p is returned unchanged when the content is not a nonzero rational number.";
dropNumericContent[p_] := Module[{ftl = Quiet[FactorTermsList[Expand[p]]], c},
  c = If[ListQ[ftl] && Length[ftl] >= 1, First[ftl], 1];
  If[NumberQ[c] && c =!= 0 && FreeQ[c, Complex], Expand[p/c], p]];

expandRootSums::usage = "expandRootSums[e] replaces every RootSum[q, body] whose polynomial q has degree <= 2 by the explicit sum of body over the roots of q.";
expandRootSums[e_] := e /. RootSum[q_, body_] :> Module[{s = Unique["s"], rts},
  If[Exponent[q[s], s] <= 2,
    rts = s /. Solve[q[s] == 0, s];
    Total[body /@ rts],
    RootSum[q, body]]];

complexSplit::usage = "complexSplit[e] is {re, im} with e = re + I im formally, treating every symbol and radical as real: I is replaced by a symbol ii, the numerator is reduced modulo ii^2 + 1 and the denominator rationalised by its conjugate.";
complexSplit[e_] := Module[{ii = Unique["ii"], t, num, den, red, nr, ni, dr, di, d2, atoms, rules, back},
  red[p_] := Expand[PolynomialRemainder[Expand[p], ii^2 + 1, ii]];
  (* Every radical (of x or of a number) is a real atom here: replace them by *)
  (* symbols so that Together is plain polynomial arithmetic, then restore.  *)
  atoms = DeleteDuplicates[Cases[e, Power[_, _Rational], {0, Infinity}]];
  rules = Thread[atoms -> Table[Unique["ra"], {Length[atoms]}]];
  back = Reverse /@ rules;
  t = Together[(e /. rules) /. Complex[a_, b_] :> a + b ii];
  num = red[Numerator[t]]; den = red[Denominator[t]];
  {nr, ni} = {Coefficient[num, ii, 0], Coefficient[num, ii, 1]};
  {dr, di} = {Coefficient[den, ii, 0], Coefficient[den, ii, 1]};
  d2 = Expand[dr^2 + di^2];
  {Together[(nr dr + ni di)/d2] /. back, Together[(ni dr - nr di)/d2] /. back}];

realFormTerm::usage = "realFormTerm[c, u, x] is the formal real part of c Log[u] (see the file header): c1/2 Log[A^2 + B^2] - c2 ArcTan[B/A] with c = c1 + I c2, u = A + I B; terms free of I are returned as c Log[tidyLog[u]].";
realFormTerm[c_, u_, x_Symbol] := Module[{c1, c2, A, B},
  If[FreeQ[{c, u}, Complex], Return[c Log[tidyLog[u, x]]]];
  {c1, c2} = complexSplit[c];
  {A, B} = complexSplit[u];
  Which[
    Together[B] === 0, c1 Log[tidyLog[A, x]],
    Together[A] === 0, c1 Log[tidyLog[B, x]],
    True, c1/2 Log[tidyLog[A^2 + B^2, x]] + arcTanTerm[c2, A, B, x]]];

arcTanTerm::usage = "arcTanTerm[c2, A, B, x] is -c2 ArcTan[B/A] up to an additive constant: either -c2 ArcTan[B/A] or c2 ArcTan[A/B], whichever quotient's denominator (in lowest terms) changes sign fewer times on a sample of the real line. The two differ by the piecewise constant c2 Pi/2 Sign[A B], so either is an antiderivative; the choice only keeps a definite integral's Newton-Leibniz evaluation away from spurious jumps.";
arcTanTerm[c2_, A_, B_, x_Symbol] := Module[{ab = togetherAtoms[A/B], ba = togetherAtoms[B/A]},
  (* ArcTan[f] jumps where f passes through infinity: where the denominator  *)
  (* of f (in lowest terms) changes sign, not where A or B do.               *)
  If[signChanges[Denominator[ab], x] < signChanges[Denominator[ba], x],
    c2 ArcTan[ab],
    -c2 ArcTan[ba]]];

signChanges::usage = "signChanges[e, x] counts the sign changes of the real-valued expression e at 15 sample points spread over [-10, 10] (denser near 0; points where e is not a real number are skipped). A heuristic for the ArcTan orientation only -- the answer is re-verified either way.";
signChanges[e_, x_Symbol] := Module[{vals, pts = {-10, -5, -3, -2, -1, -1/2, -1/5, 1/97, 1/5, 1/2, 1, 2, 3, 5, 10} + 1/193},
  (* machine numbers: substituting exact rationals into exact radicals is   *)
  (* two orders of magnitude slower, and a sign needs no precision          *)
  vals = Select[Quiet[N[e /. x -> N[#]]] & /@ pts, NumberQ[#] && Im[#] == 0 && # != 0 &];
  Count[Partition[Sign[Re[vals]], 2, 1], {a_, b_} /; a =!= b]];

realForm::usage = "realForm[F, x] rewrites the antiderivative F (a sum of algebraic terms, coef Log[u] terms and RootSums) in real form: small RootSums are expanded, every coef Log[u] term with complex numbers is replaced by its formal real part, and the remaining terms are kept. It returns F unchanged when it has no complex numbers.";
realForm[F_, x_Symbol] := Module[{G = expandRootSums[F], terms},
  If[FreeQ[G, Complex] && FreeQ[G, RootSum], Return[G /. Log[u_] :> Log[tidyLog[u, x]]]];
  terms = If[Head[Expand[G]] === Plus, List @@ Expand[G], {Expand[G]}];
  Total[Function[t,
    Which[
      MatchQ[t, Log[_]], realFormTerm[1, First[t], x],
      MatchQ[t, c_ Log[_] /; FreeQ[c, Log]],
        With[{lg = First[Cases[t, Log[_], {0, 1}]]}, realFormTerm[Together[t/lg], First[lg], x]],
      FreeQ[t, Complex], t,
      True, With[{cs = complexSplit[t]}, cs[[1]]]]] /@ terms]];
