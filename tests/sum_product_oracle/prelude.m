(* tests/sum_product_oracle/prelude.m

   Verifier for the Sum / Product oracle corpus (corpus.m).

   Every record is {label, Hold[input], kind, data, subs}.  spoCheck returns
   one verdict string per record:

     "PASS"    the answer is verified
     "WRONG"   a value came back and it is wrong (or a divergent / pole input
               came back as a finite value)
     "UNEVAL"  the input came back unevaluated where a value is expected
     "NOCHECK" the answer could not be numericalised for the comparison

   The check never trusts a closed form against another closed form:

     kind "value"     data = the exact (or 40-digit) reference value, written
                      by hand and checked against an mpmath partial sum
                      (tools in the corpus header).  The answer and the
                      reference are both taken /. subs and compared at 40
                      digits.  subs is applied AFTER evaluation, so a
                      symbolic parameter (Zeta[s], 1/(1-x)) is genuinely
                      exercised.
     "finite"         data = {n, nlo, nhi} or {n, nlo, nhi, off} (sample at
                      n = v + off).  The answer, a function of the
                      symbolic bound n, is compared at n = nlo..nhi with the
                      SAME input enumerated at that concrete integer bound
                      (a finite numeric range enumerates term by term, so
                      the reference is independent of every closed-form
                      stage).  A pole of the enumeration (ComplexInfinity)
                      must be matched by a non-finite answer.
     "indef"          data = {i, ilo, ihi}.  An indefinite Sum[f, i] /
                      Product[f, i] is an antidifference, defined up to a
                      constant: F(i+1) - F(i) == f(i) (ratio for Product) at
                      i = ilo..ihi.
     "pole"           a pole of the summand lies inside the range: the
                      answer must be ComplexInfinity / non-finite, or the
                      input held.  A finite value is WRONG.
     "divergent"      the series diverges: held or an infinity is PASS, any
                      other value (finite number, Indeterminate, a
                      symbolic expression) is WRONG.
     "zeroOrHeld"     a product that diverges to 0: 0 or held.
     "held"           must stay unevaluated (an undefined summand, a
                      malformed iterator).
*)

spoHeadQ[r_] := !FreeQ[r, Sum | Product, {0, Infinity}, Heads -> True];

spoNonFinite[x_] := !FreeQ[x, Indeterminate | ComplexInfinity | DirectedInfinity | Infinity];

(* 40-digit numeric value of an exact or inexact answer, or $Failed. *)
spoNum[x_] := Module[{v = Quiet[N[x, 40]]},
  If[NumberQ[v], v, $Failed]];

(* True / False / $Failed: does the answer a match the reference b? *)
(* No Return[] in these helpers: inside If within Module it leaves only the
   If, not the Module (see the .m pitfalls note), so every early exit is a
   Which branch. *)
spoClose[a_, b_] := Module[{na = spoNum[a], nb = spoNum[b]},
  If[na === $Failed || nb === $Failed, $Failed,
     TrueQ[Abs[na - nb] <= 10^-20 (1 + Abs[nb])]]];

(* Evaluate the held input quietly. *)
spoEval[h_Hold] := Quiet[ReleaseHold[h]];

spoVerdict[c_] := Which[c === $Failed, "NOCHECK", c, "PASS", True, "WRONG"];

spoCheckValue[r_, data_, subs_] := Module[{rs = r /. subs},
  Which[spoHeadQ[r], "UNEVAL",
        spoNonFinite[rs], "WRONG",
        True, spoVerdict[spoClose[rs, data /. subs]]]];

(* One sample point of a "finite" / "indef" record: answer value av against
   the enumerated reference ev. *)
spoPoint[av_, ev_] := Which[
  spoNonFinite[ev], If[spoNonFinite[av] || spoHeadQ[av], "PASS", "WRONG"],
  spoNonFinite[av], "WRONG",
  spoHeadQ[av], "UNEVAL",
  True, spoVerdict[spoClose[av, ev]]];

spoWorst[vs_List] := Which[
  MemberQ[vs, "WRONG"], "WRONG", MemberQ[vs, "UNEVAL"], "UNEVAL",
  MemberQ[vs, "NOCHECK"], "NOCHECK", True, "PASS"];

(* {n, nlo, nhi, off}: sample at n = v + off -- a range {k, 1/2, n} has its
   upper bound on the lattice 1/2 + Integers (the convention is that
   n - imin is an integer, as in the reference). *)
spoCheckFinite[h_Hold, r_, {n_, nlo_, nhi_, off_}, subs_] :=
  If[spoHeadQ[r], "UNEVAL",
    spoWorst[Table[
      spoPoint[Quiet[r /. subs /. n -> v + off],
               Quiet[ReleaseHold[h /. n -> v + off /. subs]]],
      {v, nlo, nhi}]]];

spoCheckFinite[h_Hold, r_, {n_, nlo_, nhi_}, subs_] :=
  If[spoHeadQ[r], "UNEVAL",
    spoWorst[Table[
      spoPoint[Quiet[r /. subs /. n -> v],
               Quiet[ReleaseHold[h /. n -> v /. subs]]],
      {v, nlo, nhi}]]];

(* Indefinite Sum[f, i] / Product[f, i]: an antidifference is defined only up
   to an additive (Sum) or multiplicative (Product) constant, so the check is
   F(v + 1) - F(v) == f(v) (resp. F(v + 1) / F(v) == f(v)) at v = ilo..ihi. *)
spoCheckIndef[h_Hold, r_, {i_, ilo_, ihi_}, subs_] := Module[{isSum, body},
  isSum = MatchQ[h, Hold[Sum[_, _]]];
  body = h /. Hold[(Sum | Product)[f_, _]] :> Hold[f];
  If[spoHeadQ[r], "UNEVAL",
    spoWorst[Table[
      spoPoint[
        Quiet[If[isSum, (r /. i -> v + 1) - (r /. i -> v),
                        (r /. i -> v + 1) / (r /. i -> v)] /. subs],
        Quiet[ReleaseHold[body /. i -> v /. subs]]],
      {v, ilo, ihi}]]]];

spoCheck[label_, h_Hold, kind_String, data_, subs_] := Module[{r = spoEval[h]},
  Switch[kind,
    "value",      spoCheckValue[r, data, subs],
    "finite",     spoCheckFinite[h, r, data, subs],
    "indef",      spoCheckIndef[h, r, data, subs],
    "pole",       If[spoHeadQ[r] || spoNonFinite[r /. subs], "PASS", "WRONG"],
    "divergent",  If[spoHeadQ[r] || MatchQ[r, Infinity | -Infinity | _DirectedInfinity | ComplexInfinity],
                     "PASS", "WRONG"],
    "zeroOrHeld", If[spoHeadQ[r] || r === 0, "PASS", "WRONG"],
    "held",       If[spoHeadQ[r], "PASS", "WRONG"],
    _,            "BADKIND"]];

(* Human-readable one-liner: the answer itself, for the driver's report. *)
spoShow[h_Hold] := ToString[InputForm[spoEval[h]]];
