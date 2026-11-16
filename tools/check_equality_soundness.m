(* check_equality_soundness.m -- differential soundness gate for the
   comparison heads, run as  ./Mathilda -file tools/check_equality_soundness.m
   (make check-equality).

   Every probe is a pair (c, r) of EXACT quantities whose order is known
   independently: r runs over the continued-fraction convergents of a constant
   c (ever closer, alternating sides), and over c shifted by 10^-k up to
   k = 400.  The sign of c - r is read from a 300-digit N (or is the sign of
   the shift), and every head must agree with it:

     ==  !=  <  <=  >  >=  PossibleZeroQ  Max  Sign[c - r]  I c == I r

   plus the invariants Equal xor Unequal and Less xor GreaterEqual.  Before
   v0.6584 263 of the 550 convergent probes were wrong (22/7 == Pi was True).
   Exits 1 on any failure, listing the first 40. *)

consts = {Pi, E, Sqrt[2], GoldenRatio, Log[2], EulerGamma, Catalan, Zeta[3],
  Pi^2, Sqrt[3] + Pi, Exp[Pi], Gamma[1/3], BesselJ[0, 1], Sin[1], ArcTan[2],
  Log[3]/Log[2], 2^(1/3), Khinchin, Glaisher, PolyLog[2, 1/3], Erf[1], Zeta[5]};

bad = 0; total = 0;

probe[c_, r_, s_] := Module[{got, want},
  total++;
  got = {c == r, c != r, c < r, c <= r, c > r, c >= r, PossibleZeroQ[c - r],
         Max[c, r] === If[s > 0, c, r], Sign[c - r] === s, I c == I r,
         Xor[TrueQ[c == r], TrueQ[c != r]], Xor[TrueQ[c < r], TrueQ[c >= r]]};
  want = {False, True, s < 0, s < 0, s > 0, s > 0, False,
          True, True, False, True, True};
  If[got =!= want,
    bad++;
    If[bad <= 40, Print["FAIL ", InputForm[c], " vs ", InputForm[r],
      "  sign ", s, "\n  got  ", got, "\n  want ", want]]]];

Do[
  With[{cf = ContinuedFraction[N[c, 400], 25]},
    Do[With[{r = FromContinuedFraction[Take[cf, k]]},
         With[{s = Sign[N[c - r, 300]]}, If[s =!= 0, probe[c, r, s]]]],
      {k, 1, Length[cf]}]];
  Do[probe[c + 10^-k, c, 1]; probe[c - 10^-k, c, -1], {k, {5, 15, 30, 60, 120, 400}}],
  {c, consts}];

Print["check-equality: ", total - bad, " / ", total, " probes agree"];
If[bad > 0, Exit[1]];
