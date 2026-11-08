(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 236 -- SIAM 100-Digit Challenge, sealed holdout B for Appendix D
   problem 1 (Kempner-type series): ten other digit-restricted series.

     bd1.1   Sum 1/n,   n with no decimal digit 9 (Kempner)       = 22.92067661...
     bd1.2   Sum 1/n,   n with no decimal digit 0                 = 23.10344790...
     bd1.3   Sum 1/n,   n whose decimal digits are all odd        = 3.171765473...
     bd1.4   Sum 1/n,   n without the decimal substring "314"     = 2299.829782...
     bd1.5   Sum 1/n,   n with exactly one decimal digit 9 (Irwin) = 23.04428708...
     bd1.6   Sum 1/n^2, n with no decimal digit 7                 = 1.614031813...
     bd1.7   Sum 1/n,   n with no "11" in binary (fibbinary)      = 3.704711752...
     bd1.8   Sum 1/n,   n with no digit 2 in base 3               = 2.682853110...
     bd1.9   Sum 2^-n,  n with no "11" in binary                  = 0.850614071...
     bd1.10  Sum 1/n,   n with no two equal adjacent decimal digits = 23.38295748...

   Idiomatic statement only, as in experiment 176: NSum to Infinity of the
   indicator-weighted term (a ?NumericQ helper, since NSum samples the summand
   at numeric n), at machine precision and at WorkingPrecision 60.  The
   true values come from the automaton power-sum recursion (Schmelzer-Baillie,
   generalised); no system has that built in, so most rows are expected to
   fail -- except bd1.6 (absolutely convergent like n^-2) and bd1.9
   (geometric), which a careful summer can do, and bd1.7/bd1.8, whose terms
   thin out geometrically (ratio phi/2 and 2/3 per digit).
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSum", "IntegerString", "IntegerDigits", "StringFreeQ", "Differences"}];

ok[1][m_] := StringFreeQ[IntegerString[m], "9"];
ok[2][m_] := StringFreeQ[IntegerString[m], "0"];
ok[3][m_] := StringFreeQ[IntegerString[m], {"0", "2", "4", "6", "8"}];
ok[4][m_] := StringFreeQ[IntegerString[m], "314"];
ok[5][m_] := Count[IntegerDigits[m], 9] == 1;
ok[6][m_] := StringFreeQ[IntegerString[m], "7"];
ok[7][m_] := StringFreeQ[IntegerString[m, 2], "11"];
ok[8][m_] := StringFreeQ[IntegerString[m, 3], "2"];
ok[9][m_] := StringFreeQ[IntegerString[m, 2], "11"];
ok[10][m_] := FreeQ[Differences[IntegerDigits[m]], 0];
wt[6][n_] := 1/n^2;
wt[9][n_] := 2^-n;
wt[_][n_] := 1/n;
term[p_][n_?NumericQ] := With[{m = Round[n]}, If[TrueQ[ok[p][m]], wt[p][n], 0]];
(* NHoldAll: NSum applies N to the summand, which would turn term[9][n] into
   term[9.][n] and select no ok/wt definition at all. *)
SetAttributes[term, NHoldAll];

Do[
  vrun["bd1." <> ToString[p] <> " machine", NSum[term[p][n], {n, 1, Infinity}], 10, 1, 60];
  vrun["bd1." <> ToString[p] <> " wp50",
       NSum[term[p][n], {n, 1, Infinity}, WorkingPrecision -> 60], 48, 1, 60],
  {p, 1, 10}];
