(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 176 -- SIAM 100-Digit Challenge, Appendix D problem 1
   (a Kempner-type series).

       Sum of 1/n over the positive integers n whose decimal representation
       does not contain the substring "42".               = 228.4463041592...

   Idiomatic statement only: NSum of the indicator-weighted harmonic term to
   Infinity, at machine precision and at WorkingPrecision 60.  The summand is
   a ?NumericQ helper (IntegerString needs an integer, and NSum samples the
   summand at numeric n).  No built-in in any system knows the
   Schmelzer-Baillie / Burnol digit-recursion that actually sums this series;
   the series converges like n^-0.0044 in its tail (the admissible n thin out
   only by a factor 9.899/10 per decade), so a generic extrapolating NSum is
   expected to FAIL everywhere -- the case records how each system fails.
   True value and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSum", "IntegerString", "StringFreeQ"}];

k42[n_?NumericQ] := Boole[StringFreeQ[IntegerString[Round[n]], "42"]]/n;

vrun["d1.machine", NSum[k42[n], {n, 1, Infinity}], 10];
vrun["d1.wp50", NSum[k42[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48];
