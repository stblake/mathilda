(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 181 -- SIAM 100-Digit Challenge, Appendix D Problem 6
   ("quadruple exponential", Trefethen).

       Int_{-1}^{1}  exp(x + sin(exp(exp(x + 1/3))))  dx  =  3.0709937884...

   Reading of the formula: the book (p. 288, PDF p. 298) typesets the
   argument as  x + sin e^{e^{x+1/3}},  a tower of TWO exponentials under the
   sine, all inside the outer exp -- four transcendental levels in all, hence
   the name.  The inner phase exp(exp(x+1/3)) runs from 1.67 to 44.3 over
   [-1, 1], so the integrand makes about 14 sign changes of its sine, crowded
   toward x = 1; it is smooth (entire) and bounded by e^2.

   Idiomatic statement only: the integral as written, at machine precision
   and at 50 digits (WorkingPrecision 60, PrecisionGoal 50).
   True value and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate"}];

vrun["d6.machine", NIntegrate[Exp[x + Sin[Exp[Exp[x + 1/3]]]], {x, -1, 1}], 10, 3];
vrun["d6.wp50", NIntegrate[Exp[x + Sin[Exp[Exp[x + 1/3]]]], {x, -1, 1},
                           WorkingPrecision -> 60, PrecisionGoal -> 50], 48];
