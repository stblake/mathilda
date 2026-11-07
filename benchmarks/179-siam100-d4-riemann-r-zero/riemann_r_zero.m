(* bench-timeout: 400 *)
(* ==========================================================================
   Experiment 179 -- SIAM 100-Digit Challenge, Appendix D Problem 4
   (Joerg Waldvogel).

       R(x) = Sum_{k>=1} mu(k)/k li(x^(1/k))   (Riemann's prime-counting
                                                 function, li a principal value)
       What is the largest positive zero of R?

       x* = 1.8286432697...e-14828   (t* = log x* = -34142.1281846...)

   "The answer is truly shocking": R > 0 for all x > x*; near 0, R(e^-u) =
   21.9/u^3 plus an oscillation of amplitude about 1e-10 u^(-1/2) driven by the
   first zeta zero, and the two first balance at u = 34142.  x* is far below
   the smallest double (it underflows), and R at x* is a cancellation of terms
   of size 1 down to about 1e-13, so the Moebius/li series as written cannot
   locate it in any precision a user would choose.

   Rows (t = log x is the variable, started at t = -34000, i.e. the user is
   told the magnitude; the answer recorded is x = Exp[t*]):
     d4.machine        FindRoot on the series as written, the infinite sum
                       by NSum, machine precision;
     d4.wp50           the same with WorkingPrecision -> 60 (NSum and
                       FindRoot);
     d4.wp50-builtin   FindRoot on the BUILT-IN RiemannR (WL has it, mpmath
                       has riemannr; Maxima and FriCAS do not), WP 60.
   Expected: the two series rows return confident wrong answers in every
   system (the measurement is whether a system SAYS so); only a built-in R
   with its own high-precision evaluation can get the digits.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindRoot", "NSum", "MoebiusMu", "LogIntegral", "RiemannR"}];

rr[t_?NumericQ] := NSum[MoebiusMu[k]/k LogIntegral[Exp[t/k]], {k, 1, Infinity}];
rr60[t_?NumericQ] := NSum[MoebiusMu[k]/k LogIntegral[Exp[t/k]], {k, 1, Infinity},
                          WorkingPrecision -> 60];

vrun["d4.machine", Exp[ruleVal[FindRoot[rr[t] == 0, {t, -34000}], t]], 10];
vrun["d4.wp50", Exp[ruleVal[FindRoot[rr60[t] == 0, {t, -34000},
                                     WorkingPrecision -> 60], t]], 48];
vrun["d4.wp50-builtin", Exp[ruleVal[FindRoot[RiemannR[Exp[t]] == 0, {t, -34000},
                                             WorkingPrecision -> 60], t]], 48];
