(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 170 -- SIAM 100-Digit Challenge, Problem 5 ("Hard Work").

       f(z) = 1/Gamma(z); p(z) is the cubic that best approximates f on the
       unit disk in the sup norm.  What is ||f - p||_oo ?
                                                       =  0.2143352345...

   Idiomatic statement, and the number that is scored.  By the maximum
   principle the sup over the disk is the max over the circle |z| = 1; the
   direct formulation a user writes is NMinimize over the coefficients of
   the max over SAMPLED circle points:

       E_M = min over complex c0..c3 of
             max_k |1/Gamma(z_k) - (c0 + c1 z_k + c2 z_k^2 + c3 z_k^3)|,
       z_k = exp(2 pi i k/M),  k = 0..M-1,  M = 128,

   with the 8 real and imaginary parts of c0..c3 as the unknowns (no symmetry
   reduction).  Both rows are scored against E_128 itself, a well-defined
   number known to 70 digits (truth.tsv):
       E_128 = 0.21432875514988691993...
   which differs from the continuous challenge answer 0.21433523459045963946
   (book App. B) by 3.0e-5 relatively: the sampled max under-estimates the
   true max between samples.  (The continuous problem is a semi-infinite
   minimax problem that no system's built-ins state directly; the book solves
   it by a 6-equation equioscillation system -- a specialised method.)
   The objective is NON-SMOOTH (a max of 128 moduli), which is what the case
   tests: a derivative-free global minimiser has to resolve a kink at the
   optimum, where 9 of the 128 moduli are equal.
     c5.machine   NMinimize at machine precision;
     c5.wp50      NMinimize with WorkingPrecision -> 60.
   mpmath has no minimiser (wp50 SKIP); Maxima: fmin_cobyla on the same max
   (machine), no bigfloat minimiser (wp50 SKIP); FriCAS: no minimiser (SKIP).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NMinimize", "Gamma"}];

c5obj[prec_] := Module[{z = Exp[2 Pi I Range[0, 127]/128], zk, fk},
  zk = N[z, prec]; fk = N[1/Gamma[z], prec];      (* 1/Gamma[-1] = 0 exactly *)
  Max[Abs[fk - ((a0 + I b0) + (a1 + I b1) zk + (a2 + I b2) zk^2 + (a3 + I b3) zk^3)]]];

vrun["c5.machine", First[NMinimize[c5obj[MachinePrecision],
                                   {a0, a1, a2, a3, b0, b1, b2, b3}]], 10];
vrun["c5.wp50", First[NMinimize[c5obj[60], {a0, a1, a2, a3, b0, b1, b2, b3},
                                WorkingPrecision -> 60]], 48];
