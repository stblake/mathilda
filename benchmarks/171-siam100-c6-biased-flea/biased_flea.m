(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 171 -- SIAM 100-Digit Challenge, Problem 6 ("Biased Random Walk").

       A flea starts at (0,0) on the 2-D integer lattice and steps N or S with
       probability 1/4 each, E with 1/4 + eps, W with 1/4 - eps.  For which
       eps is the probability that it ever returns to (0,0) exactly 1/2?
                                                       =  0.0619139544...

   Idiomatic statement: the return probability from the lattice Green's
   function, p(eps) = 1 - 1/G(eps), with G the expected number of visits to
   the origin,
       G(eps) = (2 pi)^-2 Int_{[-pi,pi]^2} Re 1/(1 - phi(k)) dk
              = pi^-2 Int_0^pi Int_0^pi (1-c)/((1-c)^2 + s^2) dk1 dk2,
       c = (cos k1 + cos k2)/2,  s = 2 eps sin k1,
   (phi the characteristic function of one step; the integrand is even in
   both variables, so the user folds it onto [0,pi]^2 and takes the real part
   explicitly), and FindRoot on p(eps) == 1/2, i.e. G(eps) == 2, bracketed in
   [0.06, 0.07].  The integrand is singular (but integrable) at the corner
   k = 0, which is what makes the 2-D quadrature hard.  Not the book's AGM
   closed form (6.1), which is the reference method.
     c6.machine   FindRoot + NIntegrate at machine precision;
     c6.wp50      the same with WorkingPrecision -> 60.
   True value and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindRoot", "NIntegrate"}];

g6[e_?NumericQ, wp_] := NIntegrate[
   With[{c = (Cos[k1] + Cos[k2])/2}, (1 - c)/((1 - c)^2 + (2 e Sin[k1])^2)],
   {k1, 0, Pi}, {k2, 0, Pi}, WorkingPrecision -> wp]/Pi^2;

vrun["c6.machine", e /. FindRoot[g6[e, MachinePrecision] == 2, {e, 0.06, 0.07}], 10, 3];
vrun["c6.wp50", e /. FindRoot[g6[e, 60] == 2, {e, 6/100, 7/100},
                              WorkingPrecision -> 60], 48];
