(* bench-timeout: 500 *)
(* ==========================================================================
   Experiment 168 -- SIAM 100-Digit Challenge, Problem 3 ("How Far Away Is
   Infinity?").

       The infinite matrix A with a_ij = 1/((i+j-1)(i+j-2)/2 + j)
       (a11 = 1, a12 = 1/2, a21 = 1/3, a13 = 1/4, a22 = 1/5, a31 = 1/6, ...)
       is a bounded operator on l^2.  What is ||A|| ?  =  1.2742241528...

   Idiomatic statement: the 2-norm (Norm, the largest singular value) of the
   n x n section A_n, for n large enough.  ||A|| - ||A_n|| ~ n^-3 (measured:
   the error falls by ~8 per doubling: 5.5e-10 at n = 2048, 6.9e-11 at
   n = 4096), so n = 4096 is the smallest power of two whose section carries
   the challenge's 10 digits; ||A_4096|| = 1.27422415275... has 10 correct
   digits of ||A||.  No system has a public built-in for extrapolating a
   list of section norms (WL's SequenceLimit is only a Method name of
   NLimit/NSum, and NLimit samples n = 16, 160, 1600, ... -- sections of
   size 10^5 and up), so the statement is the plain truncation.  Rows:
     c3.machine        Norm of A_4096 in doubles, scored against ||A||
                       (need 10: what the section can give);
     c3.wp50           Norm of A_4096 at 60 digits, scored against ||A||:
                       an INFEASIBLE tier -- no truncation reaches 50 digits
                       (that needs n ~ 10^16), and the section's own norm at
                       60 digits is a 4096^2 MPFR SVD; capped at 60 s.
                       Truth known to 40 digits, need 38;
     c3.machine-n500   the reduced instance ||A_500|| in doubles (a
                       well-defined number, known to 70 digits);
     c3.wp50-n100      the reduced instance ||A_100|| at 60 digits.
   Maxima's SVD is LAPACK dgesvd (machine only; bigfloat rows SKIP); FriCAS
   has no numerical SVD (SKIP).  True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"Norm"}];

(* the section, built directly in the working precision *)
secA[n_, prec_] := Table[N[1/((i + j - 1) (i + j - 2)/2 + j), prec], {i, n}, {j, n}];

vrun["c3.machine", Norm[secA[4096, MachinePrecision]], 10];
vrun["c3.wp50", Norm[secA[4096, 60]], 38, 1, 60];
vrun["c3.machine-n500", Norm[secA[500, MachinePrecision]], 10, 3];
vrun["c3.wp50-n100", Norm[secA[100, 60]], 48];
