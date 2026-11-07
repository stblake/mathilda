(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 177 -- SIAM 100-Digit Challenge, Appendix D Problem 2
   (David Smith; cf. Putnam 2002 A6).

       S = Sum_{n>=1} 1/f(n),   f(1) = 1,  f(2) = 2,
                                f(n) = n f(d)  (n > 2),  d = number of base-2
                                                         digits of n,
       S = 2.0863776650...

   The series converges, but extremely slowly: the terms of bit length k sum
   to about log 2 / f(k), and f(k) grows like k log k log log k ..., so the
   partial sum to n = 10^100 is still wrong in the second digit.  No
   series-acceleration method that sees only the terms can get it.  (The
   closed reference groups by bit length and turns the self-similarity into
   a linear equation for S; see truth.tsv and siam100_refs/d2_ref.*.)

   Idiomatic statement only: NSum of the direct term, at machine precision
   and with WorkingPrecision 60.  f is the definition as written (glue).
   Every system is expected to return a confident wrong number here; that is
   the measurement.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSum", "BitLength"}];

f[n_Integer] := If[n <= 2, n, n f[BitLength[n]]];

vrun["d2.machine", NSum[1/f[n], {n, 1, Infinity}], 10, 3];
vrun["d2.wp50", NSum[1/f[n], {n, 1, Infinity}, WorkingPrecision -> 60], 48];
