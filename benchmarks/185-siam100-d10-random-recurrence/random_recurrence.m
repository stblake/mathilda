(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 185 -- SIAM 100-Digit Challenge, Appendix D Problem 10
   ("A random recurrence"; problem 12 on the book's web page).

     x_0 = x_1 = 1,  x_{n+1} = 2 x_n +- x_{n-1}, each sign independently +
     or - with probability 1/2.  To what does |x_n|^(1/n) converge almost
     surely?                         Answer: 1.9159489091701204381698...

   No system has a built-in that computes a Lyapunov exponent of a random
   matrix product, so the direct statement is the experiment itself -- a
   Monte-Carlo estimate, accurate to O(n^-1/2) whatever the arithmetic:
     machine : n = 10^6 random signs (SeedRandom[1]); the ratio
               r_k = x_{k+1}/x_k obeys r_k = 2 +- 1/r_{k-1}, r_0 = 1, so
               |x_n|^(1/n) = Exp[Mean[Log[r_0 .. r_{n-1}]]]  (x_n itself
               overflows a double near n = 1090);
     wp50    : the literal |x_n|^(1/n) with x_n an EXACT integer, n = 10^4
               (x_n has ~2800 digits), evaluated by N[., 60].  Sixty-digit
               arithmetic does not help: the statistical error is ~1e-2.
   Expected correct digits for every system: ~3 (machine) and ~1-2 (wp50),
   seed-dependent; the RNG streams differ across systems.  These rows time
   the loop/bignum machinery; the digits column is honest, not a target.
   The reference (truth.tsv) is the Furstenberg integral
   lambda = Int log r d nu(r) over the stationary law nu of the ratio chain,
   computed two ways (siam100_refs/d10_ref.{wl,py}).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"RandomChoice", "FoldList", "Fold"}];

lyapM[] := (SeedRandom[1];
  Exp[Mean[Log[FoldList[2. + #2/#1 &, 1., RandomChoice[{-1, 1}, 10^6 - 1]]]]]);
lyap50[] := (SeedRandom[1];
  N[Abs[Last[Fold[{#1[[2]], 2 #1[[2]] + #2 #1[[1]]} &, {1, 1},
                   RandomChoice[{-1, 1}, 10^4 - 1]]]]^(1/10^4), 60]);

vrun["d10.machine", lyapM[], 2, 3];
vrun["d10.wp50", lyap50[], 1, 3];
