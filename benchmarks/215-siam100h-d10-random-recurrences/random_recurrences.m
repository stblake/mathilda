(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 215 -- SIAM 100-Digit Challenge HOLDOUT A, family D10
   (generalises experiment 185, "a random recurrence").

   The capability: the almost-sure growth rate lim |x_n|^(1/n) of a random
   linear recurrence (the top Lyapunov exponent of a random matrix product),
   x_0 = x_1 = 1.  Ten instances, varied in coefficients, sign law and order:

     hd10.1   x_{n+1} = c x_n,            c in {1, 2}          = Sqrt[2]  (first order; easy)
     hd10.2   x_{n+1} = 3 x_n +- x_{n-1}
     hd10.3   x_{n+1} = 4 x_n +- x_{n-1}
     hd10.4   x_{n+1} = (5/2) x_n +- x_{n-1}                   (rational x_n)
     hd10.5   x_{n+1} = 3 x_n +- 2 x_{n-1}
     hd10.6   x_{n+1} = c x_n + x_{n-1},  c in {1, 2}          (random continued fraction)
     hd10.7   x_{n+1} = 3 x_n + s x_{n-1}, s in {-1, 0, 1}
     hd10.8   x_{n+1} = 3 x_n +- x_{n-1}, P(+) = 3/4           (biased signs)
     hd10.9   x_{n+1} = c x_n + x_{n-1},  c in {1, 2, 3}
     hd10.10  x_{n+1} = 2 x_n +- x_{n-1}, P(+) = 3/4           (the original, biased;
                                                               marginal fixed point r = 1)

   Every choice is independent and uniform unless a probability is given.
   As in 185 no system has a built-in for a Lyapunov exponent, so the direct
   statement is the Monte-Carlo experiment itself:
     machine : n = 10^6 random choices (SeedRandom[1]); the ratio
               r_k = x_{k+1}/x_k obeys r_k = c + s/r_{k-1}, r_0 = 1, and
               |x_n|^(1/n) = Exp[Mean[Log[r_0 .. r_{n-1}]]];
     wp50    : the literal |x_n|^(1/n) with x_n EXACT (integer or rational),
               n = 10^4, evaluated by N[., 60].
   Expected correct digits for every system: ~2-3 (machine) and ~1-2 (wp50),
   seed-dependent; the RNG streams differ across systems.  These rows time the
   loop/bignum machinery; the digits column is honest, not a target.
   The reference (truth.tsv) is the Furstenberg integral over the stationary
   law of the ratio chain, computed two ways (siam100_refs/hd10_ref.{wl,py}).
   Caps: 60 s per run (warm-up + 3 timed).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"RandomChoice", "FoldList", "Fold"}];

nM = 10^6; n5 = 10^4;
(* machine: ratio chain r <- f[r, choice];  wp50: exact pair {x_{k-1}, x_k} *)
lyapM[f_, ch_] := (SeedRandom[1]; Exp[Mean[Log[FoldList[f, 1., RandomChoice[ch, nM - 1]]]]]);
lyap50[g_, ch_] := (SeedRandom[1];
  N[Abs[Last[Fold[{#1[[2]], g[#1, #2]} &, {1, 1}, RandomChoice[ch, n5 - 1]]]]^(1/n5), 60]);

vrun["hd10.1 machine", (SeedRandom[1]; Exp[Mean[Log[RandomChoice[{1., 2.}, nM]]]]), 2, 3, 60];
vrun["hd10.1 wp50", (SeedRandom[1]; N[Abs[Fold[#2 #1 &, 1, RandomChoice[{1, 2}, n5]]]^(1/n5), 60]), 1, 3, 60];
vrun["hd10.2 machine", lyapM[3. + #2/#1 &, {-1, 1}], 2, 3, 60];
vrun["hd10.2 wp50", lyap50[3 #1[[2]] + #2 #1[[1]] &, {-1, 1}], 1, 3, 60];
vrun["hd10.3 machine", lyapM[4. + #2/#1 &, {-1, 1}], 2, 3, 60];
vrun["hd10.3 wp50", lyap50[4 #1[[2]] + #2 #1[[1]] &, {-1, 1}], 1, 3, 60];
vrun["hd10.4 machine", lyapM[2.5 + #2/#1 &, {-1, 1}], 2, 3, 60];
vrun["hd10.4 wp50", lyap50[5/2 #1[[2]] + #2 #1[[1]] &, {-1, 1}], 1, 3, 60];
vrun["hd10.5 machine", lyapM[3. + 2 #2/#1 &, {-1, 1}], 2, 3, 60];
vrun["hd10.5 wp50", lyap50[3 #1[[2]] + 2 #2 #1[[1]] &, {-1, 1}], 1, 3, 60];
vrun["hd10.6 machine", lyapM[#2 + 1/#1 &, {1., 2.}], 2, 3, 60];
vrun["hd10.6 wp50", lyap50[#2 #1[[2]] + #1[[1]] &, {1, 2}], 1, 3, 60];
vrun["hd10.7 machine", lyapM[3. + #2/#1 &, {-1, 0, 1}], 2, 3, 60];
vrun["hd10.7 wp50", lyap50[3 #1[[2]] + #2 #1[[1]] &, {-1, 0, 1}], 1, 3, 60];
vrun["hd10.8 machine", lyapM[3. + #2/#1 &, {1, 3} -> {-1, 1}], 2, 3, 60];
vrun["hd10.8 wp50", lyap50[3 #1[[2]] + #2 #1[[1]] &, {1, 3} -> {-1, 1}], 1, 3, 60];
vrun["hd10.9 machine", lyapM[#2 + 1/#1 &, {1., 2., 3.}], 2, 3, 60];
vrun["hd10.9 wp50", lyap50[#2 #1[[2]] + #1[[1]] &, {1, 2, 3}], 1, 3, 60];
vrun["hd10.10 machine", lyapM[2. + #2/#1 &, {1, 3} -> {-1, 1}], 2, 3, 60];
vrun["hd10.10 wp50", lyap50[2 #1[[2]] + #2 #1[[1]] &, {1, 3} -> {-1, 1}], 1, 3, 60];
