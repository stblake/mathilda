(* bench-timeout: 1200 *)
(* ==========================================================================
   Experiment 245 -- SIAM 100-Digit Challenge, SEALED HOLDOUT B for Appendix D
   Problem 10 (experiment 185, "a random recurrence").

   x_0 = x_1 = 1,  x_{n+1} = a x_n + b x_{n-1}, with the pair (a, b) drawn
   independently at every step from a finite list with the given weights.
   To what does |x_n|^(1/n) converge almost surely?  Ten new recurrences:

     bd10.1   (a,b) in {(1,0),(2,0),(3,0),(4,0)}, equal weights:
              x_{n+1} = c x_n, the limit is 24^(1/4)          (easy, closed form)
     bd10.2   {(3, 1), (3,-1)}           x_{n+1} = 3 x_n +- x_{n-1}
     bd10.3   {(3, 1), (3,-1)} weights 1:2  (biased signs)
     bd10.4   {(1, 1), (1, 2)}           x_{n+1} = x_n + (1 or 2) x_{n-1}
     bd10.5   {(3,-1), (3, 0), (3, 1)}   three-way choice
     bd10.6   {(4, 2), (4,-2)}
     bd10.7   {(1, 1), (2, 1)}           random continued-fraction digits 1, 2
     bd10.8   {(5, 3), (5,-3)}
     bd10.9   {(3, 2), (1, 2)}           random leading coefficient
     bd10.10  {(1, 1), (2, 1), (3, 1)}   digits 1, 2, 3

   As in experiment 185, no system has a built-in for the Lyapunov exponent of
   a random matrix product, so the direct statement is the Monte-Carlo
   experiment itself:
     machine : n = 10^6 draws (SeedRandom[1]); the ratio r_k = x_{k+1}/x_k
               obeys r_k = a + b/r_{k-1}, r_0 = 1, so |x_n|^(1/n) =
               Exp[Mean[Log[r_0 .. r_{n-1}]]] (x_n itself overflows a double);
     wp50    : the literal |x_n|^(1/n) with x_n an EXACT integer, n = 10^4,
               evaluated by N[., 60].
   The digits are limited by the O(n^-1/2) statistical error for ANY system
   and seed (need = what 4 standard errors allow; truth.tsv); the rows time
   the loop / bignum machinery.  The reference is the Furstenberg integral
   lambda = Int log r dnu(r) over the stationary law of the ratio chain
   (siam100_refs/bd10_ref.{py,gp}).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"RandomChoice", "FoldList", "Fold"}];

cs10[1] = {{1, 1, 1, 1}, {{1, 0}, {2, 0}, {3, 0}, {4, 0}}};
cs10[2] = {{1, 1}, {{3, 1}, {3, -1}}};
cs10[3] = {{1, 2}, {{3, 1}, {3, -1}}};
cs10[4] = {{1, 1}, {{1, 1}, {1, 2}}};
cs10[5] = {{1, 1, 1}, {{3, -1}, {3, 0}, {3, 1}}};
cs10[6] = {{1, 1}, {{4, 2}, {4, -2}}};
cs10[7] = {{1, 1}, {{1, 1}, {2, 1}}};
cs10[8] = {{1, 1}, {{5, 3}, {5, -3}}};
cs10[9] = {{1, 1}, {{3, 2}, {1, 2}}};
cs10[10] = {{1, 1, 1}, {{1, 1}, {2, 1}, {3, 1}}};

lyM[{w_, cs_}] := (SeedRandom[1];
  Exp[Mean[Log[FoldList[#2[[1]] + #2[[2]]/#1 &, 1., RandomChoice[w -> cs, 10^6 - 1]]]]]);
ly50[{w_, cs_}] := (SeedRandom[1];
  N[Abs[Last[Fold[{#1[[2]], #2[[1]] #1[[2]] + #2[[2]] #1[[1]]} &, {1, 1},
                   RandomChoice[w -> cs, 10^4 - 1]]]]^(1/10^4), 60]);

needM = {2, 3, 3, 3, 3, 3, 3, 3, 2, 2};
need50 = {1, 2, 2, 2, 2, 2, 2, 2, 1, 1};
Do[With[{c = cs10[k], dm = needM[[k]], dw = need50[[k]]},
  vrun["bd10." <> ToString[k] <> " machine", lyM[c], dm, 3, 60];
  vrun["bd10." <> ToString[k] <> " wp50", ly50[c], dw, 3, 60]],
  {k, 1, 10}];
