(* Holdout A, family hd2, method B (independent of hd2_ref.py's exact
   fixed-point anchor): direct summation of the first b^K - 1 terms, the
   digit-length blocks K < k <= kmax summed with exact block sums
   (PolyGamma / HurwitzZeta), and the remainder k > kmax replaced by
   Log[b] * (a tail of the next-level series) -- the dropped O(b^-kmax) block
   corrections are below 10^-75.  When the next level is the SAME series the
   remainder is S minus a known partial sum, and S solves a linear equation.
   Run: wolframscript -file hd2_ref.wl *)
$MaxExtraPrecision = 500;
prec = 100;
(* {b, c, seed, w, e, q, shift, p, weight} *)
probs = <|
  0 -> {2, 2, {1, 2}, 1, 0, 1, 0, 1, "all"},
  1 -> {2, 2, {1, 2}, 2, 0, 1, 0, 1, "all"},
  2 -> {2, 2, {1, 3}, 1, 0, 1, 0, 1, "all"},
  3 -> {2, 2, {1, 2}, 1, 0, 1, 0, 1, "alt"},
  4 -> {10, 9, Range[9], 1, 1, 1, 0, 1, "all"},
  5 -> {3, 2, {1, 2}, 1, 0, 2, 0, 1, "all"},
  6 -> {2, 2, {1, 2}, 1, 0, 1, 0, 2, "all"},
  7 -> {2, 1, {1}, 1, 0, 1, 1, 1, "all"},
  8 -> {10, 9, Range[9], 1, 1, 1, 0, 1, "alt"},
  9 -> {2, 2, {1, 2}, 1, 1/2, 1, 0, 1, "all"},
  10 -> {2, 2, {1, 2}, 1, 0, 1, 0, 1, "odd"}|>;

solve[{b_, c_, seed_, w_, e_, q_, sh_, p0_, wt_}] := Module[
  {f, L, blk, T, K, kmax, Sv, sol, top, pg, hz},
  L[n_] := IntegerLength[n, b] - sh;
  f[n_] := f[n] = If[n <= c, seed[[n]], n^w L[n]^e f[L[n]]^q];
  (* block sum of n^-sig over b^(k-1) <= n <= b^k - 1 in class cl *)
  (* numeric arguments: PolyGamma of an exact integer would expand to a HarmonicNumber *)
  pg[x_] := PolyGamma[N[x, prec]]; hz[sg_, x_] := HurwitzZeta[sg, N[x, prec]];
  blk[sig_, k_, cl_] := Module[{a = b^(k - 1), z = b^k - 1},
    Switch[cl,
      "all", If[sig == 1, pg[z + 1] - pg[a],
                 hz[sig, a] - hz[sig, z + 1]],
      "odd", Module[{m0 = Ceiling[(a - 1)/2], m1 = Floor[(z - 1)/2]},
                If[sig == 1, (pg[m1 + 3/2] - pg[m0 + 1/2])/2,
                   (hz[sig, m0 + 1/2] - hz[sig, m1 + 3/2])/2^sig]],
      "even", Module[{m0 = Ceiling[a/2], m1 = Floor[z/2]},
                If[sig == 1, (pg[m1 + 1] - pg[m0])/2,
                   (hz[sig, m0] - hz[sig, m1 + 1])/2^sig]]]];
  (* T[p, s, N, unknown]: Sum_{n >= N} n^-s f(n)^-p as (number) + coeff * unknown,
     where unknown stands for the top-level total when the recursion closes *)
  T[p_, s_, n0_, depth_] := Module[{K0, km, sig, dir, bl, rem},
    sig = s + w p;
    K0 = Max[IntegerLength[n0, b] + 2, If[b == 2, 12, If[b == 3, 8, 4]]];
    km = If[sig == 1, If[b == 2, 260, If[b == 3, 170, 90]], Ceiling[90/((sig - 1) Log10[b])] + K0 + 5];
    dir = Sum[N[n^-s f[n]^-p, prec], {n, n0, b^(K0 - 1) - 1}];
    bl = Sum[N[(k - sh)^(-e p) f[k - sh]^(-q p) blk[sig, k, "all"], prec], {k, K0, km}];
    If[sig != 1, dir + bl,
      (* remainder: Log[b] Sum_{k > km} (k - sh)^(-ep) f(k - sh)^(-qp) *)
      If[{q p, e p} === {p, s} && depth > 0,
        dir + bl + N[Log[b], prec] (Sv - Sum[N[n^-s f[n]^-p, prec], {n, 1, km - sh}]),
        dir + bl + N[Log[b], prec] T[q p, e p, km + 1 - sh, depth + 1]]]];
  top = If[wt === "all",
     Sum[N[f[n]^-p0, prec], {n, 1, c}] + T[p0, 0, c + 1, 1],
     Module[{K0 = If[b == 2, 12, 4], km, cls, sg, dir, bl, sig = w p0},
       sg[n_] := If[wt === "alt", (-1)^(n + 1), Mod[n, 2]];
       km = If[b == 2, 260, 90];
       dir = Sum[N[sg[n] f[n]^-p0, prec], {n, 1, b^(K0 - 1) - 1}];
       bl = Sum[N[(k - sh)^(-e p0) f[k - sh]^(-q p0) If[wt === "alt",
                blk[sig, k, "odd"] - blk[sig, k, "even"], blk[sig, k, "odd"]], prec], {k, K0, km}];
       dir + bl + If[wt === "alt", 0, N[Log[b]/2, prec] (S0 - Sum[N[f[n]^-1, prec], {n, 1, km - sh}])]]];
  sol = Sv /. Solve[Sv == top, Sv];
  First[sol]];
(* the odd-n series (problem 10) recurses into the plain series of problem 0 *)
S0 = solve[probs[0]];
Do[Print[k, " ", NumberForm[solve[probs[k]], 70]], {k, 0, 10}];
