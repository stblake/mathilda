(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 225 -- SIAM 100-Digit Challenge HOLDOUT A, family hd20
   (generalises exp 195, App. D Problem 20, "half-exponential").

   Exp 195 asks for h(0) where h(h(z)) = exp(z); its live rows are a finite
   linear system for a functional-equation series (Robbins' truncated Abel
   system), solved by LinearSolve, followed by FindRoot.  This family asks the
   same capability -- truncate a functional equation to a linear system for
   Taylor coefficients, LinearSolve, FindRoot -- on half-iterates whose TRUE
   value is known:

     f(p) = p, f'(p) = lam (0 < lam != 1).  The regular half-iterate at p is
     h = Psi(Sqrt[lam] sigma(z)), sigma the Schroeder function:
        sigma(f(z)) = lam sigma(z), sigma(p) = 0, sigma'(p) = 1.
     Statement: sigma_n(z) = Sum_{k=1}^n c_k z^k (z measured from p), c_1 = 1,
        [z^j] ( sigma_n(f(z)) - lam sigma_n(z) ) / lam^j = 0,   j = 2..n
     (an n x n lower-triangular system, A_jk = ([z^j] f(z)^k - lam delta_jk)
     / lam^j; equation j is divided by lam^j so that the diagonal is O(1) --
     unscaled, the rows span ~lam^n and partial pivoting in a general
     LinearSolve loses 6-10 digits), then
        h(z0) = the root x of sigma_n(x) == Sqrt[lam] sigma_n(z0).
     n is chosen per tier so that the truncation error is < 1e-13 (machine)
     / < 1e-57 (wp50): the rows score against the EXACT half-iterate.

   k  f (fixed point p)            lam    z0     h(z0) (truth.tsv)
   1  16z/(4+3z)         (0)       4      1/2    sigma = 4z/(4-z): 8/9 (Moebius)
   2  tan(4 atan z) =              4      1/10   tan(2 atan z) = 20/99
      (4z-4z^3)/(1-6z^2+z^4)
   3  4z(1-z)  logistic  (0)       4      1/10   sin^2(Sqrt2 asin Sqrt z0)
   4  4z(1+z)            (0)       4      1/10   sinh^2(Sqrt2 asinh Sqrt z0)
   5  Chebyshev T3 at x=1: in u = 1-x, 9   u0=1/10  h = cos(Sqrt3 acos(9/10))
      f(u) = 1 - T3(1-u)                      (the answer is x = 1 - u)
   6  2z/(1+z^2) = tanh(2 atanh z) 2  1/5    tanh(Sqrt2 atanh z0)
   7  2 sin z                      2      1/5    no closed form (Koenigs limit)
   8  2 z e^z                      2      1/10   no closed form (Koenigs limit)
   9  Sqrt2^z at its ATTRACTING fixed point 2 (w = x - 2):
      f(w) = Sqrt2^(w+2) - 2,     Log 2   x0=3/2 no closed form (Koenigs);
      the answer is x = 2 + w; Taylor coefficients involve Log 2
   10 exp 195's own statement for base 2: Robbins' truncated Abel system,
      n = 16 at both tiers, h_16(0) for h(h(z)) = 2^z (a FINITE proxy: its
      truth is the exact n = 16 value)
   Rows hd20.k machine (A in machine numbers, machine FindRoot) and hd20.k
   wp50 (A at 60 digits, FindRoot WorkingPrecision -> 60).  Python: mpmath.taylor coefficients,
   power by truncated convolution (glue), numpy.linalg.solve + brentq /
   mpmath lu_solve + findroot.  Maxima: taylor + linsolve_by_lu (floatfield /
   bigfloatfield) + find_root / bf_find_root.  FriCAS: rows 1-5 only (all
   data rational): exact solve + real-root isolation solve(p = 0, eps);
   rows 6-10 need an irrational Sqrt[lam] / Log 2 in the polynomial, which
   FriCAS's root isolation does not accept: SKIP.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"Series", "CoefficientList", "LinearSolve", "FindRoot"}];

(* the n x n Schroeder matrix at precision prec: A_jk = ([z^j] f^k - lam delta_jk)/lam^j,
   row 1 normalises c_1 = 1.  The Taylor coefficients of f^k come from an exact
   Series (rows 1-8: rational coefficients), then N[A, prec]; for hd20.9, whose
   coefficients are powers of Log 2, from Series of N[f, prec] (an exact Series
   of Sqrt2^z powers carries Log 2 symbolically and costs ~40 s in WL at n = 110;
   hd20.9 is the only row with lam = Log 2). *)
schroederA[f_, lam_, n_, prec_] := Module[{pw, fn = If[lam === Log[2], N[f, prec], f]},
  pw = Table[PadRight[CoefficientList[Normal[Series[fn^k, {z, 0, n}]], z], n + 1], {k, n}];
  N[Table[If[j == 1, KroneckerDelta[k, 1], (pw[[k, j + 1]] - lam KroneckerDelta[j, k])/lam^j], {j, n}, {k, n}], prec]];

halfIter[f_, lam_, z0_, n_] := Module[{c},
  c = LinearSolve[schroederA[f, lam, n, MachinePrecision], UnitVector[n, 1]];
  x /. FindRoot[c . x^Range[n] == Sqrt[N[lam]] c . N[z0]^Range[n], {x, N[z0]}]];

halfIterWP[f_, lam_, z0_, n_] := Module[{c},
  c = LinearSolve[schroederA[f, lam, n, 60], UnitVector[n, 1]];
  x /. FindRoot[c . x^Range[n] == Sqrt[N[lam, 60]] c . z0^Range[n], {x, z0}, WorkingPrecision -> 60]];

robbinsH0[b_, n_, wp_] := Module[{A, v},
  A = Table[(k Log[b])^j/j! - KroneckerDelta[j, k], {j, 0, n - 1}, {k, 1, n}];
  v = LinearSolve[If[wp === MachinePrecision, N[A], N[A, wp]], UnitVector[n, 1]];
  x /. If[wp === MachinePrecision, FindRoot[-1 + v . x^Range[n] == -1/2, {x, 1/2}],
     FindRoot[-1 + v . x^Range[n] == -1/2, {x, 1/2}, WorkingPrecision -> wp]]];

(* {f, lam, z0, n machine, n wp50, post} *)
probs = {
  {16 z/(4 + 3 z), 4, 1/2, 25, 90, Identity},
  {(4 z - 4 z^3)/(1 - 6 z^2 + z^4), 4, 1/10, 25, 90, Identity},
  {4 z (1 - z), 4, 1/10, 25, 90, Identity},
  {4 z (1 + z), 4, 1/10, 25, 90, Identity},
  {1 - (4 (1 - z)^3 - 3 (1 - z)), 9, 1/10, 25, 90, (1 - #) &},
  {2 z/(1 + z^2), 2, 1/5, 30, 110, Identity},
  {2 Sin[z], 2, 1/5, 30, 90, Identity},
  {2 z Exp[z], 2, 1/10, 30, 110, Identity},
  {Sqrt[2]^(z + 2) - 2, Log[2], -1/2, 30, 110, (2 + #) &}};

Do[With[{p = probs[[k]], lab = "hd20." <> ToString[k]},
   vrun[lab <> " machine", p[[6]][halfIter[p[[1]], p[[2]], p[[3]], p[[4]]]], 10, 1, 60];
   vrun[lab <> " wp50", p[[6]][halfIterWP[p[[1]], p[[2]], p[[3]], p[[5]]]], 48, 1, 120]],
  {k, Length[probs]}];
vrun["hd20.10 machine", robbinsH0[2, 16, MachinePrecision], 10, 1, 60];
vrun["hd20.10 wp50", robbinsH0[2, 16, 60], 48, 1, 60];
