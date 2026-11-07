\\ Holdout B, family D10 -- method B: PARI/gp, realprecision 100.  Solve the
\\ Poisson equation  h(r) - sum_i p_i h(a_i + b_i/r) = log r - lambda  by
\\ collocation at N Chebyshev points of the FIRST kind on the invariant
\\ interval (barycentric interpolation matrix), with sum h = 0; lambda is
\\ the extra unknown.  (Method A, bd10_ref.py, is power iteration at
\\ second-kind points.)  Answer exp(lambda).
default(realprecision, 100);
interval(cs) = {
  my(lo = 1., hi = 1., ims, nlo, nhi);
  for (it = 1, 4000,
    ims = concat(vector(#cs, i, [cs[i][2] + cs[i][3]/lo, cs[i][2] + cs[i][3]/hi]));
    nlo = vecmin(ims); nhi = vecmax(ims);
    if (abs(nlo - lo) + abs(nhi - hi) < 1e-95, break);
    lo = nlo; hi = nhi);
  [lo, hi];
}
lyap(cs, N) = {
  my(I = interval(cs), lo = I[1], hi = I[2], th, xs, w, A, rhs, y, d, row, s);
  th = vector(N, k, (2*k - 1)*Pi/(2*N));
  xs = vector(N, k, (lo + hi)/2 + (hi - lo)/2*cos(th[k]));
  w = vector(N, k, (-1)^k * sin(th[k]));
  A = matrix(N + 1, N + 1); rhs = vectorv(N + 1);
  for (j = 1, N,
    A[j, j] += 1; A[j, N + 1] = 1; rhs[j] = log(xs[j]);
    for (i = 1, #cs,
      y = cs[i][2] + cs[i][3]/xs[j];
      row = vector(N, k, w[k]/(y - xs[k])); s = vecsum(row);
      for (k = 1, N, A[j, k] -= cs[i][1] * row[k]/s)));
  for (k = 1, N, A[N + 1, k] = 1);
  exp(matsolve(A, rhs)[N + 1]);
}
pr(lab, v) = printf("%s %s\n", lab, Strprintf("%.70g", v));
N = 110;
pr("bd10.2", lyap([[1/2, 3, 1], [1/2, 3, -1]], N));
pr("bd10.3", lyap([[1/3, 3, 1], [2/3, 3, -1]], N));
pr("bd10.4", lyap([[1/2, 1, 1], [1/2, 1, 2]], N));
pr("bd10.5", lyap([[1/3, 3, -1], [1/3, 3, 0], [1/3, 3, 1]], N));
pr("bd10.6", lyap([[1/2, 4, 2], [1/2, 4, -2]], N));
pr("bd10.7", lyap([[1/2, 1, 1], [1/2, 2, 1]], N));
pr("bd10.8", lyap([[1/2, 5, 3], [1/2, 5, -3]], N));
pr("bd10.9", lyap([[1/2, 3, 2], [1/2, 1, 2]], N));
pr("bd10.10", lyap([[1/3, 1, 1], [1/3, 2, 1], [1/3, 3, 1]], N));
pr("bd10.1", (1*2*3*4)^(1/4));
