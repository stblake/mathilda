\\ Holdout B, family C5 (experiment 230) -- route B (PARI/gp), realprecision 120.
\\ Real: solve the levelled system  p(x_r) + s_r E = f(x_r)  on the LP's
\\ alternation set (signs s_r from the LP residual), then check
\\ max_k |f(x_k) - p(x_k)| <= E (1 + 1e-100) over all M points.
\\ Complex: Newton (derivnum Jacobian) on |f(z_k) - p(z_k)|^2 = E^2 over the
\\ SLSQP active set, real coefficients, then the same all-points check.
\p 120
read("bc5_data.gp");
setv(w, j, t) = my(u = w); u[j] = t; u;
{
for (q = 1, #RJ,
  my([f, a, b, n, M, alt, sg] = RJ[q], xs, A, rhs, s, c, E, mx);
  xs = vector(M, k, a + (b - a) * (k - 1) / (M - 1));
  A = matrix(n + 2, n + 2, i, j, if (j <= n + 1, xs[alt[i] + 1]^(j - 1) * 1., sg[i]));
  rhs = vector(n + 2, i, f(xs[alt[i] + 1] * 1.))~;
  s = matsolve(A, rhs); c = s[1..n+1]; E = abs(s[n + 2]);
  mx = vecmax(vector(M, k, abs(f(xs[k] * 1.) - sum(j = 1, n + 1, c[j] * (xs[k] * 1.)^(j - 1)))));
  printf("bc5.%d B E = %.80g\n        max_all/E - 1 = %.3g\n", q, E, mx / E - 1));
for (q = 1, #CJ,
  my([f, n, M, act, x0] = CJ[q], zs, fs, e2, G, u, J, d, mx, E);
  zs = vector(M, k, exp(2 * Pi * I * (k - 1) / M));
  fs = vector(M, k, f(zs[k]));
  e2 = ((cf, k) -> norm(fs[k + 1] - sum(j = 1, n + 1, cf[j] * zs[k + 1]^(j - 1))));
  G = (w -> vector(#act, i, e2(w[1..n+1], act[i]) - w[n + 2]^2));
  u = x0 * 1.;
  for (it = 1, 80,
    J = matrix(#act, n + 2, i, j, derivnum(t = u[j], G(setv(u, j, t))[i]));
    d = matsolve(J, G(u)~); u -= d~;
    if (normlp(d, oo) < 1e-105, break));
  E = abs(u[n + 2]);
  mx = vecmax(vector(M, k, sqrt(e2(u[1..n+1], k - 1))));
  printf("bc5.%d B E = %.80g\n        max_all/E - 1 = %.3g\n", #RJ + q, E, mx / E - 1));
}
quit
