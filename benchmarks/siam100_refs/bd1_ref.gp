\\ Holdout B, family D1 reference, method B (PARI/gp, independent code from
\\ bd1_ref.py): the same automaton power-sum recursion, but K0 = 3 explicit
\\ decimal digits (base 2: 12 bits, base 3: 7 trits), J = 45 moments, 110
\\ digits, the tail by matsolve; bd1.9 by the functional equation of the
\\ fibbinary generating function F(x) = F(x^2) + x F(x^4) (F(0) = 1), which
\\ never enumerates n.   Run:  gp -q bd1_ref.gp
default(parisizemax, 2000000000);
default(realprecision, 110);
\\ step(p, q, d): next state or -1 (rejected); states 0..ns-1; acc(p, q)
step(p, q, d) = {
  if (p == 1, return(if (d == 9, -1, 0)));
  if (p == 2, return(if (d == 0, -1, 0)));
  if (p == 3, return(if (d % 2 == 1, 0, -1)));
  if (p == 4,                                 \\ avoid "314": q = matched prefix length
    if (q == 0, return(if (d == 3, 1, 0)));
    if (q == 1, return(if (d == 1, 2, if (d == 3, 1, 0))));
    if (q == 2, return(if (d == 4, -1, if (d == 3, 1, 0)))));
  if (p == 5, return(if (d == 9, if (q == 0, 1, -1), q)));
  if (p == 6, return(if (d == 7, -1, 0)));
  if (p == 7, return(if (q == 1 && d == 1, -1, d)));
  if (p == 8, return(if (d == 2, -1, 0)));
  if (p == 10, return(if (q == d, -1, d)));
}
nstates(p) = [1, 1, 1, 3, 2, 1, 2, 1, 0, 10][p];
acc(p, q) = if (p == 5, q == 1, 1);
expo(p) = if (p == 6, 2, 1);
baseof(p) = if (p == 7, 2, if (p == 8, 3, 10));
q0(p) = if (p == 10, 10, 0);
run(p) = {
  my(b = baseof(p), s = expo(p), K0 = [3, 12, 7][if (b == 10, 1, if (b == 2, 2, 3))],
     J = 45, ns = nstates(p), nj = J - s + 1, N = ns * nj, S = 0., psi, M, q, dg, x, v, y);
  psi = matrix(ns, nj);
  for (n = 1, b^K0 - 1,
    dg = digits(n, b); q = q0(p);
    for (i = 1, #dg, q = step(p, q, dg[i]); if (q < 0, break));
    if (q < 0, next);
    x = 1. / n;
    if (acc(p, q), S += x^s);
    if (#dg == K0, for (j = s, J, psi[q + 1, j - s + 1] += x^j)));
  M = matrix(N, N);
  for (q = 0, ns - 1, for (d = 0, b - 1,
    my(l = step(p, q, d));
    if (l < 0, next);
    for (j = s, J, for (w = 0, J - j,
      M[l * nj + (j - s) + 1, q * nj + (j + w - s) + 1] +=
        (-1)^w * binomial(j + w - 1, w) * d^w / (1. * b)^(j + w)))));
  v = vector(N, i, psi[(i - 1) \ nj + 1, (i - 1) % nj + 1])~;
  y = matsolve(matid(N) - M, M * v);
  for (q = 0, ns - 1, if (acc(p, q), S += y[q * nj + 1]));
  S;
}
foreach ([1, 2, 3, 4, 5, 6, 7, 8, 10], p, printf("bd1.%d %.95g\n", p, run(p)));
\\ bd1.9: F(x) = sum over fibbinary a >= 0 of x^a = F(x^2) + x F(x^4); answer F(1/2) - 1
F(x) = if (x < 1e-120, 1. + x, F(x^2) + x * F(x^4));
printf("bd1.9 %.95g\n", F(1/2.) - 1);
