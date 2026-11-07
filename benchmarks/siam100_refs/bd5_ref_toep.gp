\\ Holdout B, family D5 reference, method B for the NONNORMAL TOEPLITZ rows
\\ (PARI/gp, 110 digits; independent code from bd5_ref_toep.py): the same KKT
\\ characterisation (M(c) = p(A) J symmetric, eigenvalues +- singular values;
\\ 'pair': lmax = -lmin with dual weight w in (0,1); 'smooth': grad lmax = 0),
\\ but eigenpairs from PARI's qfjacobi and the Jacobian by CENTRAL differences
\\ (h = 1e-40), started from the stage-1 point bd5_<p>.start.
\\ Run: gp -q bd5_ref_toep.gp   (after bd5_ref_toep.py p for p = 3, 6, 7, 8)
default(parisizemax, 2000000000);
default(realprecision, 110);
pr(name, v) = printf("%s %.95g\n", name, v);
toep(n, L) = matrix(n, n, i, j, my(r = 0); for (t = 1, #L, if (L[t][1] == j - i, r = L[t][2])); r);
GR = [[-1, -1], [0, 1], [1, 1], [2, 1], [3, 1]];
ext(cc) = {
  my(M = gMj[gk + 1] + sum(j = 1, gk, cc[j] * gMj[j]), E, im = 1, in = 1);
  M = (M + M~) / 2; E = qfjacobi(M);
  for (i = 2, gn, if (E[1][i] > E[1][im], im = i); if (E[1][i] < E[1][in], in = i));
  [E[1][im], E[2][, im], E[1][in], E[2][, in]];
}
F(y) = {
  my(e = ext(y), gM, gm, k = gk);
  gM = vector(k, j, e[2]~ * gMj[j] * e[2]); gm = vector(k, j, -(e[4]~ * gMj[j] * e[4]));
  if (gkind == "pair",
    [concat([e[1] + e[3]], vector(k, j, y[k + 1] * gM[j] + (1 - y[k + 1]) * gm[j])), max(e[1], -e[3])],
    if (gkind == "max", [gM, e[1]], [gm, -e[3]]));
}
solvecase(p, n, L, k, kind) = {
  my(A = toep(n, L) * 1., J = matrix(n, n, i, j, i + j == n + 1), Pk = matid(n),
     c0 = readvec(Str("bd5_", p, ".start")), x, m, f0, Jac, h = 1e-40, nf, lam);
  gn = n; gk = k; gkind = kind; gMj = vector(k + 1);
  for (j = 1, k + 1, gMj[j] = Pk * J; Pk = Pk * A);
  x = vector(k, j, c0[j]); if (kind == "pair", x = concat(x, [0.5]));
  m = #x;
  for (it = 0, 30,
    f0 = F(x); nf = sqrt(f0[1] * f0[1]~); lam = f0[2];
    printf("# Newton %d |F| = %.3g  norm = %.50g\n", it, nf, lam);
    if (nf < 1e-90, break);
    Jac = matrix(m, m);
    for (kk = 1, m, my(xp = x, xm = x, fp, fm);
      xp[kk] += h; xm[kk] -= h; fp = F(xp)[1]; fm = F(xm)[1];
      for (i = 1, m, Jac[i, kk] = (fp[i] - fm[i]) / (2 * h)));
    x = x - (matsolve(Jac, f0[1]~))~);
  if (kind == "pair", printf("# dual weight w = %.20g\n", x[m]));
  pr(Str("bd5.", p), lam);
}
\\ structures as found by stage 1 (bd5_ref_toep.py): 'max' = smooth, top eigenvalue
solvecase(7, 48, GR, 1, "max");
solvecase(6, 40, [[-1, 1], [1, 1/4]], 3, "max");
solvecase(10, 32, [[-1, 1], [1, 1], [2, 1]], 3, "max");
solvecase(3, 48, GR, 4, "max");
solvecase(8, 24, GR, 6, "pair");
