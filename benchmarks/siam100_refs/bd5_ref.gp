\\ Holdout B, family D5 reference, method B for the SYMMETRIC and SMOOTH rows
\\ (PARI/gp, 110 digits): independent of bd5_ref_sym.py in the eigenvalues and
\\ in the ALGORITHM -- a discrete Remez single-exchange iteration (no LP) for
\\ min over monic p of max_i |p(l_i)|, which stops when the levelled error on
\\ the reference equals the max over all eigenvalues (alternation theorem).
\\ Eigenvalues: bd5.1 cos(j pi/6); bd5.4 polroots of the characteristic
\\ polynomial of the Kac matrix (off-diagonal sqrt(j (11 - j)), 110 digits);
\\ bd5.5 polroots of charpoly of tridiag(-1, 2, -1); bd5.9 polroots of the
\\ EXACT rational characteristic polynomial of the Hilbert matrix.
\\ bd5.2: exact rational normal equations by matsolve.   Run: gp -q bd5_ref.gp
default(parisizemax, 2000000000);
default(realprecision, 110);
pr(name, v) = printf("%s %.95g\n", name, v);
remez(L, k) = {
  my(m = #L, R, M, x, H, e, im, pos, s, err);
  L = vecsort(L);
  R = vector(k + 1, t, 1 + round((t - 1) * (m - 1) / k));
  for (it = 1, 200,
    M = matrix(k + 1, k + 1, t, j, if (j <= k, L[R[t]]^(j - 1), -(-1)^(t - 1)));
    x = matsolve(M, vector(k + 1, t, -L[R[t]]^k)~);
    H = x[k + 1];
    e = vector(m, i, L[i]^k + sum(j = 1, k, x[j] * L[i]^(j - 1)));
    im = 1; for (i = 2, m, if (abs(e[i]) > abs(e[im]), im = i));
    err = abs(e[im]) - abs(H);
    if (err <= abs(H) * 1e-100, return(abs(H)));
    if (setsearch(Set(R), im), return(abs(H)));
    s = sign(e[im]);
    \\ single exchange keeping alternation
    pos = 0; for (t = 1, k + 1, if (R[t] < im, pos = t));
    if (pos == 0,
      if (sign(e[R[1]]) == s, R[1] = im, R = concat([im], R[1..k])),
    if (pos == k + 1,
      if (sign(e[R[k + 1]]) == s, R[k + 1] = im, R = concat(R[2..k + 1], [im])),
      if (sign(e[R[pos]]) == s, R[pos] = im, R[pos + 1] = im))));
  error("remez did not converge");
}
pr("bd5.1", remez(vector(7, j, cos((j - 1) * Pi / 6)), 6));
{ my(n = 11, K = matrix(n, n)); for (j = 1, n - 1, K[j, j + 1] = K[j + 1, j] = sqrt(j * (n - j)));
  pr("bd5.4", remez(real(polroots(charpoly(K))), 4)); }
{ my(n = 20, T = matrix(n, n, i, j, if (i == j, 2, if (abs(i - j) == 1, -1, 0))));
  pr("bd5.5", remez(real(polroots(charpoly(T))), 4)); }
pr("bd5.9", remez(real(polroots(charpoly(mathilbert(8)))), 3));
\\ bd5.2: Grcar 48, degree 8, Frobenius norm, exact
ip(X, Y) = trace(X~ * Y);
GA = matrix(48, 48, i, j, if (i == j || (j - i >= 1 && j - i <= 3), 1, if (j == i - 1, -1, 0)));
GP = vector(9); GP[1] = matid(48); for (j = 2, 9, GP[j] = GP[j - 1] * GA);
GG = matrix(8, 8, i, j, ip(GP[i], GP[j])); Gb = vector(8, i, ip(GP[i], GP[9]))~;
Gc = matsolve(GG, -Gb); Gv2 = ip(GP[9], GP[9]) + Gb~ * Gc;
printf("# bd5.2 squared minimum (exact) = %s\n", Gv2);
pr("bd5.2", sqrt(Gv2));
