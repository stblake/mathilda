\\ Holdout B, family D3 reference, method B (PARI/gp, 110 digits).
\\ Independent code AND a different coordinate chart from bd3_ref.py: the double
\\ minimiser written by bd3_ref.py (bd3_<p>.start) is mapped to Cartesian
\\ points, rigidly rotated so that point N-1 sits at the north pole and point
\\ N-2 at phi = 0, and Newton is run in THOSE spherical coordinates (analytic
\\ gradient, Hessian by central differences h = 1e-45, matsolve).  LJ13
\\ (bd3.8): Cartesian coordinates with atom 0 at the origin, atom 1 on the z
\\ axis, atom 2 in the xz plane (33 unknowns), same Newton.  Closed forms:
\\ bd3.1 (triangular bipyramid), bd3.2 (icosahedron), bd3.3 (regular 30-gon).
\\ Run: gp -q bd3_ref.gp   (after bd3_ref.py p for p = 4..10)
default(parisizemax, 2000000000);
default(realprecision, 110);
pr(name, v) = printf("%s %.95g\n", name, v);
pr("bd3.1", sqrt(3) + 3 * sqrt(2) + 1/2);
{ my(a = 1 / sin(2 * Pi / 5), ph = (1 + sqrt(5)) / 2); pr("bd3.2", 30 / a + 30 / (ph * a) + 6 / 2); }
pr("bd3.3", 30 / 2 * sum(k = 1, 29, 1 / (2 * sin(Pi * k / 30))));
\\ kernels: kind 1 Coulomb, 2 log, 3 r^-2, 4 -r ;  K(r), K'(r)
Kf(kind, r) = [1 / r, -log(r), 1 / r^2, -r][kind];
dKf(kind, r) = [-1 / r^2, -1 / r, -2 / r^3, -1][kind];
\\ point from (theta, phi)
P(t, f) = [sin(t) * cos(f), sin(t) * sin(f), cos(t)];
\\ chart: point N (1-based) at the pole, point N-1 at phi 0; u = [theta_1..theta_{N-1}, phi_1..phi_{N-2}]
pts(u, N) = {
  my(X = vector(N));
  X[N] = [0, 0, 1];
  for (i = 1, N - 1, X[i] = P(u[i], if (i == N - 1, 0, u[N - 1 + i])));
  X;
}
eg(u, N, kind) = {
  my(X = pts(u, N), G = vector(N, i, [0, 0, 0]), E = 0, d, r, w, g, th, f);
  for (i = 1, N, for (j = i + 1, N,
    d = X[i] - X[j]; r = sqrt(d * d~); E += Kf(kind, r); w = dKf(kind, r) / r;
    G[i] += w * d; G[j] -= w * d));
  g = vector(2 * N - 3);
  for (i = 1, N - 1,
    th = u[i]; f = if (i == N - 1, 0, u[N - 1 + i]);
    g[i] = G[i] * [cos(th) * cos(f), cos(th) * sin(f), -sin(th)]~;
    if (i < N - 1, g[N - 1 + i] = G[i] * [-sin(th) * sin(f), sin(th) * cos(f), 0]~));
  [E, g];
}
\\ rotation taking unit vector a to e_z
rotto(a) = {
  my(z = [0, 0, 1], v, c, s, K);
  v = [a[2] * z[3] - a[3] * z[2], a[3] * z[1] - a[1] * z[3], a[1] * z[2] - a[2] * z[1]];
  c = a * z~; s = sqrt(v * v~);
  if (s < 1e-30, return(matid(3)));
  v = v / s;
  K = [0, -v[3], v[2]; v[3], 0, -v[1]; -v[2], v[1], 0];
  matid(3) + s * K + (1 - c) * K^2;
}
\\ free: coordinates Newton moves (the phi of a point at the south pole of the
\\ chart is a coordinate singularity -- a zero Hessian mode -- and is frozen)
newton(F, u, label, free) = {
  my(n = #free, h = 1e-45, r, H, g, gn, d);
  for (it = 0, 14,
    r = F(u); g = vector(n, i, r[2][free[i]]); gn = sqrt(g * g~);
    printf("# Newton %d |grad| = %.3g  E = %.40g\n", it, gn, r[1]);
    if (gn < 1e-80, break);
    H = matrix(n, n);
    for (k = 1, n, my(up = u, um = u, gp, gm);
      up[free[k]] += h; um[free[k]] -= h; gp = F(up)[2]; gm = F(um)[2];
      for (i = 1, n, H[i, k] = (gp[free[i]] - gm[free[i]]) / (2 * h)));
    d = matsolve(H, g~);
    for (i = 1, n, u[free[i]] -= d[i]));
  printf("# Hessian min eigenvalue %.4g\n", vecmin(real(qfjacobi((H + H~) / 2)[1])));
  pr(label, r[1]);
}
sphere(p, N, kind) = {
  my(v = readvec(Str("bd3_", p, ".start")), ph1 = Pi * (3 - sqrt(5)), th, f, X, R, R2, Y, a, u);
  th = concat([0], vector(N - 1, i, v[i])); f = concat([0, ph1], vector(N - 2, i, v[N - 1 + i]));
  X = vector(N, i, P(th[i], f[i]));
  R = rotto(X[N]); Y = vector(N, i, (R * X[i]~)~);
  a = atan2(Y[N - 1][2], Y[N - 1][1]);
  R2 = [cos(a), sin(a), 0; -sin(a), cos(a), 0; 0, 0, 1];
  Y = vector(N, i, (R2 * Y[i]~)~);
  u = vector(2 * N - 3);
  for (i = 1, N - 1, u[i] = acos(Y[i][3]); if (i < N - 1, u[N - 1 + i] = atan2(Y[i][2], Y[i][1])));
  my(free = []);
  for (i = 1, 2 * N - 3, if (i < N || abs(sin(u[i - N + 1])) > 1e-4, free = concat(free, [i])));
  printf("# bd3.%d: %d of %d coordinates free\n", p, #free, 2 * N - 3);
  newton(w -> eg(w, N, kind), u, Str("bd3.", p), free);
}
atan2(y, x) = arg(x + I * y);
sphere(4, 24, 1);
sphere(5, 32, 1);
sphere(6, 20, 2);
sphere(7, 20, 3);
sphere(10, 30, 4);
sphere(9, 60, 1);
\\ LJ13
ljx(u) = { my(X = vector(13)); X[1] = [0, 0, 0]; X[2] = [0, 0, u[1]]; X[3] = [u[2], 0, u[3]];
  for (i = 4, 13, X[i] = [u[3 * i - 8], u[3 * i - 7], u[3 * i - 6]]); X; }
ljeg(u) = {
  my(X = ljx(u), G = vector(13, i, [0, 0, 0]), E = 0, d, r2, i6, w, g);
  for (i = 1, 13, for (j = i + 1, 13,
    d = X[i] - X[j]; r2 = d * d~; i6 = 1 / r2^3; E += 4 * (i6^2 - i6);
    w = 4 * (-12 * i6^2 + 6 * i6) / r2; G[i] += w * d; G[j] -= w * d));
  g = vector(33); g[1] = G[2][3]; g[2] = G[3][1]; g[3] = G[3][3];
  for (i = 4, 13, for (c = 1, 3, g[3 * i - 9 + c] = G[i][c]));
  [E, g];
}
{ my(v = readvec("bd3_8.start"), X, R, Y, a, R2, u);
  X = vector(13, i, [v[3 * i - 2], v[3 * i - 1], v[3 * i]]);
  X = vector(13, i, X[i] - X[1]);
  R = rotto(X[2] / sqrt(X[2] * X[2]~)); Y = vector(13, i, (R * X[i]~)~);
  a = atan2(Y[3][2], Y[3][1]); R2 = [cos(a), sin(a), 0; -sin(a), cos(a), 0; 0, 0, 1];
  Y = vector(13, i, (R2 * Y[i]~)~);
  u = concat([Y[2][3], Y[3][1], Y[3][3]], concat(vector(10, i, Y[i + 3])));
  newton(ljeg, u, "bd3.8", vector(33, i, i)); }
