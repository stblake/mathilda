\\ Holdout B, family C2 (experiment 227) -- reference values, route B (PARI/gp).
\\ An independent implementation of the same billiards: mirrors as
\\ [kind, a, b, c, id] (kind 0: disk centre (a,b) radius c; 1: wall a x + b y = c),
\\ the next hit as the smallest positive root among ALL candidate mirrors
\\ (lattice disks: every disk in the box around the remaining segment, no
\\ line-distance prefilter -- polroots of the hit polynomial, real roots only),
\\ reflection v' = v - 2 (v.n) n.  realprecision 150 (and 200 for the check).
s3() = sqrt(3);
box2(x0, x1, y0, y1) = [[1,1,0,x1,-1],[1,-1,0,-x0,-2],[1,0,1,y1,-3],[1,0,-1,-y0,-4]];
sinai() = concat(box2(-1,1,-1,1), [[0,0,0,1/2,-5]]);
\\ [lattice kind (0 none, 1 square, 2 triangular), r, mirrors, p0, v0, T (0: none), reflections (0: none)]
sc(k) = {
  if (k == 1, [1, 1/3, [], [1/2, 1/10], [3/5, 4/5], 6, 0],
  k == 2, [1, 1/4, [], [1/2, 1/2], [4/5, 3/5], 5, 0],
  k == 3, [2, 1/3, [], [1/2, 1/5], [1, 0], 4, 0],
  k == 4, [0, 0, sinai(), [3/4, 0], [3/5, 4/5], 12, 0],
  k == 5, [0, 0, concat(box2(0,2,0,1), [[0,1/2,1/2,1/4,-5],[0,3/2,1/2,1/4,-6]]), [1, 1/10], [3/5, 4/5], 6, 0],
  k == 6, [0, 0, [[0,0,6/5,1,-1],[0,-3*s3()/5,-3/5,1,-2],[0,3*s3()/5,-3/5,1,-3]], [0, 0], [12/13, -5/13], 20, 0],
  k == 7, [0, 0, sinai(), [-1/2, -3/4], [5/13, 12/13], 0, 16],
  k == 8, [0, 0, box2(0,3,0,2), [1/3, 1/7], [12/13, 5/13], 100, 0],
  k == 9, [1, 9/20, [], [1/2, 1/2], [3/5, 4/5], 3, 0],
  k == 10, [0, 0, sinai(), [3/4, 0], [3/5, 4/5], 25, 0]);
}
latd(kind, r, p, q) = {
  my(L = List(), hh = s3()/2, x0 = min(p[1], q[1]), x1 = max(p[1], q[1]), y0 = min(p[2], q[2]), y1 = max(p[2], q[2]));
  if (kind == 1, for (i = floor(x0) - 1, floor(x1) + 2, for (j = floor(y0) - 1, floor(y1) + 2,
      listput(L, [0, i, j, r, 10^6 * (i + 5000) + j + 5000]))));
  if (kind == 2, for (j = floor(y0/hh) - 1, floor(y1/hh) + 2, for (i = floor(x0 - j/2) - 1, floor(x1 - j/2) + 2,
      listput(L, [0, i + j/2, j*hh, r, 10^6 * (i + 5000) + j + 5000]))));
  Vec(L);
}
posroots(P) = my(R = polroots(P), out = List()); for (i = 1, #R, if (abs(imag(R[i])) < 10^(-default(realprecision)/2) && real(R[i]) > 0, listput(out, real(R[i])))); Vec(out);
fly(k) = {
  my(S = sc(k), p = S[4] * 1., v = S[5] * 1., hasT = S[6] > 0, left = S[6] * 1., nmax = S[7],
     last = 0, clock = 0., nb = 0, C, bt, bm, P, R, d, n);
  while (1,
    C = concat(S[3], latd(S[1], S[2], p, p + min(if (hasT, left, 50), 50) * v));
    bt = -1; bm = 0;
    for (q = 1, #C, my(m = C[q]);
      if (m[5] == last, next);
      if (m[1] == 0,
        P = (p[1] + 'x*v[1] - m[2])^2 + (p[2] + 'x*v[2] - m[3])^2 - m[4]^2,
        P = m[2] * (p[1] + 'x*v[1]) + m[3] * (p[2] + 'x*v[2]) - m[4]);
      if (poldegree(P) < 1, next);
      R = posroots(P);
      \\ a disk is hit at the SMALLER root (entering); the larger is the exit
      if (#R && m[1] == 0, R = if (#R == 2, [vecmin(R)], if (#R == 1 && poldegree(P) == 2 && #select(t -> real(t) < 0, polroots(P)) == 0, [], R)));
      for (i = 1, #R, if (bt < 0 || R[i] < bt, bt = R[i]; bm = m)));
    if (bt < 0 || (hasT && bt >= left), p += left * v; return (sqrt(p[1]^2 + p[2]^2)));
    p += bt * v; clock += bt; nb++; last = bm[5];
    if (hasT, left -= bt);
    n = if (bm[1] == 0, [p[1] - bm[2], p[2] - bm[3]] / bm[4], [bm[2], bm[3]]);
    d = v[1]*n[1] + v[2]*n[2]; v -= 2 * d * n;
    if (nmax && nb == nmax, return (clock)));
}
\p 150
r150 = vector(10, k, fly(k));
\p 200
r200 = vector(10, k, fly(k));
for (k = 1, 10, printf("bc2.%d %.80g  (150 vs 200 agree %.1f)\n", k, r200[k], -log(abs(r150[k] - r200[k]) / abs(r200[k]) + 10^-190) / log(10)));
quit
