\\ Holdout B, family D9 -- method B: PARI/gp, graphs built here independently
\\ of bd9_common.py; exact rational propagation of the start distribution.
\\ Each walk: nbr(v) = list of neighbour vertices (uniform choice).
prop(N, nbr, s, tg, n) = {
  my(v = vector(N), w, l, x);
  v[s] = 1;
  for (k = 1, n, w = vector(N);
    for (i = 1, N, if (v[i], l = nbr(i); for (t = 1, #l, w[l[t]] += v[i]/#l)));
    v = w);
  x = sum(t = 1, #tg, v[tg[t]]);
  x;
}
pr(lab, x) = printf("%s %s\n", lab, Strprintf("%.66g", x * 1.));
default(realprecision, 80);
\\ grids: vertex = mixed-radix index of coordinates (row-major, last fastest)
gidx(p, d) = my(i = 0); for (t = 1, #d, i = i*d[t] + p[t]); i + 1;
gco(i, d) = my(p = vector(#d), j = i - 1); forstep (t = #d, 1, -1, p[t] = j % d[t]; j \= d[t]); p;
gnbr(i, d, mv) = {
  my(p = gco(i, d), L = List(), q, ok);
  for (m = 1, #mv, q = p + mv[m]; ok = 1;
    for (t = 1, #d, if (q[t] < 0 || q[t] >= d[t], ok = 0));
    if (ok, listput(L, gidx(q, d))));
  Vec(L);
}
\\ triangle with n per side: (r, c), 0 <= c <= r < n, index r(r+1)/2 + c + 1
tco(i) = my(r = 0); while ((r+1)*(r+2)/2 < i, r++); [r, i - 1 - r*(r+1)/2];
tnbr(i, n) = {
  my(p = tco(i), r = p[1], c = p[2], L = List(), mv = [[0,-1],[0,1],[-1,-1],[-1,0],[1,0],[1,1]], a, b);
  for (m = 1, 6, a = r + mv[m][1]; b = c + mv[m][2];
    if (0 <= b && b <= a && a < n, listput(L, a*(a+1)/2 + b + 1)));
  Vec(L);
}
D1 = [21];          pr("bd9.1", prop(21, i -> gnbr(i, D1, [[1], [-1]]), 1, vector(10, k, k + 11), 40));
D2 = [10, 10];      pr("bd9.2", prop(100, i -> gnbr(i, D2, [[1,0],[-1,0],[0,1],[0,-1]]), 1, vector(10, j, gidx([9, j-1], D2)), 30));
pr("bd9.3", prop(210, i -> tnbr(i, 20), 1, vector(20, c, 19*20/2 + c), 40));
D4 = [6, 6, 6];     pr("bd9.4", prop(216, i -> gnbr(i, D4, [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]), 1, vector(36, k, gidx([5, (k-1)\6, (k-1)%6], D4)), 24));
D5 = [8, 8];        pr("bd9.5", prop(64, i -> gnbr(i, D5, [[1,2],[2,1],[-1,2],[-2,1],[1,-2],[2,-1],[-1,-2],[-2,-1]]), 1, [gidx([7,7], D5)], 20));
KM = [[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]];
pr("bd9.6", prop(64, i -> gnbr(i, D5, KM), 1, vector(8, k, gidx([k-1, 7], D5)), 30));
pr("bd9.7", prop(1024, i -> vector(10, t, bitxor(i - 1, 2^(t-1)) + 1), 1, [1024], 40));
pr("bd9.8", prop(17, i -> [(i % 17) + 1, ((i - 2) % 17) + 1], 1, [1], 60));
\\ 9: gambler's ruin -- not uniform: own propagation, states 0..30 -> 1..31
{ my(v = vector(31), w); v[16] = 1;
  for (k = 1, 200, w = vector(31); w[1] += v[1]; w[31] += v[31];
    for (i = 2, 30, w[i+1] += 2/5*v[i]; w[i-1] += 3/5*v[i]); v = w);
  pr("bd9.9", v[31]); }
pr("bd9.10", prop(820, i -> tnbr(i, 40), 1, vector(40, c, 39*40/2 + c), 100));
