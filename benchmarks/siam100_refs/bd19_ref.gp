\\ Reference B for holdout-B family D19 (exp 254): a duality certificate in PARI/gp.
\\ For ANY vector V: J(c) = max_W [2 F.W - W' K(c) W] >= 2 F.V - Sum_e c_e E_e(V),
\\ E_e(V) = V_e' K_e V_e, so  min_c J  >=  LB(V) = 2 F.V - max_{c feasible} Sum_e c_e E_e(V),
\\ the inner max a fractional knapsack (exact: sort E_e, fill hi greedily).
\\ UB = J(c*) at the design c* of bd19_ref.py (rounded to 85 digits, made exactly
\\ feasible, assembled and solved here at 200 digits).  UB - LB bounds the error of UB.
\\ Usage: python3.11 bd19_ref.py && gp -q bd19_ref.gp   (same folder: bd19_ref.py writes bd19_<k>.gpdat)
default(realprecision, 200);
cert(k) = {
  read(Str("bd19_", k, ".gpdat"));
  my(c = vector(m, i, round(eval(CS[i]) * 10^85) / 10^85), K = matrix(N, N), U, E, ord, R, cm, ub, lb, s, ix);
  c[FREE[1]] += V - vecsum(c);
  for (e = 1, m, ix = IDX[e]; for (r = 1, #ix, if (ix[r] > 0, for (q = 1, #ix, if (ix[q] > 0,
      K[ix[r], ix[q]] += c[e] * KE[e][r, q])))));
  U = matsolve(1.0 * K, 1.0 * F);
  E = vector(m, e, ix = IDX[e]; my(ue = vector(#ix, r, if (ix[r] > 0, U[ix[r]], 0))); ue * KE[e] * ue~);
  ub = F~ * U;
  ord = vecsort(E, , 4 + 1);   \\ indices, decreasing
  cm = vector(m, i, lo); R = V - m * lo;
  for (i = 1, m, s = min(hi - lo, R); cm[ord[i]] += s; R -= s);
  lb = 2 * F~ * U - sum(e = 1, m, cm[e] * E[e]);
  print("bd19.", k, "  UB = ", precision(ub, 75));
  print("        UB - LB = ", precision(ub - lb, 5), "   min c = ", precision(vecmin(c), 10), "  max c = ", precision(vecmax(c), 10));
}
for (k = 1, 10, iferr(cert(k), err, print("bd19.", k, ": ", err)));
quit;
