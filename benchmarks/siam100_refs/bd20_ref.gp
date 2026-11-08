\\ Reference B for holdout-B family D20 (exp 255), PARI/gp at 220 digits.
\\ Robbins rows: matsolve of the n x n system, root of slog_n(x) = level by
\\ solve() on a bracket.  Carleman rows bd20.7 / bd20.9: principal square root
\\ through the eigendecomposition C = P D P^-1 (distinct positive eigenvalues
\\ 2^k), first row of P D^(1/2) P^-1.  (bd20.8 is unipotent, not
\\ diagonalisable: its method B is mpmath.sqrtm in bd20_ref.py.)
\\ Usage: gp -q bd20_ref.gp
default(realprecision, 220);
robbins(n, L) = {
  my(A = matrix(n, n, j, k, (k*L)^(j-1)/(j-1)! - (j-1 == k)), b = vector(n, j, j == 1)~);
  matsolve(A, b);
}
slog(v, x) = -1 + sum(k = 1, #v, v[k]*x^k);
root(v, lev, a, b) = solve(x = a, b, slog(v, x) - lev);
pr(lab, val) = print(lab, "  ", precision(val, 70));
pr("bd20.1", root(robbins(10, 1), -1/2, 0.4, 0.6));
pr("bd20.2", root(robbins(20, log(2)), -1/2, 0.45, 0.65));
pr("bd20.3", root(robbins(30, log(2)/2), -1/2, 0.55, 0.7));
pr("bd20.4", root(robbins(15, log(10)), -1/2, 0.3, 0.5));
pr("bd20.5", root(robbins(30, 1), -3/4, 0.2, 0.35));
pr("bd20.6", robbins(24, 1)[1]);
pr("bd20.10", root(robbins(45, 1), -1/2, 0.4, 0.6));
sqrtrow(C) = {
  my(P = mateigen(C), D, S);
  D = P^-1 * C * P;
  S = P * matdiagonal(vector(#C, i, sqrt(D[i, i]))) * P^-1;
  S[1, ];
}
n = 20; C7 = matrix(n, n, j, k, if(k >= j && k <= 2*j, binomial(j, k-j)*2^(2*j-k), 0));
r = sqrtrow(C7); pr("bd20.7", sum(k = 1, n, r[k]*(1/2)^k));
n = 25; C9 = matrix(n, n, j, k, if(k >= j && (k-j)%2 == 0, 2^j*binomial(j + (k-j)/2 - 1, (k-j)/2), 0));
r = sqrtrow(C9); pr("bd20.9", sum(k = 1, n, r[k]*(3/10)^k));
quit;
