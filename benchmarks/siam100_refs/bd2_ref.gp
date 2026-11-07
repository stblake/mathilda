\\ Holdout B, family D2 reference, method B (PARI/gp, 110 digits): independent
\\ code and decomposition from bd2_ref.py.  Every series is summed DIRECTLY
\\ over a head n < b^K (K = 5..20), and only the tail is grouped by digit
\\ length (block sums from psi); the log b part of the tail closes with a
\\ Hurwitz zeta, an alternating-harmonic remainder, or the self-similar
\\ linear equation for S solved from the head rather than from n = 1.
\\ bd2.3/bd2.4: head to N = 10^4, then Euler-Maclaurin with the EXACT tail
\\ integral (1/log N, resp. 1/loglog N) and exact Taylor derivatives at N
\\ (power-series arithmetic, 40 Bernoulli terms).   Run: gp -q bd2_ref.gp
default(parisizemax, 2000000000);
default(realprecision, 110);
nd(n, b) = #digits(n, b);
Hb(b, k) = psi(1. * b^k) - psi(1. * b^(k - 1));
pr(name, v) = printf("%s %.95g\n", name, v);
\\ bd2.1
{ my(S = sum(n = 1, 10^5 - 1, 1. / (n * nd(n, 10)^2)));
  S += sum(k = 6, 150, Hb(10, k) / k^2) + log(10) * zetahurwitz(2, 151);
  pr("bd2.1", S); s1 = S; }
\\ bd2.2
{ my(S = sum(n = 1, 2^20 - 1, 1. / (n * nd(n, 2)^2)));
  S += sum(k = 21, 400, Hb(2, k) / k^2) + log(2) * zetahurwitz(2, 401);
  pr("bd2.2", S); }
\\ bd2.3, bd2.4: Euler-Maclaurin from N with exact derivatives
em(f, a, N, I) = {
  my(S = sum(n = a, N - 1, f(1. * n)), t = f(N + x + O(x^90)));
  S += I + polcoef(t, 0) / 2;
  for (k = 1, 40, S -= bernfrac(2 * k) / (2 * k) * polcoef(t, 2 * k - 1));
  S;
}
pr("bd2.3", em(z -> 1 / (z * log(z)^2), 2, 10^4, 1 / log(10.^4)));
pr("bd2.4", em(z -> 1 / (z * log(z) * log(log(z))^2), 3, 10^4, 1 / log(log(10.^4))));
\\ bd2.5: g(1) = 1, g(n) = n g(floor(log2 n)) = n g(bitlength(n) - 1)
g5(n) = if (n == 1, 1, n * g5(nd(n, 2) - 1));
{ my(K = 16, head = sum(n = 1, 2^K - 1, 1. / g5(n)), h = sum(m = 1, K - 1, 1. / g5(m)),
     E = sum(m = K, 400, (Hb(2, m + 1) - log(2)) / g5(m)));
  \\ blocks k >= K+1 contribute B_k / g(k-1):  S = head + log2 (S - h) + E
  pr("bd2.5", (head - log(2) * h + E) / (1 - log(2))); }
\\ bd2.6: f(1) = 1, f(n) = n f(d(n))^2
f6(n) = if (n == 1, 1, n * f6(nd(n, 10))^2);
{ my(S = sum(n = 1, 10^5 - 1, 1. / f6(n)));
  S += sum(k = 6, 999, Hb(10, k) / f6(k)^2);
  \\ k >= 1000: H_k = log 10 to 1e-1000; group k by its digit count m >= 4
  S += log(10) * sum(m = 4, 200, (psi(1. * 10^(m - 1), 1) - psi(1. * 10^m, 1)) / f6(m)^4);
  pr("bd2.6", S); }
\\ bd2.7: sum 1/(n d(n) d(d(n))^2)
{ my(S = sum(n = 1, 10^5 - 1, 1. / (n * nd(n, 10) * nd(nd(n, 10), 10)^2)));
  S += sum(k = 6, 150, Hb(10, k) / (k * nd(k, 10)^2));
  S += log(10) * sum(k = 151, 9999, 1. / (k * nd(k, 10)^2));
  \\ k >= 10^4: group k by digit count m >= 5
  S += log(10) * (sum(m = 5, 150, Hb(10, m) / m^2) + log(10) * zetahurwitz(2, 151));
  pr("bd2.7", S); }
\\ bd2.8: sum 1/(n d(n)!)
{ my(S = sum(n = 1, 10^6 - 1, 1. / (n * nd(n, 10)!)));
  S += sum(k = 7, 150, Hb(10, k) / k!);
  pr("bd2.8", S); }
\\ bd2.9: sum (-1)^b(n) / (n b(n))
{ my(S = sum(n = 1, 2^20 - 1, (-1.)^nd(n, 2) / (n * nd(n, 2))));
  S += sum(k = 21, 400, (-1)^k * (Hb(2, k) - log(2)) / k);
  S += log(2) * (-log(2) - sum(k = 1, 20, (-1.)^k / k));
  pr("bd2.9", S); }
\\ bd2.10: f(1) = 1, f(2) = 2, f(n) = 2 n f(t(n)), t = ternary length
f10(n) = if (n <= 2, n, 2 * n * f10(nd(n, 3)));
{ my(K = 10, head = sum(n = 1, 3^K - 1, 1. / f10(n)), h = sum(k = 1, K, 1. / f10(k)),
     E = sum(k = K + 1, 300, (Hb(3, k) - log(3)) / (2 * f10(k))));
  \\ blocks k >= K+1 contribute T_k / (2 f(k)):  S = head + (log3/2)(S - h) + E
  pr("bd2.10", (head - log(3) / 2 * h + E) / (1 - log(3) / 2)); }
