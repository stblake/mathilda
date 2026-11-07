\\ Holdout B, family D4 reference, method B (PARI/gp, 110 digits): independent
\\ formulas from bd4_ref.py (mpmath).   Run: gp -q bd4_ref.gp
\\  bd4.1  li(x) = gamma + log log x + sum (log x)^k/(k k!)  (series, not eint1)
\\  bd4.2  zeta(s) = 2                (PARI zeta)
\\  bd4.3  sum n^-s/log n = int_s^oo (zeta(t) - 1) dt = 1   (intnum; method A sums the series)
\\  bd4.4, bd4.9  R(e^-u) by the zeta-zero (Mellin-Barnes residue) expansion
\\         sum_n (-1)^(n+1) 2 (2pi)^(2n)/((2n+1) zeta(2n+1) u^(2n+1))
\\           - 2 Re sum_rho Gamma(1-rho) u^(rho-1)/((1-rho) zeta'(rho))
\\         (zeros from lfunzeros, zeta' from lfun(.,.,1)); root in t = -u
\\  bd4.5  R(x) by the Gram series 1 + sum (log x)^k/(k k! zeta(k+1))
\\  bd4.6  sum_{k<=4000} moebius(k) x^k = 0      (direct, PARI moebius)
\\  bd4.7  first zero from lfunzeros(1, [10, 20])
\\  bd4.8  P(s) = sum_k mu(k)/k log zeta(k s) = 1  (Moebius inversion; method A: mpmath.primezeta)
\\  bd4.10 zeta'(s) = lfun(1, s, 1) = -1
default(parisizemax, 2000000000);
default(realprecision, 110);
pr(name, v) = printf("%s %.95g\n", name, v);
li(x) = { my(y = log(x)); Euler + log(y) + suminf(k = 1, y^k / (k * k!)); }
pr("bd4.1", solve(x = 1.3, 1.6, li(x)));
pr("bd4.2", solve(s = 1.5, 2, zeta(s) - 2));
I3(s) = intnum(t = s, [oo, log(2)], zeta(t) - 1);
pr("bd4.3", solve(s = 1.2, 2.5, I3(s) - 1));
\\ zeta-zero expansion
Z = lfunzeros(1, [0, 260]);
ZP = vector(#Z, j, my(r = 1/2 + I * Z[j]); [r, lfun(1, r, 1)]);
Rasym(u) = {
  my(s = 0., n = 1, t);
  while (1, t = (-1)^(n + 1) * 2 * (2 * Pi)^(2 * n) / ((2 * n + 1) * zeta(2 * n + 1) * u^(2 * n + 1));
         s += t; if (abs(t) < 10^-130, break); n++);
  for (j = 1, #ZP, my(r = ZP[j][1]); s -= 2 * real(gamma(1 - r) * u^(r - 1) / ((1 - r) * ZP[j][2])));
  s;
}
printf("# zeros used: %d, last |Gamma(1-rho)/zeta'(rho)| = %.3g\n", #Z, abs(gamma(1 - ZP[#Z][1]) / ZP[#Z][2]));
pr("bd4.4", -solve(u = 20, 40, Rasym(u) - 1/1000));
gram(t) = 1 + suminf(k = 1, t^k / (k * k! * zeta(k + 1)));
pr("bd4.5", solve(t = -12, -3, gram(t) - 1/10));
MU = vector(4000, k, moebius(k));
f6(x) = sum(k = 1, 4000, MU[k] * x^k);
pr("bd4.6", solve(x = 0.5, 0.7, f6(x)));
pr("bd4.7", lfunzeros(1, [10, 20])[1]);
P(s) = sum(k = 1, 400, moebius(k) / k * log(zeta(k * s)));
pr("bd4.8", solve(s = 1.1, 2, P(s) - 1));
pr("bd4.9", -solve(u = 2000, 4000, Rasym(u) - 10^-9));
pr("bd4.10", solve(s = 1.5, 3, lfun(1, s, 1) + 1));
