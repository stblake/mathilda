\\ Holdout B, family D6 -- method B: PARI/gp intnum (double-exponential) at
\\ realprecision 80, split at the points where the phase phi(x) = k*pi + pi/2
\\ (a different split from bd6_ref.py), phi^{-1} by solve() on [a,b].
default(realprecision, 80);
sp(f, a, b, phi) = {
  my(lo = phi(a), hi = phi(b), pts = List([a]), s = 0.);
  for (k = ceil(lo/Pi - 1/2), floor(hi/Pi - 1/2),
    my(y = k*Pi + Pi/2);
    if (y > lo && y < hi, listput(pts, solve(x = a, b, phi(x) - y))));
  listput(pts, b);
  for (i = 1, #pts - 1, s += intnum(x = pts[i], pts[i+1], f(x)));
  s;
}
pr(lab, v) = print(lab, " ", v);
pr("bd6.1",  sp(x -> exp(sin(exp(x))), 0, 1, x -> exp(x)));
pr("bd6.2",  sp(x -> cos(x + exp(exp(x + 1/2))), -1, 1, x -> exp(exp(x + 1/2))));
pr("bd6.3",  sp(x -> exp(cos(x^3 + exp(x))), 0, 2, x -> x^3 + exp(x)));
pr("bd6.4",  sp(x -> log(2 + sin(exp(exp(x)))), -1, 1, x -> exp(exp(x))));
\\ near-poles at exp(4x) = (2j+1) pi: split 8x finer (phase scaled by 8)
pr("bd6.5",  sp(x -> 1/(5/4 + cos(exp(4*x))), 0, 1, x -> 8*exp(4*x)));
pr("bd6.6",  sp(x -> sin(exp(exp(2*x))), 0, 1, x -> exp(exp(2*x))));
pr("bd6.7",  sp(x -> cos(50*sin(x) - 3*x), 0, Pi, x -> 50*sin(x/2)));
pr("bd6.7 closed", Pi*besselj(3, 50));
pr("bd6.8",  sp(x -> exp(x)*cos(sinh(exp(x + 1/2))), -1, 1, x -> sinh(exp(x + 1/2))));
pr("bd6.9",  sp(x -> cos(exp(exp(exp(x)))), -1, 1/2, x -> exp(exp(exp(x)))));
pr("bd6.10", sp(x -> exp(x + sin(exp(exp(x + 1/2)))), -1, 1, x -> exp(exp(x + 1/2))));
