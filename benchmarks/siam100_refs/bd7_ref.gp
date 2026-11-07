\\ Holdout B, family D7 -- method B: PARI/gp, realprecision 90.  Roots by
\\ solve() on sign-change brackets (not Newton from a start), integrals by
\\ intnum (double exponential), breakpoints at roots and kinks; cube roots
\\ real: cb(t) = sign(t) |t|^(1/3).
default(realprecision, 90);
cb(t) = t = real(t); if (t == 0, 0, sign(t) * abs(t)^(1/3));
pr(lab, v) = print(lab, " ", real(v));
x0 = solve(x = 0.5, 1, cos(x) - x^2);
pr("bd7.1", intnum(x = -x0, x0, cos(x) - x^2));
a = solve(x = -2, -1.5, exp(x) + x^2 - 3); b = solve(x = 0.5, 1, exp(x) + x^2 - 3);
pr("bd7.2", intnum(x = a, b, 3 - x^2 - exp(x)));
pr("bd7.3", intnum(x = Pi/6, Pi/2, Pi - 2*asin(1/(2*sin(x)))) + intnum(x = Pi/2, 5*Pi/6, Pi - 2*asin(1/(2*sin(x)))));
g4(y) = 1/2 + y^2/8;
s4(y) = sqrt(-log(g4(y)));
yA = solve(y = -1, -0.5, -log(g4(y)) - y^2); yB = solve(y = 0.5, 1, -log(g4(y)) - y^2);
pr("bd7.4", intnum(y = yA, 0, 2*sqrt(y + s4(y))) + intnum(y = 0, yB, 2*sqrt(y + s4(y))) + intnum(y = yB, 2, 2*(sqrt(y + s4(y)) - sqrt(y - s4(y)))));
r1 = solve(x = -2, -1.5, x^3 - 3*x + 1); r2 = solve(x = 0, 1, x^3 - 3*x + 1);
pr("bd7.5", intnum(x = r1, 0, 2*sqrt(x^3 - 3*x + 1)) + intnum(x = 0, r2, 2*sqrt(x^3 - 3*x + 1)));
pr("bd7.6", 2*intnum(x = 0, x0, 2*(cos(x) - x^2)^(1/4)));
h7(u) = 1 + sin(sqrt(2)*u) - u^2;
u1 = solve(u = -1, 0, h7(u)); u2 = solve(u = 1, 2, h7(u));
pr("bd7.7", intnum(u = u1, 0, 2*sqrt(h7(u))) + intnum(u = 0, u2, 2*sqrt(h7(u))));
g8(y) = y^2/16 + exp(cos(y) - 1);
s8(y) = sqrt(-log(g8(y)));
w8(y) = cb(s8(y) - 2*y) - cb(-s8(y) - 2*y);
p8 = solve(y = 3, 4, g8(y) - 1); m8 = solve(y = -4, -3, g8(y) - 1);
pr("bd7.8", intnum(y = m8, 0, w8(y)) + intnum(y = 0, p8, w8(y)));
c9(x) = sqrt(9/4 - x^2);
a9(x) = acos(1 - cos(x));
x1 = solve(x = 0.6, 1, c9(x) - a9(x)); x2 = solve(x = 1.1, 1.45, c9(x) - a9(x));
pr("bd7.9", 4*(intnum(x = 0, x1, c9(x)) + intnum(x = x1, x2, a9(x)) + intnum(x = x2, 3/2, c9(x))));
g10(y) = y^2/100 + exp(sin(y));
s10(y) = sqrt(-log(g10(y)));
w10(y) = cb(-y + s10(y)) - cb(-y - s10(y));
G(y) = g10(y) - 1;
y1 = solve(y = -8.2, -7.5, G(y)); y2 = solve(y = -7.5, -6.5, G(y)); y3 = solve(y = -3.5, -2.5, G(y));
yc = solve(y = -1, -0.7, -log(g10(y)) - y^2); y4 = solve(y = 3, 4, G(y)); y5 = solve(y = 5.5, 6.2, G(y));
pr("bd7.10", intnum(y = y1, y2, w10(y)) + intnum(y = y3, yc, w10(y)) + intnum(y = yc, 0, w10(y)) + intnum(y = y4, y5, w10(y)));
