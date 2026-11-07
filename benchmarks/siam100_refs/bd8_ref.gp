\\ Holdout B, family D8 -- method B: PARI/gp, realprecision 100.  Cartesian
\\ iterated integrals: the x1-integral of 1/P by residues in t = x1^2
\\ (Int_R dx/B(x^2) = Pi sum_i (-t_i)^(-1/2)/B'(t_i) over the roots of B),
\\ the remaining coordinates by intnum over [0, oo) (P is even in each).
\\ P = sum_{|k|<=m} x^(2k) = sum_a t^a S_{m-a}(others), S_j = sum_{i<=j} h_i.
default(realprecision, 100);   \\ 75 loses ~20 digits on the 3-D cases
inner(B) = my(r = polroots(B), dB = deriv(B)); real(Pi * sum(i = 1, #r, (-r[i])^(-1/2) / subst(dB, 'T, r[i])));
\\ h_i of the squares of the other coordinates
h1(i, y) = y^(2*i);
h2(i, y, z) = sum(j = 0, i, y^(2*j) * z^(2*(i-j)));
B1(m) = sum(a = 0, m, 'T^a);
B2(m, y) = sum(a = 0, m, 'T^a * sum(i = 0, m - a, h1(i, y)));
B3(m, y, z) = sum(a = 0, m, 'T^a * sum(i = 0, m - a, h2(i, y, z)));
c1(m) = sqrt(inner(B1(m)) / (2*Pi));
c2(m) = sqrt(2 * intnum(y = 0, oo, inner(B2(m, y))) / (2*Pi)^2);
\\ N = 3: the Cartesian residue sum cancels catastrophically for large y, z
\\ (tried: ~1e-17 error), so here spherical coordinates over the octant,
\\ t = cos(theta), u = squared direction cosines, angular integrals by intnum:
\\  c_2: radial integral in closed form,
\\       Int_0^oo r^2/(1 + r^2 + q r^4) dr = Pi/(2 sqrt(q) sqrt(1 + 2 sqrt(q))),
\\       q = h_2(u) = (1 + |u|^2)/2;
\\  c_3: radial integral by residues in gp (t = r^2, s = 3/2):
\\       (1/2) Pi/sin(3 Pi/2) sum_i (-t_i)^(1/2)/A'(t_i).
sph(R) = 8 * intnum(t = 0, 1, intnum(p = 0, Pi/2, R([(1 - t^2)*cos(p)^2, (1 - t^2)*sin(p)^2, t^2])));
R2(u) = my(q = (1 + u[1]^2 + u[2]^2 + u[3]^2)/2); Pi/(2*sqrt(q)*sqrt(1 + 2*sqrt(q)));
hh(u, j) = if (j == 0, 1, my(s = 0); forvec(e = vector(3, i, [0, j]), if (vecsum(e) == j, s += u[1]^e[1]*u[2]^e[2]*u[3]^e[3])); s);
R3(u) = my(A = sum(j = 0, 3, hh(u, j)*'T^j), r = polroots(A), dA = deriv(A)); real(-Pi/2 * sum(i = 1, #r, (-r[i])^(1/2)/subst(dA, 'T, r[i])));
c3(m) = sqrt(sph(if (m == 2, R2, R3)) / (2*Pi)^3);
pr(lab, v) = printf("%s %s\n", lab, Strprintf("%.66g", v));
pr("bd8.1", c1(2));
pr("bd8.2", c1(6));
pr("bd8.3", c2(3));
pr("bd8.4", c2(4));
\\ bd8.7 / bd8.10: Int_{R^2} 1/(1 + x^4 + y^4), 1/(1 + x^4 + x^2 y^2 + y^4)
pr("bd8.7", 2 * intnum(y = 0, oo, inner('T^2 + 1 + y^4)));
pr("bd8.10", 2 * intnum(y = 0, oo, inner('T^2 + y^2*'T + 1 + y^4)));
\\ bd8.8 / bd8.9: radial integrals numerically
pr("bd8.8", 4*Pi * intnum(r = 0, oo, r^2/(1 + r^2)^3));
pr("bd8.9", 2*Pi^2 * intnum(r = 0, oo, r^3/(1 + r^2)^3));
pr("bd8.6", c3(2));
pr("bd8.5", c3(3));
