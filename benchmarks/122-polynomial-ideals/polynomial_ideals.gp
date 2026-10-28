\\ Experiment 122 -- polynomial gcds, lcms and ideals (PARI/GP column).
\\ See polynomial_ideals.m for the cases and the checks.
\\
\\ PARI has native multivariate gcd and lcm, over Q and over F_p (coefficients
\\ times Mod(1, p)).  It has no Groebner bases, no multivariate division by a
\\ list, no elimination beyond the two-polynomial resultant and no cylindrical
\\ decomposition, so those cases are SKIP.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");

require(["gcd", "lcm", "poldegree", "substvec"]);

P = 32003;
pt1(n) = vector(n, i, (2*i + 1)/(3*i + 5));
pt2(n) = vector(n, i, (3*i + 2)/(5*i + 7));
q1(n) = vector(n, i, i + 2);
q2(n) = vector(n, i, 2*i + 5);
ratio(g, vs) = substvec(g, vs, pt1(#vs)) / substvec(g, vs, pt2(#vs));
ratioMod(g, vs) = lift(Mod(substvec(lift(g), vs, q1(#vs)), P) / Mod(substvec(lift(g), vs, q2(#vs)), P));
polyChk(g, vs) = 10^40 * poldegree(g, vs[1]) + round(10^6 * ratio(g, vs));
polyChkMod(g, vs) = 10^12 * poldegree(lift(g), vs[1]) + ratioMod(g, vs);

denseUPoly(n) = sum(k = 0, n, ((k^2 + 3*k + 1) % 17) * x^k) + 1;
factorableUPoly(k) = prod(i = 1, k, x^2 + i*x + (i + 1));

xs = [x]; v3 = [x, y, z]; v6 = [x, y, z, u, v, w];

\\ ---- PolynomialGCD -----------------------------------------------------------

sh100 = factorableUPoly(100);
g400a = sh100 * denseUPoly(200);
g400b = sh100 * (denseUPoly(199) + 1);
benchLoop("PolynomialGCD, deg-400, deg-200 gcd", () -> gcd(g400a, g400b));
check("PolynomialGCD, deg-400, deg-200 gcd", polyChk(gcd(g400a, g400b), xs));

sh250 = factorableUPoly(250);
g1000a = sh250 * denseUPoly(500);
g1000b = sh250 * (denseUPoly(499) + 1);
benchLoop("PolynomialGCD, deg-1000, deg-500 gcd", () -> gcd(g1000a, g1000b));
check("PolynomialGCD, deg-1000, deg-500 gcd", polyChk(gcd(g1000a, g1000b), xs));

cp1000a = denseUPoly(1000);
cp1000b = denseUPoly(999) + x^500 + 2;
benchLoop("PolynomialGCD, deg-1000 coprime", () -> gcd(cp1000a, cp1000b));
check("PolynomialGCD, deg-1000 coprime", polyChk(gcd(cp1000a, cp1000b), xs));

shBig = (10^50 + 3)*x^100 + denseUPoly(99);
big200a = shBig * (denseUPoly(100) + 7^60);
big200b = shBig * (denseUPoly(101) - 11^50);
benchLoop("PolynomialGCD, deg-200, 60-digit coefficients", () -> gcd(big200a, big200b));
check("PolynomialGCD, deg-200, 60-digit coefficients", polyChk(gcd(big200a, big200b), xs));

g2 = (x^2 + 3*x*y + y^3 + 1)^5;
bi40a = g2 * (x - y^2 + 3)^8;
bi40b = g2 * (x^3 + y + 5)^6;
benchLoop("PolynomialGCD, dense bivariate, deg-15 gcd", () -> gcd(bi40a, bi40b));
check("PolynomialGCD, dense bivariate, deg-15 gcd", polyChk(gcd(bi40a, bi40b), v3));

g3 = (1 + x + y + z)^6;
tri6a = g3 * (x + 2*y - z + 3)^4;
tri6b = g3 * (x - y + 2*z + 1)^4;
benchLoop("PolynomialGCD, dense trivariate, deg-6 gcd", () -> gcd(tri6a, tri6b));
check("PolynomialGCD, dense trivariate, deg-6 gcd", polyChk(gcd(tri6a, tri6b), v3));

g6 = 1 + x*y^3 + z^2*u*v^4 + w^5 + x*w;
sp6a = g6^2 * (x^2 + z*v + 2) * (y*u + w^2 + 1);
sp6b = g6^2 * (x*z + u^3 + 3) * (v*w + y^2 + 1);
benchLoop("PolynomialGCD, sparse 6-variable", () -> gcd(sp6a, sp6b));
check("PolynomialGCD, sparse 6-variable", polyChk(gcd(sp6a, sp6b), v6));

ps10 = vector(10, i, (x^2 + y^2 + x*y + 1)^2 * (x + y + i)^3 * (x - i*y + 1)^2);
benchLoop("PolynomialGCD, 10 bivariate arguments", () -> gcd(ps10));
check("PolynomialGCD, 10 bivariate arguments", polyChk(gcd(ps10), v3));

\ The reduction mod p is part of the operation (the other systems take the
\ integer polynomials and reduce inside the call), so it is inside the timing.
modGcd(a, b) = gcd(a * Mod(1, P), b * Mod(1, P));
modLcm(a, b) = lcm(a * Mod(1, P), b * Mod(1, P));
benchLoop("PolynomialGCD, Modulus 32003, deg-1000", () -> modGcd(g1000a, g1000b));
check("PolynomialGCD, Modulus 32003, deg-1000", polyChkMod(modGcd(g1000a, g1000b), xs));

benchLoop("PolynomialGCD, Modulus 32003, trivariate", () -> modGcd(tri6a, tri6b));
check("PolynomialGCD, Modulus 32003, trivariate", polyChkMod(modGcd(tri6a, tri6b), v3));

\\ ---- PolynomialLCM -----------------------------------------------------------

benchLoop("PolynomialLCM, deg-400, deg-200 gcd", () -> lcm(g400a, g400b));
check("PolynomialLCM, deg-400, deg-200 gcd", polyChk(lcm(g400a, g400b), xs));

benchLoop("PolynomialLCM, deg-1000, deg-500 gcd", () -> lcm(g1000a, g1000b));
check("PolynomialLCM, deg-1000, deg-500 gcd", polyChk(lcm(g1000a, g1000b), xs));

benchLoop("PolynomialLCM, dense trivariate, deg-6 gcd", () -> lcm(tri6a, tri6b));
check("PolynomialLCM, dense trivariate, deg-6 gcd", polyChk(lcm(tri6a, tri6b), v3));

benchLoop("PolynomialLCM, sparse 6-variable", () -> lcm(sp6a, sp6b));
check("PolynomialLCM, sparse 6-variable", polyChk(lcm(sp6a, sp6b), v6));

benchLoop("PolynomialLCM, 10 bivariate arguments", () -> lcm(ps10));
check("PolynomialLCM, 10 bivariate arguments", polyChk(lcm(ps10), v3));

benchLoop("PolynomialLCM, Modulus 32003, trivariate", () -> modLcm(tri6a, tri6b));
check("PolynomialLCM, Modulus 32003, trivariate", polyChkMod(modLcm(tri6a, tri6b), v3));

\\ ---- no Groebner bases, division by a list, elimination or CAD in PARI -------

skip("PolynomialReduce, by the cyclic-5 DRL basis", "multivariate division");
skip("PolynomialReduce, by the katsura-4 DRL basis", "multivariate division");
skip("PolynomialReduce, by a 3-variable Lex basis", "multivariate division");
skip("PolynomialReduce, by a 3-variable DegLex basis", "multivariate division");
skip("Eliminate, 2 variables from 3 quadratics", "elimination");
skip("Eliminate, plane curve of degree 5", "elimination");
skip("Eliminate, parametric surface", "elimination");
skip("Eliminate, 3 variables from a 4-variable chain", "elimination");
skip("CylindricalDecomposition, unit disc", "cylindrical decomposition");
skip("CylindricalDecomposition, two ellipses", "cylindrical decomposition");
skip("CylindricalDecomposition, quartic curve", "cylindrical decomposition");
skip("CylindricalDecomposition, ball and cylinder", "cylindrical decomposition");
skip("CylindricalDecomposition, Whitney umbrella in a ball", "cylindrical decomposition");
skip("GroebnerBasis, 3-variable quadratics, Lex", "groebner basis");
skip("GroebnerBasis, cyclic-5, DRL", "groebner basis");
skip("GroebnerBasis, katsura-4, DRL", "groebner basis");
skip("GroebnerBasis, elimination form, plane curve", "groebner basis");
skip("GroebnerBasis, cyclic-5, Modulus 32003", "groebner basis");
skip("GroebnerBasis, katsura-5, DRL", "groebner basis");
skip("GroebnerBasis, cyclic-5, Lex", "groebner basis");
skip("GroebnerBasis, cyclic-6, DRL, Modulus 32003", "groebner basis");
skip("GroebnerBasis, cyclic-6, DRL", "groebner basis");
skip("GroebnerBasis, katsura-4, Lex", "groebner basis");
skip("GroebnerBasis, katsura-6, DRL", "groebner basis");

quit();
