\\ Experiment 121 -- Together, Cancel and Apart at stress sizes (PARI/GP column).
\\ See rational_functions.m.
\\
\\ PARI's rational functions are reduced on construction, so Together is the
\\ sum itself and Cancel is the quotient: building the value IS the operation.
\\ PARI has no general partial-fraction decomposition, so the Apart cases are
\\ SKIP.  The checks are the exact integers of the .m half.
\\ One statement per top-level line: gp's file reader runs a line at a time.

read("../harness.gp");

require(["sum", "denominator", "poldegree", "substvec"]);

\\ Variables in this order so x has the highest priority.
vv = [x, y, z, u, v, w];
pv = [11/10, 1/3, 2/7, 3/11, 5/13, 7/17];
chk(t) = 10^40 * poldegree(denominator(t), x) + round(10^6 * substvec(t, vv, pv));
denseUPoly(n) = sum(k = 0, n, ((k^2 + 3*k + 1) % 17) * x^k) + 1;
factorableUPoly(k) = prod(i = 1, k, x^2 + i*x + (i + 1));

\\ ---- Together --------------------------------------------------------------

benchLoop("Together, 50 linear summands", () -> sum(i = 1, 50, (i^2 + 1)/(x - i)));
check("Together, 50 linear summands", chk(sum(i = 1, 50, (i^2 + 1)/(x - i))));

benchLoop("Together, 200 linear summands", () -> sum(i = 1, 200, (i + 1)/(x - i)));
check("Together, 200 linear summands", chk(sum(i = 1, 200, (i + 1)/(x - i))));

togMV() = sum(i = 1, 12, (x + i*y)/((x + i*y + z)*(y - i*z + 1)));
benchLoop("Together, 12 trivariate summands", togMV);
check("Together, 12 trivariate summands", chk(togMV()));

nest32() = my(t = x); for (k = 1, 32, t = x + 1/(1 + t)); t;
benchLoop("Together, nested depth 32", nest32);
check("Together, nested depth 32", chk(nest32()));

\\ ---- Cancel ----------------------------------------------------------------

sh200 = factorableUPoly(100);
c400n = sh200 * denseUPoly(200);
c400d = sh200 * (denseUPoly(199) + 1);
benchLoop("Cancel deg-400, deg-200 gcd", () -> c400n / c400d);
check("Cancel deg-400, deg-200 gcd", chk(c400n / c400d));

sh500 = factorableUPoly(250);
c1000n = sh500 * denseUPoly(500);
c1000d = sh500 * (denseUPoly(499) + 1);
benchLoop("Cancel deg-1000, deg-500 gcd", () -> c1000n / c1000d);
check("Cancel deg-1000, deg-500 gcd", chk(c1000n / c1000d));

g3 = (1 + x + y + z)^6;
c3n = g3 * (x + 2*y - z + 3)^4;
c3d = g3 * (x - y + 2*z + 1)^4;
benchLoop("Cancel dense trivariate, deg-6 gcd", () -> c3n / c3d);
check("Cancel dense trivariate, deg-6 gcd", chk(c3n / c3d));

g6 = 1 + x*y^3 + z^2*u*v^4 + w^5 + x*w;
c6n = g6^2 * (x^2 + z*v + 2) * (y*u + w^2 + 1);
c6d = g6^2 * (x*z + u^3 + 3) * (v*w + y^2 + 1);
benchLoop("Cancel sparse 6-variable", () -> c6n / c6d);
check("Cancel sparse 6-variable", chk(c6n / c6d));

\\ ---- Apart -----------------------------------------------------------------

skip("Apart, 40 linear factors", "partial fractions");
skip("Apart, multiplicity 10 and 10", "partial fractions");
skip("Apart, 6 irreducible quadratics", "partial fractions");
skip("Apart, 5 parameters", "partial fractions");

quit();
