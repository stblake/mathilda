\\ Holdout B, family C1 (experiment 226) -- reference values, method B (PARI/gp).
\\ bc1.3 and bc1.10 have no closed form: here each is a sum of Gauss-Legendre
\\ (intnum) pieces between consecutive zeros, accelerated by sumalt
\\ (Cohen-Villegas-Zagier) -- independent of mpmath's quadosc (nsum/Levin).
\\ Also re-derives the closed forms of the other eight with PARI's own special functions.
\p 90
\\ 3. Int_0^1 cos(log(x)/x) dx = Int_0^oo cos(t e^t) e^-t dt; zeros t e^t = (n-1/2) Pi
z3(n) = lambertw((n - 1/2) * Pi);
f3(t) = cos(t * exp(t)) * exp(-t);
v3 = intnum(t = 0, z3(1), f3(t)) + sumalt(n = 1, intnum(t = z3(n), z3(n+1), f3(t)));
\\ 10. Int_2^oo sin(x)/log(x) dx; zeros k Pi
f10(x) = sin(x) / log(x);
v10 = intnum(x = 2, Pi, f10(x)) + sumalt(k = 1, intnum(x = k*Pi, (k+1)*Pi, f10(x)));
ci1 = -eint1(-I) ; \\ E1(-i) = -Ci(1) + i(Si(1) - pi/2)... take parts below
Ci1 = -real(eint1(I)); Si1 = Pi/2 + imag(eint1(I));
print("bc1.1  ", sin(1) - Ci1);
print("bc1.2  ", Pi/2 - Si1);
print("bc1.3  ", v3);
print("bc1.4  ", sqrt(Pi/8));
print("bc1.5  ", sqrt(Pi/2));
\\ 6: Int_1^oo sin t t^(-3/2) dt = Im( Gamma(-1/2, -i) * (-i)^(1/2) )  (incomplete gamma)
print("bc1.6  ", imag(incgam(-1/2, -I) * (-I)^(1/2)));
print("bc1.7  ", Ci1 * sin(1) + (Pi/2 - Si1) * cos(1));
print("bc1.8  ", gamma(3/4) / gamma(5/4) / sqrt(2));
print("bc1.9  ", gamma(4/3) * cos(Pi/6));
print("bc1.10 ", v10);
quit
