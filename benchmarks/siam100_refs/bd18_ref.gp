\\ Reference B for holdout-B family D18 (exp 253), PARI/gp at 100 digits:
\\ Bessel zeros by solve() on a sign-change bracket of besselj / besseln.
default(realprecision, 100);
pr(lab, v) = print(lab, "  ", precision(v, 70));
j01 = solve(x = 2, 3, besselj(0, x));
pr("bd18.1", 2*Pi^2);
pr("bd18.2", 5*Pi^2);
pr("bd18.3", 112*Pi^2/9);
pr("bd18.4", j01^2);
pr("bd18.5", solve(x = 2.5, 3.5, besselj(2/3, x))^2);
pr("bd18.6", solve(x = 5.5, 7, besselj(0, x/2)*besseln(0, x) - besselj(0, x)*besseln(0, x/2))^2);
pr("bd18.7", 3*Pi^2);
pr("bd18.8", solve(x = 3, 3.3, besselj(1/2, x))^2);
pr("bd18.9", j01^2 + Pi^2);
pr("bd18.10", solve(x = 2.5, 3.5, besselj(4/7, x))^2);
quit;
