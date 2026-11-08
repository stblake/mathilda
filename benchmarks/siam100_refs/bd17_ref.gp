\\ Reference B for holdout-B family D17 (exp 252), PARI/gp at 100 digits.
\\ Bessel zeros by solve() on brackets; bd17.6 as twice the half-curve
\\ x in [2/3, 1] of cos(pi x) cos(pi y) = -1/4 (symmetric under (x,y) -> (1-y,1-x)).
default(realprecision, 100);
pr(lab, v) = print(lab, "  ", precision(v, 70));
J(m, a, b) = solve(x = a, b, besselj(m, x));
j01 = J(0, 2, 3); j02 = J(0, 5, 6); j03 = J(0, 8.5, 9);
j11 = J(1, 3.5, 4); j12 = J(1, 6.8, 7.2); j21 = J(2, 5, 5.3); j22 = J(2, 8.2, 8.6);
yp(x) = my(u = -1/(4*cos(Pi*x))); sin(Pi*x)/(4*cos(Pi*x)^2*sqrt(1 - u^2));
pr("bd17.1", 1);
pr("bd17.2", 2);
pr("bd17.3", 2*Pi*j01/j02);
pr("bd17.4", Pi*j11/j12);
pr("bd17.5", sqrt(2)/2);
pr("bd17.6", 2*intnum(x = 2/3, 1, sqrt(1 + yp(x)^2)));
pr("bd17.7", Pi/2*j21/j22);
pr("bd17.8", 3/2);
pr("bd17.9", 3/2);
pr("bd17.10", 2*Pi*(j01 + j02)/j03);
quit;
