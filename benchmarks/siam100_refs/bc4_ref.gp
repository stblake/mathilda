\\ Holdout B, family C4 (experiment 229) -- reference values, route B (PARI/gp).
\\ Objectives written out afresh in GP.  Starting points: the best three of the
\\ L-BFGS-B grid multistart (bc4_ref.py B -> bc4_starts.gp).  Each start is
\\ polished by Newton's method on the gradient, gradient and Hessian by
\\ derivnum, Hessian by central differences of it, at realprecision 110; bc4.6 keeps x = 512 fixed (edge
\\ minimum) and polishes y only.  Prints the lowest value of the three.
\p 110
g5(t) = sum(i = 1, 5, i * cos((i+1)*t + i));
F = vector(10);
F[1] = v -> sin(10*v[1]) + cos(11*v[2]) + (v[1]^2 + v[2]^2)/4;
F[2] = v -> (v[1]^2 + v[2] - 11)^2 + (v[1] + v[2]^2 - 7)^2 + 3*v[1];
F[3] = v -> exp(sin(30*v[1])) + sin(35*exp(v[2])) + sin(sin(40*v[3])) - sin(10*(v[1]+v[2]+v[3])) + (v[1]^2+v[2]^2+v[3]^2)/4;
F[4] = v -> -20*exp(-sqrt(((v[1] - 3/10)^2 + (v[2] + 2/5)^2)/2)/5) - exp((cos(2*Pi*v[1]) + cos(2*Pi*v[2]))/2) + exp(1) + 20;
F[5] = v -> g5(v[1]) * g5(v[2]);
F[6] = v -> -(v[2] + 47)*sin(sqrt(abs(v[1]/2 + v[2] + 47))) - v[1]*sin(sqrt(abs(v[1] - (v[2] + 47))));
F[7] = v -> cos(v[1]) + cos(sqrt(2)*v[1]) + cos(sqrt(3)*v[1]);
F[8] = v -> (sin(25*v[1]) + sin(25*v[2]))/2 + (v[1]^2 + v[2]^2)/10 - 2*exp(-300*((v[1] - 3/5)^2 + (v[2] + 7/20)^2));
F[9] = v -> 100*(v[2] - v[1]^2)^2 + (1 - v[1])^2 + 3/10*sin(15*v[1])*sin(15*v[2]);
F[10] = v -> sum(i = 1, 4, v[i]^2 - 10*cos(2*Pi*v[i])) + (v[1]*v[2] + v[2]*v[3] + v[3]*v[4])/2 + 3/10*v[1];
setc(v, i, t) = my(w = v); w[i] = t; w;
grad(f, v, act) = vector(#act, a, derivnum(t = v[act[a]], f(setc(v, act[a], t))));
\\ Hessian: central differences (h = 1e-25) of the full-precision derivnum gradient
hess(f, v, act) = my(h = 1e-25, H = matrix(#act, #act)); for (b = 1, #act, my(gp = grad(f, setc(v, act[b], v[act[b]] + h), act), gm = grad(f, setc(v, act[b], v[act[b]] - h), act)); for (a = 1, #act, H[a, b] = (gp[a] - gm[a]) / (2*h))); H;
newton(f, v0, act) = {
  my(v = v0, d);
  for (it = 1, 60,
    d = matsolve(hess(f, v, act), grad(f, v, act)~);
    for (a = 1, #act, v[act[a]] = precision(v[act[a]] - d[a], 115));   \\ keep full precision
    if (normlp(d, oo) < 1e-100, break));
  [f(v), v];
}
read("bc4_starts.gp");
{
for (k = 1, 10,
  my(best = [10^9], act = if (k == 6, [2], vector(#starts[k][1], i, i)));
  for (s = 1, #starts[k],
    my(v = starts[k][s] * 1., r);
    if (k == 6, v[1] = 512);
    r = newton(F[k], v, act);
    if (r[1] < best[1], best = r));
  printf("bc4.%d B %.80g\n   at %s\n", k, best[1], apply(t -> precision(t, 20), best[2])));
}
\\ bc4.6 KKT: df/dx at the edge must be negative (f still decreasing outward)
print("bc4.6 df/dx at edge: ", derivnum(t = 512, F[6]([t, newton(F[6], [512., 404.2318], [2])[2][2]])) * 1.);
quit
