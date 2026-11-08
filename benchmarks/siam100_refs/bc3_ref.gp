\\ Holdout B, family C3 (experiment 228) -- reference values, method A (PARI/gp).
\\ sigma_max of each section = sqrt(lambda_max(A^H A)):
\\   * well-separated spectra (sigma2/sigma1 < 0.5): power iteration x <- A^H (A x)
\\     (no Gram matrix formed), until the Rayleigh quotient is stable to 1e-100;
\\   * clustered spectra (bc3.3, 3.4, 3.5, 3.7: sigma2/sigma1 = 0.64 .. 0.9998):
\\     shifted inverse iteration on M = A^H A with the shift mu = (double
\\     estimate)^2 (numpy/LAPACK, ~15 digits), invert once, iterate.
\\ realprecision 120.  Entries built exactly (rationals / exact phases) then
\\ converted.  Method B: siam100_refs/bc3_ref.py (gmpy2, independent code).
default(parisizemax, 4000000000);
\p 120
ent(k, i, j, n) = {
  if (k == 1, 1/(i+j-1),
  k == 2, 1/(i+j-1)!,
  k == 3, 1/(1+(i-j)^2),
  k == 4, if (j <= i, 1/i, 0),
  k == 5, if (i >= j, 2^-(i-j), (-1/3)^(j-i)),
  k == 6, 1/((i+j-1)*(i+j-2)/2 + j),
  k == 7, exp(2*Pi*I*((i*j)%7)/7)/(i+j),
  k == 8, if (j >= i-1, n+1-max(i,j), 0),
  k == 9, 1/(2^i + 3^j),
  k == 10, 1/sqrt(i+j-1));
}
sec(k, m, n) = matrix(m, n, i, j, ent(k, i, j, n) * 1.);
rq(M, x) = real(conj(x~) * M * x) / norml2(x);
powsig(A) = {
  my(x = vector(matsize(A)[2], i, 1.)~, B = conj(A~), l0 = 0., l1, it = 0);
  while (1, my(y = A * x, z = B * y);
    l1 = norml2(y) / norml2(x);
    x = z / sqrt(norml2(z));
    it++; if (abs(l1 - l0) < 1e-105 * l1 && it > 5, break); l0 = l1);
  [sqrt(l1), it];
}
invsig(A, s0) = {
  my(M = conj(A~) * A, n = #M, Fi, x, l0 = 0., l1, it = 0);
  Fi = (M - s0^2 * matid(n))^-1;
  x = vector(n, i, 1.)~;
  while (1, x = Fi * x; x = x / sqrt(norml2(x));
    l1 = rq(M, x); it++;
    if (abs(l1 - l0) < 1e-105 * l1 && it > 2, break); l0 = l1);
  [sqrt(l1), it];
}
\\ [k, tag, m, n, method, double estimate]
jobs = [  [1,"m",500,500,0,0], [1,"w",60,60,0,0],  [2,"m",40,40,0,0],  [3,"m",400,400,1,3.135326220838264], [3,"w",60,60,1,3.038012952371806],  [4,"m",500,500,1,1.7192159353130747], [4,"w",80,80,1,1.6144112097670473],  [5,"m",400,400,1,1.7499030315692388], [5,"w",60,60,1,1.7459121257721564],  [6,"m",500,150,0,0], [6,"w",100,30,0,0],  [7,"m",300,300,1,1.07493889956485], [7,"w",50,50,1,1.007912280679732],  [8,"m",300,300,0,0], [8,"w",60,60,0,0],  [9,"m",150,150,0,0],  [10,"m",1000,1000,0,0], [10,"w",120,120,0,0]];
{
for (q = 1, #jobs, my(J = jobs[q], A = sec(J[1], J[3], J[4]), r);
  r = if (J[5], invsig(A, J[6]), powsig(A));
  printf("bc3.%d %s %d %d it=%d %.100g\n", J[1], J[2], J[3], J[4], r[2], r[1]));
}
\\ convergence of the two infinite-operator rows: the section norm at a larger n
printf("bc3.2 n=60 %.100g\n", powsig(sec(2, 60, 60))[1]);
printf("bc3.9 n=220 %.100g\n", powsig(sec(9, 220, 220))[1]);
quit
