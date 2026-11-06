(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 136 -- ISSAC '97 holdout, family 1: ill-conditioned matrices
   ==========================================================================

       cd benchmarks/136-issac97h-illcond-matrices
       ../../Mathilda -file illcond_matrices.m

   HOLDOUT.  Ten problems probing what experiment 126 (the Hilbert condition
   number) probes -- extreme eigenvalues and singular values of exactly
   specified, severely ill-conditioned matrices at a precision far beyond a
   double -- on different matrices, so a fix tuned to the Hilbert matrix
   shows up here as a failure.  See benchmarks/ISSAC97_HOLDOUT.md.

   Every arbitrary-precision case carries the exact matrix to N[m, d] with d
   at least the decimal exponent range of the answer plus 20 digits.  Every
   machine case is a REDUCED instance (smaller n, or milder parameter) whose
   condition number a double still resolves to the 4 digits checked.

   PROBLEMS AND TRUE VALUES.  Each value was computed twice, by two different
   algorithms at two different precisions, and the two agree to every digit
   shown (and, for 1.4, 1.6 and 1.10, to the last digit carried):
     (A) eigenvalues (or singular values) at 400 digits;
     (B) a different route at 600 digits: singular values for eigenvalue
         problems, eigenvalues of the inverse / of A^T A, or the roots of the
         exact characteristic polynomial (1.5, 1.6).
   h1.1  kappa_2 of the 100 x 100 symmetric Pascal matrix Binomial[i+j-2, j-1]
         = 9.191581826323644019411559631e116        machine n = 12: 8.763949109e11
   h1.2  smallest singular value of the 50 x 50 Vandermonde matrix on the
         equispaced nodes (i-1)/49 = 1.214052974061843949823397917e-42
                                                    machine n = 12: 5.670141114e-9
   h1.3  smallest eigenvalue of the 60 x 60 Cauchy matrix 1/(i^2 + j^2)
         = 1.580934514544763948518980013e-63        machine n = 8:  1.147829916e-8
   h1.4  smallest singular value of the 120 x 120 Kahan matrix, theta = 6/5
         (diag Sin[t]^(i-1), -Cos[t] Sin[t]^(i-1) above the diagonal)
         = 4.488013531293049259832894189e-20        machine n = 40: 6.928154746e-7
   h1.5  smallest eigenvalue of the 40 x 40 Frank matrix (n + 1 - Max[i, j]
         for j >= i - 1) = 0.007509065981719966651839414854
                                                    machine n = 12: 0.03102806064
   h1.6  gap between the two largest eigenvalues of the Wilkinson matrix
         W_41^+ (diag |i - 21|, unit off-diagonals)
         = 1.272072120490893356864979269e-37        machine n = 13: 2.021759535e-6
   h1.7  smallest eigenvalue of the 64 x 64 prolate matrix, w = 1/4
         (1/2 on the diagonal, Sin[Pi k/2]/(Pi k) at offset k)
         = 6.738206080100924396271878802e-48        machine n = 12: 1.778272860e-8
   h1.8  smallest eigenvalue of the 80 x 80 Moler matrix (i on the diagonal,
         Min[i, j] - 2 off it) = 6.158049892052418768707796020e-48
                                                    machine n = 12: 5.364405880e-7
   h1.9  kappa_2 of the 60 x 60 Kac-Murdock-Szego matrix rho^|i-j|,
         rho = 1 - 10^-30 = 1.199177720852744324270695235e32
                                    machine rho = 1 - 10^-4: 1.196723905e6
   h1.10 smallest eigenvalue of the 40 x 40 Gaussian kernel matrix
         Exp[-((i - j)/5)^2] = 5.423169387591532736872088203e-21
                                                    machine n = 10: 8.651257513e-10
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Eigenvalues", "SingularValueList"}];

pascal[n_] := Table[Binomial[i + j - 2, j - 1], {i, n}, {j, n}];
vander[n_] := Table[If[j == 1, 1, ((i - 1)/(n - 1))^(j - 1)], {i, n}, {j, n}];
cauchy[n_] := Table[1/(i^2 + j^2), {i, n}, {j, n}];
kahan[n_, t_] := Table[Which[i == j, Sin[t]^(i - 1), j > i, -Cos[t] Sin[t]^(i - 1), True, 0], {i, n}, {j, n}];
frank[n_] := Table[If[j >= i - 1, n + 1 - Max[i, j], 0], {i, n}, {j, n}];
wilk[n_] := Table[Which[i == j, Abs[i - (n + 1)/2], Abs[i - j] == 1, 1, True, 0], {i, n}, {j, n}];
prolate[n_] := Table[If[i == j, 1/2, Sin[Pi (i - j)/2]/(Pi (i - j))], {i, n}, {j, n}];
moler[n_] := Table[If[i == j, i, Min[i, j] - 2], {i, n}, {j, n}];
kms[n_, r_] := Table[r^Abs[i - j], {i, n}, {j, n}];
gk[n_, h_] := Table[Exp[-((i - j) h)^2], {i, n}, {j, n}];

(* Sorted ascending; Re drops a 0``k I a nonsymmetric solver may attach. *)
ev[m_, p_] := Sort[Re[Eigenvalues[N[m, p]]]];
sv[m_, p_] := Sort[SingularValueList[N[m, p]]];
cond[m_, p_] := Module[{e = ev[m, p]}, Last[e]/First[e]];
mp = MachinePrecision;

run["h1.1 eigenvalues machine n=12",  cond[pascal[12], mp], sig[#, 4] &, 3];
run["h1.1 eigenvalues wp250 n=100",   cond[pascal[100], 250], sig[#, 12] &];
run["h1.2 singular values machine n=12", First[sv[vander[12], mp]], sig[#, 4] &, 3];
run["h1.2 singular values wp100 n=50",   First[sv[vander[50], 100]], sig[#, 12] &];
run["h1.3 eigenvalues machine n=8",   First[ev[cauchy[8], mp]], sig[#, 4] &, 3];
run["h1.3 eigenvalues wp120 n=60",    First[ev[cauchy[60], 120]], sig[#, 12] &];
run["h1.4 singular values machine n=40",  First[sv[kahan[40, 6/5], mp]], sig[#, 4] &, 3];
run["h1.4 singular values wp60 n=120",    First[sv[kahan[120, 6/5], 60]], sig[#, 12] &];
run["h1.5 eigenvalues machine n=12",  First[ev[frank[12], mp]], sig[#, 4] &, 3];
run["h1.5 eigenvalues wp100 n=40",    First[ev[frank[40], 100]], sig[#, 12] &];
run["h1.6 eigenvalue gap machine n=13", (#[[-1]] - #[[-2]] &)[ev[wilk[13], mp]], sig[#, 4] &, 3];
run["h1.6 eigenvalue gap wp60 n=41",    (#[[-1]] - #[[-2]] &)[ev[wilk[41], 60]], sig[#, 12] &];
run["h1.7 eigenvalues machine n=12",  First[ev[prolate[12], mp]], sig[#, 4] &, 3];
run["h1.7 eigenvalues wp80 n=64",     First[ev[prolate[64], 80]], sig[#, 12] &];
run["h1.8 eigenvalues machine n=12",  First[ev[moler[12], mp]], sig[#, 4] &, 3];
run["h1.8 eigenvalues wp80 n=80",     First[ev[moler[80], 80]], sig[#, 12] &];
run["h1.9 eigenvalues machine rho=1-1e-4", cond[kms[60, 1 - 10^-4], mp], sig[#, 4] &, 3];
run["h1.9 eigenvalues wp60 rho=1-1e-30",   cond[kms[60, 1 - 10^-30], 60], sig[#, 12] &];
run["h1.10 eigenvalues machine n=10", First[ev[gk[10, 1/5], mp]], sig[#, 4] &, 3];
run["h1.10 eigenvalues wp50 n=40",    First[ev[gk[40, 1/5], 50]], sig[#, 12] &];
