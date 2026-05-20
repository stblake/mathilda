#include "info.h"
#include "symtab.h"

void info_init(void) {
    // Arithmetic
    symtab_set_docstring("Plus", "x + y + ... or Plus[x, y, ...] represents a sum of terms.");
    symtab_set_docstring("Times", "x * y * ... or Times[x, y, ...] represents a product of terms.");
    symtab_set_docstring("Power", "x ^ y or Power[x, y] represents x to the power y.");
    symtab_set_docstring("Subtract", "x - y or Subtract[x, y] represents x - y.");
    symtab_set_docstring("Divide", "x / y or Divide[x, y] represents x / y.");
    symtab_set_docstring("Sqrt", "Sqrt[z] represents the square root of z.");
    symtab_set_docstring("Complex", "Complex[re, im] represents a complex number with real part re and imaginary part im.");
    symtab_set_docstring("Rational", "Rational[n, d] represents a rational number with numerator n and denominator d. Automatically simplifies if arguments are integers.");
    symtab_set_docstring("GCD", "GCD[n1, n2, ...] gives the greatest common divisor of the integers ni.");
    symtab_set_docstring("LCM", "LCM[n1, n2, ...] gives the least common multiple of the integers ni.");
    symtab_set_docstring("PowerMod", "PowerMod[a, b, m] gives a^b mod m.\nPowerMod[a, -1, m] finds the modular inverse of a modulo m.\nPowerMod[a, 1/r, m] finds a modular r-th root of a.");
    symtab_set_docstring("Factorial", "n! gives the factorial of n.\nFor integers and half integers, Factorial automatically evaluates to exact values.");
    symtab_set_docstring("Factorial2", "Factorial2[n] (also typeset n!!) gives the double factorial of n.\nFor non-negative integer n: n!! = n * (n-2) * (n-4) * ... down to 2 (n even) or 1 (n odd).\nSpecial values: 0!! = 1, (-1)!! = 1.\nNegative even integers and negative odd integers below -1 give ComplexInfinity.\nFactorial2 stays unevaluated on symbolic arguments.");
    symtab_set_docstring("Binomial", "Binomial[n, m] gives the binomial coefficient.");

    // Structural Manipulation
    symtab_set_docstring("Part",
        "expr[[i]] or Part[expr, i]\n"
        "\tgives the i-th part of expr.\n"
        "expr[[-i]]\n"
        "\tcounts from the end.\n"
        "expr[[0]]\n"
        "\tgives the head of expr.\n"
        "expr[[i, j, ...]] or Part[expr, i, j, ...]\n"
        "\tis equivalent to expr[[i]][[j]]..., descending into nested parts.\n"
        "expr[[{i1, i2, ...}]]\n"
        "\tgives a list of the parts i1, i2, ... of expr (wrapped in the head of expr).\n"
        "expr[[m;;n]] / expr[[m;;n;;s]]\n"
        "\tgives the span of parts m through n (with optional step s); ;; alone or All means all parts.\n\n"
        "Part is treated as atomic on Integer, Real, String, Symbol, Rational[n, d], and Complex[re, im]; "
        "Part[atom, i] for i != 0 stays unevaluated.\n"
        "Indices are 1-based and may be negative; out-of-range indices leave the expression unevaluated.");
    symtab_set_docstring("Extract",
        "Extract[expr, pos]\n"
        "\textracts the part of expr at the position specified by pos.\n"
        "Extract[expr, {pos1, pos2, ...}]\n"
        "\textracts a list of parts of expr.\n"
        "Extract[expr, pos, h]\n"
        "\textracts parts of expr, wrapping each of them with head h before evaluation.\n"
        "Extract[pos]\n"
        "\trepresents an operator form of Extract that can be applied to an expression.\n\n"
        "The position pos has the same form as a Position result: a list of indices "
        "{i1, i2, ...} that descends i1 into expr, then i2, etc. A scalar index n is "
        "treated as the path {n}, so Extract[expr, n] is equivalent to Extract[expr, {n}]; "
        "in particular Extract[expr, 0] gives Head[expr].\n"
        "Indices are 1-based and may be negative; index 0 selects the head. "
        "Extract is treated as atomic on Integer, Real, String, Symbol, "
        "Rational[n, d], and Complex[re, im].");
    symtab_set_docstring("Span", "i;;j represents a span of elements i through j. i;;j;;k represents a span in steps of k.");
    symtab_set_docstring("UpTo",
        "UpTo[n]\n"
        "\tis a symbolic specification that represents up to n objects or positions. If n objects or positions are available, all are used. If fewer are available, only those available are used.");


    // Linear Algebra
    symtab_set_docstring("Dot", "a.b.c or Dot[a, b, c] gives products of vectors, matrices, and tensors.");
    symtab_set_docstring("Det", "Det[m] gives the determinant of the square matrix m.");
    symtab_set_docstring("Cross", "Cross[a, b, ...] gives the vector cross product of the arguments.");
    symtab_set_docstring("Norm", "Norm[expr] gives the norm of a number, vector, or matrix.\nNorm[expr, p] gives the p-norm.");
    symtab_set_docstring("Tr", "Tr[list] finds the trace of the matrix or tensor list.\nTr[list, f] finds a generalized trace, combining terms with f instead of Plus.\nTr[list, f, n] goes down to level n in list.");
    symtab_set_docstring("RowReduce",
        "RowReduce[m]\n"
        "\tgives the row-reduced form of the matrix m.\n"
        "RowReduce[m, Method -> \"<name>\"]\n"
        "\truns a specific elimination algorithm.  Accepted method names:\n"
        "\t  \"Automatic\"                 — alias for \"DivisionFreeRowReduction\" (default)\n"
        "\t  \"DivisionFreeRowReduction\"  — Bareiss-like fraction-free Gauss-Jordan\n"
        "\t  \"OneStepRowReduction\"       — classical Gauss-Jordan with division per pivot\n"
        "\t  \"CofactorExpansion\"         — identity-if-invertible via Laplace cofactor\n"
        "\t                                 Det[m] (singular / rectangular m falls back\n"
        "\t                                 to \"DivisionFreeRowReduction\")");
    symtab_set_docstring("IdentityMatrix", "IdentityMatrix[n] gives the n x n identity matrix.\nIdentityMatrix[{m, n}] gives the m x n identity matrix.");
    symtab_set_docstring("DiagonalMatrix", "DiagonalMatrix[list] gives a matrix with the elements of list on the leading diagonal, and zero elsewhere.\nDiagonalMatrix[list, k] gives a matrix with the elements of list on the k-th diagonal.\nDiagonalMatrix[list, k, n] pads with zeros to create an n x n matrix.");
    symtab_set_docstring("Inverse",
        "Inverse[m]\n"
        "\tgives the inverse of a square matrix m.\n"
        "Inverse[m, Method -> \"<name>\"]\n"
        "\truns a specific inversion algorithm.\n"
        "\n"
        "Inverse works on both symbolic and numerical matrices.\n"
        "For matrices with approximate real or complex numbers, the\n"
        "inverse is generated to the maximum possible precision given the\n"
        "input.  Inverse::sing is issued for singular matrices and\n"
        "Inverse::matsq for non-square / empty input; in either case the\n"
        "call is returned unevaluated.\n"
        "\n"
        "Accepted method names:\n"
        "  \"Automatic\"                 — alias for \"DivisionFreeRowReduction\" (default)\n"
        "  \"DivisionFreeRowReduction\"  — Bareiss-like fraction-free Gauss-Jordan on [m | I]\n"
        "  \"OneStepRowReduction\"       — classical Gauss-Jordan with division per pivot\n"
        "  \"CofactorExpansion\"         — adjugate / determinant formula via Laplace expansion\n"
        "\n"
        "An unknown method name emits Inverse::method and leaves the call\n"
        "unevaluated.  Method -> Automatic (the symbol) is also accepted.");
    symtab_set_docstring("PseudoInverse",
        "PseudoInverse[m]\n"
        "\tfinds the Moore-Penrose pseudoinverse of a rectangular matrix m.\n"
        "PseudoInverse[m, Tolerance -> t]\n"
        "\tspecifies that singular values smaller than t times the maximum\n"
        "\tsingular value should be dropped.  With the default setting\n"
        "\tTolerance -> Automatic, the rationalisation precision of the\n"
        "\tinput is used (Real -> 53 bits, MPFR -> input precision).\n"
        "\n"
        "For non-singular square matrices m, the pseudoinverse coincides\n"
        "with the standard inverse: PseudoInverse[m] == Inverse[m].\n"
        "\n"
        "PseudoInverse works on exact (Integer / Rational / Complex)\n"
        "matrices and on approximate (Real / MPFR) matrices.  For exact\n"
        "input the result is exact; for inexact input the input is\n"
        "rationalised, the pseudoinverse is computed in exact arithmetic\n"
        "via a full-rank decomposition, and the result is numericalised\n"
        "back to the input precision.\n"
        "\n"
        "Algorithm: row-reduce m to identify rank r and a full-rank\n"
        "decomposition m = B . C with B m x r and C r x n.  Then\n"
        "    PseudoInverse[m] = ConjugateTranspose[C] . Inverse[C . ConjugateTranspose[C]]\n"
        "                                            . Inverse[ConjugateTranspose[B] . B]\n"
        "                                            . ConjugateTranspose[B].\n"
        "When m is the zero matrix the pseudoinverse is the corresponding\n"
        "zero matrix of transposed shape.");
    symtab_set_docstring("MatrixPower", "MatrixPower[m, n]\n\tgives the n-th matrix power of the square matrix m.\n\tMatrixPower[m, n, v] gives the n-th matrix power of the matrix m applied to the vector v.\n\tWhen n is negative, MatrixPower finds powers of the inverse of the matrix m.\n\tMatrixPower[m, 0] gives IdentityMatrix[Length[m]].\n\tFractional matrix powers are not currently supported.");
    symtab_set_docstring("Eigenvalues",
        "Eigenvalues[m]\n"
        "\tgives a list of the eigenvalues of the square matrix m.\n"
        "Eigenvalues[{m, a}]\n"
        "\tgives the generalized eigenvalues of m with respect to a.\n"
        "Eigenvalues[m, k]\n"
        "\tgives the first k eigenvalues (largest by absolute value).\n"
        "Eigenvalues[m, -k]\n"
        "\tgives the k eigenvalues smallest in absolute value.\n"
        "Eigenvalues[m, UpTo[k]]\n"
        "\tgives k eigenvalues, or as many as are available.\n"
        "\n"
        "Eigenvalues are computed from the roots of the characteristic\n"
        "polynomial Det[m - lambda I] (or Det[m - lambda a] for the\n"
        "generalised case). Approximate (Real / MPFR) matrices flow through\n"
        "the Solve rationalise -> solve -> numericalize pipeline and yield\n"
        "numerical eigenvalues sorted in order of decreasing absolute value.\n"
        "Repeated eigenvalues appear with their algebraic multiplicity.\n"
        "\n"
        "Options:\n"
        "    Cubics    -> True       (use radicals to solve cubics)\n"
        "    Quartics  -> True       (use radicals to solve quartics)\n"
        "    Method    -> Automatic  (numeric-matrix method dispatch)\n"
        "\n"
        "Method values for approximate-numeric matrices:\n"
        "  Automatic    selects Direct unless k is small (-> Arnoldi)\n"
        "               or the matrix is Hermitian-banded (-> Banded).\n"
        "  \"Direct\"     Hessenberg + implicit shifted QR (general); for\n"
        "               Hermitian inputs tridiagonalisation + Wilkinson-\n"
        "               shift symmetric QR.  Returns all eigenvalues.\n"
        "  \"Arnoldi\"    Krylov-subspace iteration for the k extreme\n"
        "               eigenvalues; accepts Method -> {\"Arnoldi\",\n"
        "               MaxIterations -> n, Tolerance -> t, BasisSize -> m}.\n"
        "  \"Banded\"     Hermitian only; auto-detects band structure and\n"
        "               reduces to tridiagonal before symmetric QR.\n"
        "  \"FEAST\"      Hermitian only; eigenvalues in a user-specified\n"
        "               Interval -> {a, b}; accepts Method ->\n"
        "               {\"FEAST\", \"Interval\" -> {a, b},\n"
        "                \"ContourPoints\" -> Ne, \"SubspaceSize\" -> m0,\n"
        "                \"MaxIterations\" -> k, \"Tolerance\" -> t}.\n"
        "Non-numeric matrices ignore Method and use the symbolic\n"
        "characteristic-polynomial pipeline.\n"
        "\n"
        "Implementation status: \"Direct\" runs the hand-rolled Householder\n"
        "tridiagonalisation + Wilkinson-shift symmetric QR kernel at\n"
        "machine precision for real symmetric matrices, the Hessenberg\n"
        "+ implicit double-shift Francis QR kernel for real non-symmetric\n"
        "matrices, a complex Householder tridiagonalisation + diagonal-\n"
        "phase correction + symmetric QR kernel for complex Hermitian\n"
        "matrices (returns real eigenvalues sorted by |lambda|\n"
        "descending), and a real-block-embedding kernel for complex\n"
        "non-Hermitian matrices (M = [[Re A, -Im A], [Im A, Re A]]\n"
        "routed through real Hessenberg + Francis QR with grouped\n"
        "complex Gram-Schmidt disambiguation of M's spec to recover\n"
        "spec(A)).  Automatic routes here too.  Arbitrary-precision\n"
        "(MPFR) inputs go through a parallel \"Direct\" kernel at the\n"
        "input's combined precision: all four shapes -- real symmetric\n"
        "(step 2d-A), real non-symmetric (step 2d-B), complex Hermitian\n"
        "(step 2d-C), and complex non-Hermitian (step 2d-D) -- return\n"
        "eigenvalues / eigenvectors carrying full input precision.\n"
        "\"Arnoldi\" is implemented in Phase 3 at both machine and MPFR\n"
        "precision: m-step classical Gram-Schmidt with one re-orthog-\n"
        "onalisation pass builds the Krylov basis V_m and the m x m\n"
        "upper Hessenberg H_m; H_m is diagonalised by reusing the\n"
        "\"Direct\" Francis QR pipeline (real machine, real MPFR, or via\n"
        "a 2mu x 2mu real-block embedding for complex H_m), and Ritz\n"
        "vectors V_m y_i lift back to A-eigenvectors.  Complex inputs\n"
        "use paired re/im storage for V_m and H_m.  Automatic routes to\n"
        "Arnoldi when n > 32 and k_spec is given with k <= max(20, n/10);\n"
        "small matrices always go through Direct (faster + exact).\n"
        "Default BasisSize is max(2k, 20) capped at n; on lucky breakdown\n"
        "(||w|| below tolerance at some step j) Arnoldi terminates early\n"
        "with j+1 exact eigenpairs.  MPFR Arnoldi carries through the\n"
        "input's combined precision via the same scratch-pool discipline\n"
        "as the Direct MPFR kernels.\n"
        "\"Banded\" (Phase 4, machine + MPFR) handles real symmetric and\n"
        "complex Hermitian matrices.  It auto-detects the half-bandwidth\n"
        "and reduces to symmetric tridiagonal form via Schwarz-style\n"
        "two-sided Givens rotations with bulge chasing (one off-band\n"
        "entry zeroed per Givens; the introduced bulge is chased b\n"
        "columns at a time until it falls past the matrix edge); the\n"
        "resulting tridiagonal eigenproblem reuses the Phase 2 symmetric\n"
        "QR.  Complex Hermitian banded uses paired re/im Givens with a\n"
        "real-c / complex-s parameterisation and the same phase-\n"
        "correction step as the Direct Hermitian kernel.  Banded refuses\n"
        "(returns NULL, falls back to Direct) when the matrix isn't\n"
        "Hermitian or when it's fully dense (b == n - 1).  Automatic\n"
        "routes here when the matrix is Hermitian, n > 8, and the half-\n"
        "bandwidth is at most max(8, n/4); narrower bands save more flops\n"
        "than wider ones.\n"
        "\"FEAST\" (Phase 5, machine + MPFR) handles Hermitian (real\n"
        "symmetric or complex Hermitian) input and returns only the\n"
        "eigenvalues in a user-supplied real Interval -> {a, b} -- a\n"
        "spectral-slice query rather than a full decomposition.  Uses\n"
        "Ne-point Gauss-Legendre quadrature (default Ne = 8; supported:\n"
        "2, 4, 8, 16) on the upper half of the elliptic contour through\n"
        "(a, 0) and (b, 0) to approximate the spectral projector\n"
        "P_[a,b](A); Schwarz symmetry halves the number of complex\n"
        "linear solves.  A Rayleigh-Ritz reduction with Cholesky\n"
        "B_q = L L^* extracts the in-interval eigenpairs by reusing the\n"
        "Direct symmetric / Hermitian kernel on L^-1 A_q L^-*.  Output\n"
        "is sorted by |lambda| descending so k_spec composes naturally\n"
        "with the in-interval filter.  Automatic never routes to FEAST\n"
        "(it requires the user to commit to an interval).  FEAST falls\n"
        "back to Direct (NULL return + one-shot stderr warning tagged\n"
        "with the reason) on: non-Hermitian input, missing Interval,\n"
        "degenerate or invalid {a, b} (interval_high <= interval_low),\n"
        "generalised eigenproblem, Cholesky failure on B_q (subspace\n"
        "too small for the spectral count), LU singular at any quad-\n"
        "rature node, or non-convergence within MaxIterations.");
    symtab_set_docstring("MatrixRank",
        "MatrixRank[m]\n"
        "\tgives the rank of the matrix m -- the number of linearly\n"
        "\tindependent rows (equivalently, of linearly independent\n"
        "\tcolumns).\n"
        "MatrixRank[m, Method -> \"<name>\"]\n"
        "\truns a specific elimination algorithm for the exact path.\n"
        "\tAccepted method names match NullSpace / RowReduce /\n"
        "\tLinearSolve / Inverse:\n"
        "\t  \"Automatic\"                 -- alias for \"DivisionFreeRowReduction\" (default)\n"
        "\t  \"DivisionFreeRowReduction\"  -- Bareiss-like fraction-free Gauss-Jordan\n"
        "\t  \"OneStepRowReduction\"       -- classical Gauss-Jordan with division per pivot\n"
        "\t  \"CofactorExpansion\"         -- identity-if-invertible (falls back to\n"
        "\t                                 DivisionFreeRowReduction on singular /\n"
        "\t                                 rectangular m)\n"
        "MatrixRank[m, Tolerance -> t]\n"
        "\ttreats |entry| <= t as zero during pivot selection.  With\n"
        "\tTolerance -> 0 even arbitrarily small entries count; the\n"
        "\tdefault, Tolerance -> Automatic, applies\n"
        "\tmax(rows, cols) * MachineEpsilon * Max[|entries|] for\n"
        "\tmachine-precision (Real / MPFR) matrices and 0 otherwise.\n"
        "\n"
        "MatrixRank works on both numerical and symbolic matrices and\n"
        "on square or rectangular matrices.  The default exact path\n"
        "routes through RowReduce and counts the non-zero rows; the\n"
        "numerical path (triggered by inexact leaves or an explicit\n"
        "Tolerance) runs partial-pivot Gaussian elimination over\n"
        "double-precision complex.\n"
        "\n"
        "An unknown Method value or Tolerance form emits\n"
        "MatrixRank::opt and leaves the call unevaluated.  A non-rank-2\n"
        "or empty matrix emits MatrixRank::matrix and the call is left\n"
        "unevaluated.");
    symtab_set_docstring("NullSpace",
        "NullSpace[m]\n"
        "\tgives a list of vectors that forms a basis for the null\n"
        "\tspace of the matrix m (i.e. vectors v such that m . v == 0).\n"
        "NullSpace[m, Method -> \"<name>\"]\n"
        "\truns a specific elimination algorithm.  Accepted method\n"
        "\tnames are the same as RowReduce / LinearSolve / Inverse:\n"
        "\t  \"Automatic\"                 — alias for \"DivisionFreeRowReduction\" (default)\n"
        "\t  \"DivisionFreeRowReduction\"  — Bareiss-like fraction-free Gauss-Jordan\n"
        "\t  \"OneStepRowReduction\"       — classical Gauss-Jordan with division per pivot\n"
        "\t  \"CofactorExpansion\"         — identity-if-invertible (falls back to\n"
        "\t                                 DivisionFreeRowReduction on singular /\n"
        "\t                                 rectangular m)\n"
        "\n"
        "NullSpace works on both numerical and symbolic matrices.  The\n"
        "matrix m may be square or rectangular.  When m has full column\n"
        "rank the result is the empty list `{}`.\n"
        "\n"
        "Basis vectors are returned with the rightmost free column\n"
        "first.  For exact integer / rational input each basis vector\n"
        "is scaled to clear integer denominators, so the result is\n"
        "integer-valued whenever the input is integer-valued.  For\n"
        "symbolic input the basis vectors are left in their natural\n"
        "rational form.\n"
        "\n"
        "An unknown method name emits NullSpace::method and leaves the\n"
        "call unevaluated.  A non-rank-2 / empty matrix emits\n"
        "NullSpace::matrix and the call is left unevaluated.");
    symtab_set_docstring("LinearSolve",
        "LinearSolve[m, b]\n"
        "\tfinds an x that solves the matrix equation m . x == b.\n"
        "LinearSolve[m, b, Method -> \"<name>\"]\n"
        "\truns a specific elimination algorithm.\n"
        "\n"
        "LinearSolve works on both numerical and symbolic matrices.\n"
        "The matrix m may be square or rectangular.\n"
        "The argument b may be a vector or a matrix; when b is a matrix\n"
        "(one column per RHS) LinearSolve returns a matrix of solutions.\n"
        "Higher-rank tensor inputs are also supported: when m has\n"
        "dimensions {d1, ..., d(N-1), n}, b may have dimensions\n"
        "{d1, ..., d(N-1), e1, ..., ep} and the result has dimensions\n"
        "{n, e1, ..., ep}.\n"
        "\n"
        "For under-determined systems LinearSolve returns a particular\n"
        "solution in which the free (non-pivot) variables are taken to be\n"
        "0; Solve returns the general solution.  When the equation has no\n"
        "solution LinearSolve emits LinearSolve::nosol and returns\n"
        "unevaluated.\n"
        "\n"
        "Accepted method names:\n"
        "  \"Automatic\"                 — alias for \"DivisionFreeRowReduction\" (default)\n"
        "  \"DivisionFreeRowReduction\"  — Bareiss-like fraction-free Gauss-Jordan on [m | b]\n"
        "  \"OneStepRowReduction\"       — classical Gauss-Jordan with division per pivot\n"
        "  \"CofactorExpansion\"         — Cramer's rule via Laplace cofactor expansion\n"
        "                                 (square non-singular m only; LinearSolve::cofnsq\n"
        "                                 / ::cofsng on shape / singularity errors)\n"
        "\n"
        "Default implementation: fraction-free Gauss-Jordan elimination on\n"
        "the augmented matrix [m | b] (the Bareiss-like algorithm shared\n"
        "with RowReduce and Inverse), so exact integer / rational /\n"
        "symbolic inputs flow through without any spurious denominator\n"
        "blow-up.");
    symtab_set_docstring("LeastSquares",
        "LeastSquares[m, b]\n"
        "\tfinds an x that solves the linear least-squares problem\n"
        "\tfor the matrix equation m . x == b, i.e. an x minimising\n"
        "\tNorm[m . x - b].\n"
        "LeastSquares[m, b, Method -> \"<name>\"]\n"
        "\tselects an explicit solver.\n"
        "LeastSquares[m, b, Tolerance -> t]\n"
        "\tspecifies the singular-value truncation tolerance forwarded\n"
        "\tto the underlying PseudoInverse call (Tolerance -> Automatic\n"
        "\tby default).  Method and Tolerance options may appear in any\n"
        "\torder.\n"
        "\n"
        "LeastSquares works on every input family supported by\n"
        "PseudoInverse: exact (Integer / Rational), symbolic, inexact\n"
        "(Real / MPFR), and complex.  The matrix m may be square or\n"
        "rectangular and of any rank.  When m is rank-deficient the\n"
        "result is the minimum-norm minimiser -- the Moore-Penrose\n"
        "pseudoinverse solution PseudoInverse[m] . b.\n"
        "\n"
        "The right-hand side b may be a vector or a matrix.  When b is\n"
        "a matrix (one column per RHS), LeastSquares returns a matrix\n"
        "of solutions, the j-th column of which is the least-squares\n"
        "solution for the j-th column of b -- minimising\n"
        "Norm[m . x - b, \"Frobenius\"] over the multi-RHS system.\n"
        "\n"
        "Accepted Method names:\n"
        "  \"Automatic\"           — alias for \"Direct\" (default)\n"
        "  \"Direct\"              — PseudoInverse[m] . b; works for all\n"
        "                          input families (dense or sparse,\n"
        "                          exact or numeric, real or complex)\n"
        "  \"IterativeRefinement\" — residual-correction loop on Direct,\n"
        "                          x <- x + PseudoInverse[m] . (b - m.x),\n"
        "                          terminated when ||dx||^2 <= Tolerance^2\n"
        "                          or at a 50-iteration cap.  Exact inputs\n"
        "                          converge in one pass; inexact inputs\n"
        "                          drive round-off down to Tolerance.\n"
        "  \"Krylov\"              — Conjugate-Gradient-on-Least-Squares\n"
        "                          (Hestenes-Stiefel CG on the normal\n"
        "                          equations) with x_0 = 0.  Converges\n"
        "                          to the minimum-norm LS solution for\n"
        "                          rank-deficient m, capped at 2 cols(m)\n"
        "                          + 10 iterations.  Free symbols fall\n"
        "                          back to Direct.\n"
        "  \"LSQR\"                — Paige-Saunders LSQR: Lanczos\n"
        "                          bidiagonalisation with Givens rotations.\n"
        "                          Free symbols fall back to Direct; exact\n"
        "                          rationals and complex inputs fall back\n"
        "                          to Krylov / CGLS (equivalent without\n"
        "                          square-root growth); pure-real inputs\n"
        "                          with at least one Real entry run the\n"
        "                          double-precision algorithm.\n"
        "\n"
        "An unknown Method or Tolerance value leaves the call\n"
        "unevaluated.  When m . x == b has an exact solution,\n"
        "LeastSquares[m, b] coincides with LinearSolve[m, b].");
    symtab_set_docstring("Eigenvectors",
        "Eigenvectors[m]\n"
        "\tgives a list of the eigenvectors of the square matrix m.\n"
        "Eigenvectors[{m, a}]\n"
        "\tgives the generalized eigenvectors of m with respect to a.\n"
        "Eigenvectors[m, k]\n"
        "\tgives the first k eigenvectors.\n"
        "Eigenvectors[m, UpTo[k]]\n"
        "\tgives k eigenvectors, or as many as are available.\n"
        "\n"
        "For an n x n matrix Eigenvectors always returns a list of length n.\n"
        "If a matrix is defective for some eigenvalue, the corresponding\n"
        "shortfall is padded with zero vectors. For approximate numerical\n"
        "matrices the eigenvectors are normalised to unit Norm; for exact\n"
        "or symbolic matrices the eigenvectors are not normalised.\n"
        "\n"
        "Options:\n"
        "    Cubics    -> True       (use radicals to solve cubics)\n"
        "    Quartics  -> True       (use radicals to solve quartics)\n"
        "    Method    -> Automatic  (numeric-matrix method dispatch)\n"
        "\n"
        "Method values for approximate-numeric matrices mirror Eigenvalues:\n"
        "Automatic, \"Direct\", \"Arnoldi\", \"Banded\", and \"FEAST\".  Each\n"
        "method returns the eigenvectors associated with the eigenvalues\n"
        "it would compute.  See ?Eigenvalues for the per-method semantics\n"
        "and sub-option grammar.  Non-numeric matrices ignore Method and\n"
        "use the symbolic null-space pipeline.\n"
        "\n"
        "Implementation status: \"Direct\" yields orthonormal eigenvectors\n"
        "for real symmetric matrices at machine precision (Householder +\n"
        "symmetric QR with accumulated rotations), unit-norm eigenvectors\n"
        "for real non-symmetric matrices via Hessenberg + Francis double-\n"
        "shift QR with accumulated Q followed by Schur-form back-\n"
        "substitution (complex eigenvalues yield complex eigenvectors\n"
        "emitted as Complex[re, im] entries), unitary orthonormal complex\n"
        "eigenvectors for complex Hermitian matrices via complex\n"
        "Householder tridiagonalisation + diagonal-phase correction +\n"
        "symmetric QR with composed complex Q, and unit-norm complex\n"
        "eigenvectors for complex non-Hermitian matrices via real block\n"
        "embedding into a 2n x 2n general matrix followed by grouped\n"
        "complex Gram-Schmidt extraction.  Automatic routes here.\n"
        "Arbitrary-precision (MPFR) inputs run a parallel \"Direct\" kernel\n"
        "at the input's combined precision: real symmetric (step 2d-A),\n"
        "real non-symmetric (step 2d-B), complex Hermitian (step 2d-C),\n"
        "and complex non-Hermitian (step 2d-D) MPFR all yield eigenvectors\n"
        "carrying full input precision -- orthonormal for the Hermitian /\n"
        "symmetric paths, unit 2-norm for the general paths.\n"
        "\"Arnoldi\" (Phase 3, machine + MPFR) returns Ritz vectors\n"
        "V_m y_i where V_m is the orthonormal Arnoldi basis and y_i\n"
        "diagonalises the small m x m Hessenberg H_m.  Ritz vectors are\n"
        "unit 2-norm; for ill-conditioned matrices or m close to the\n"
        "spectral diameter they may need refinement (single inverse\n"
        "iteration is sufficient in practice).  MPFR Arnoldi carries\n"
        "input precision through to all output components.\n"
        "\"Banded\" (Phase 4, machine + MPFR) returns orthonormal real\n"
        "eigenvectors for real symmetric banded inputs and unitary\n"
        "complex eigenvectors for complex Hermitian banded inputs.  The\n"
        "band-Givens reduction accumulates an orthogonal (resp. unitary)\n"
        "Q during the chase; the final Z from the symmetric tridiag QR\n"
        "is composed against Q exactly as in the Direct Hermitian path.\n"
        "Banded refuses (falls back to Direct) on non-Hermitian or fully\n"
        "dense matrices.\n"
        "\"FEAST\" (Phase 5, machine + MPFR) returns the eigenvectors\n"
        "whose eigenvalues lie in the user-supplied real Interval ->\n"
        "{a, b} -- orthonormal for real symmetric input, unitary for\n"
        "complex Hermitian input.  Sub-option grammar mirrors Eigen-\n"
        "values: Method -> {\"FEAST\", \"Interval\" -> {a, b},\n"
        "\"ContourPoints\" -> Ne, \"SubspaceSize\" -> m0,\n"
        "\"MaxIterations\" -> k, \"Tolerance\" -> t}.  Same fail-soft\n"
        "cascade as Eigenvalues -- non-Hermitian, missing / degenerate\n"
        "Interval, generalised problem, Cholesky / LU failure, or\n"
        "non-convergence all fall back to Direct with a one-shot\n"
        "stderr warning.");
    symtab_set_docstring("FullForm", "FullForm[expr] prints as the full internal structure of expr, without any special formatting.");
    symtab_set_docstring("Head",
        "Head[expr]\n"
        "    Gives the head of expr.\n"
        "Head[expr, h]\n"
        "    Wraps the result with h, i.e. returns h[Head[expr]].\n"
        "\n"
        "For atoms, Head returns Integer, Real, Symbol, or String. For a\n"
        "compound expression f[...], Head returns f.\n"
        "\n"
        "Examples:\n"
        "    Head[a + b]      -> Plus\n"
        "    Head[{a, b, c}]  -> List\n"
        "    Head[3.14]       -> Real\n"
        "    Head[x]          -> Symbol\n"
        "    Head[a + b, f]   -> f[Plus]");
    symtab_set_docstring("Length", "Length[expr] gives the number of elements in expr.");
    symtab_set_docstring("Dimensions",
        "Dimensions[expr]\n"
        "    Gives a list of the dimensions of expr.\n"
        "Dimensions[expr, n]\n"
        "    Gives a list of the dimensions of expr down to level n.\n"
        "\n"
        "expr is treated as a full array only at levels where every sub-piece\n"
        "shares the same head and length; ragged levels are not counted.\n"
        "Dimensions always returns a List, including the empty List {} for\n"
        "atomic expressions.\n"
        "\n"
        "Examples:\n"
        "    Dimensions[{{a, b, c}, {d, e, f}}]    -> {2, 3}\n"
        "    Dimensions[{{a, b, c}, {d, e}, {f}}]  -> {3}      (ragged at level 2)\n"
        "    Dimensions[{{{{a, b}}}}]              -> {1, 1, 1, 2}\n"
        "    Dimensions[{{{{a, b}}}}, 2]           -> {1, 1}\n"
        "    Dimensions[1]                         -> {}");
    symtab_set_docstring("First", "First[expr] gives the first element of expr.");
    symtab_set_docstring("Last", "Last[expr] gives the last element of expr.");
    symtab_set_docstring("Most", "Most[expr] gives all but the last element of expr.");
    symtab_set_docstring("Rest", "Rest[expr] gives all but the first element of expr.");
    symtab_set_docstring("Append", "Append[expr, elem] adds elem to the end of expr.");
    symtab_set_docstring("Prepend", "Prepend[expr, elem] adds elem to the beginning of expr.");
    symtab_set_docstring("Insert", "Insert[expr, elem, n] inserts elem at position n in expr.");
    symtab_set_docstring("Delete", "Delete[expr, n] deletes the element at position n in expr.");
    symtab_set_docstring("Reverse", "Reverse[expr] reverses the order of elements in expr.");
    symtab_set_docstring("RotateLeft", "RotateLeft[expr, n] rotates the elements of expr n positions to the left.");
    symtab_set_docstring("RotateRight", "RotateRight[expr, n] rotates the elements of expr n positions to the right.");
    symtab_set_docstring("Transpose",
        "Transpose[list]\n"
        "\tTransposes the first two levels of list (swaps rows and columns of a matrix).\n"
        "Transpose[list, {n1, n2, ...}]\n"
        "\tGives the transpose of list so that level k in list is level nk in the result.\n"
        "\tThe spec must be a permutation of {1, ..., r} where r is the depth of list.\n"
        "\tA repeated index (e.g. {1, 1}) selects the corresponding diagonal.\n"
        "\tlist must be a rectangular array.");
    symtab_set_docstring("ConjugateTranspose",
        "ConjugateTranspose[m]\n"
        "\tGives the conjugate transpose of m, equivalent to Conjugate[Transpose[m]].\n"
        "ConjugateTranspose[m, spec]\n"
        "\tGives Conjugate[Transpose[m, spec]], permuting the levels of m according to\n"
        "\tthe spec list and then conjugating every entry.\n"
        "\tOn a 1-D vector, ConjugateTranspose[vec] conjugates the entries without\n"
        "\tchanging the shape of vec.");
    symtab_set_docstring("Flatten", "Flatten[list] flattens all levels of list.");
    symtab_set_docstring("Partition", "Partition[list, n] partitions list into sublists of length n.");

    // List Operations
    symtab_set_docstring("Table", "Table[expr, n]\n\tgenerates a list of n copies of expr.\nTable[expr, {i, imax}]\n\tgenerates a list of the values of expr with i running from 1 to imax.");
    symtab_set_docstring("Range", "Range[n]\n\tgenerates the list {1, 2, 3, ..., n}.\nRange[n, m]\n\tgenerates the list {n, n + 1, ..., m - 1, m}.\nRange[n, m, d]\n\tuses step d.");
    symtab_set_docstring("Array", "Array[f, n] generates a list of length n with elements f[1], f[2], ..., f[n].");
    symtab_set_docstring("Union", "Union[list] gives a sorted list of all distinct elements in list.");
    symtab_set_docstring("Tally", "Tally[list] counts the number of occurrences of each distinct element in list.");
    symtab_set_docstring("Commonest", "Commonest[list] gives a list of the elements that are the most common in list.\nCommonest[list, n] gives a list of the n most common elements in list.");
    symtab_set_docstring("DeleteDuplicates", "DeleteDuplicates[list] removes duplicate elements from list.");
    symtab_set_docstring("Split", "Split[list] splits list into sublists of identical adjacent elements.");

    // Statistics
    symtab_set_docstring("Mean", "Mean[data] gives the mean estimate of the elements in data.");
    symtab_set_docstring("RootMeanSquare", "RootMeanSquare[list] gives the root mean square of values in list.");
    symtab_set_docstring("Variance", "Variance[data] gives the unbiased variance estimate of the elements in data.");
    symtab_set_docstring("StandardDeviation", "StandardDeviation[data] gives the standard deviation estimate of the elements in data.");
    symtab_set_docstring("MovingAverage",
        "MovingAverage[list, r]\n"
        "\tgives the moving average of list, computed by averaging runs of r elements.\n"
        "MovingAverage[list, {w_1, w_2, ..., w_r}]\n"
        "\tgives the weighted moving average of list with weights w_i (effective weights w_i / Sum[w_i]).\n"
        "MovingAverage returns a list of length Length[list] - r + 1, and stays unevaluated when r < 1 or r > Length[list].");
    symtab_set_docstring("MovingMedian",
        "MovingMedian[list, r]\n"
        "\tgives the moving median of list, computed using spans of r elements.\n"
        "MovingMedian returns a list of length Length[list] - r + 1; for matrix input the medians are taken column-wise within each row-window.\n"
        "MovingMedian requires real-valued numeric data and stays unevaluated when r < 1 or r > Length[list].");
    symtab_set_docstring("ExponentialMovingAverage",
        "ExponentialMovingAverage[list, alpha]\n"
        "\tgives the exponential moving average of list with smoothing constant alpha.\n"
        "Defined by the recurrence y_1 = x_1, y_{i+1} = y_i + alpha (x_{i+1} - y_i).\n"
        "The output has the same length as list. The smoothing constant alpha is typically a number between 0 and 1, but may be any expression; ExponentialMovingAverage handles both numerical (machine and arbitrary precision) and symbolic data.");

    // Functional Programming
    symtab_set_docstring("Map", "f /@ expr or Map[f, expr] applies f to each element of expr.");
    symtab_set_docstring("Apply", "f @@ expr or Apply[f, expr] replaces the head of expr with f.");
    symtab_set_docstring("MapAll", "f //@ expr or MapAll[f, expr] applies f to every subexpression in expr.");
    symtab_set_docstring("Through", "Through[p[f, g][x]] gives p[f[x], g[x]].");
    symtab_set_docstring("Thread",
        "Thread[f[args]]\n"
        "\t\"threads\" f over any lists that appear in args.\n"
        "Thread[f[args], h]\n"
        "\tthreads f over any objects with head h that appear in args.\n"
        "Thread[f[args], h, n]\n"
        "\tthreads f over objects with head h that appear in the first n args.\n"
        "\n"
        "Functions with attribute Listable are automatically threaded over\n"
        "lists. All the elements in the specified args whose heads are h must\n"
        "be of the same length. Arguments that do not have head h are copied\n"
        "as many times as there are elements in the arguments that do have\n"
        "head h.\n"
        "\n"
        "Thread specifies argument positions using the standard sequence\n"
        "specification:\n"
        "\tAll       all elements\n"
        "\tNone      no elements\n"
        "\tn         elements 1 through n\n"
        "\t-n        last n elements\n"
        "\t{n}       element n only\n"
        "\t{m, n}    elements m through n inclusive\n"
        "\t{m, n, s} elements m through n in steps of s\n"
        "\n"
        "Examples:\n"
        "\tThread[f[{a, b, c}]]              -> {f[a], f[b], f[c]}\n"
        "\tThread[f[{a, b, c}, x]]           -> {f[a, x], f[b, x], f[c, x]}\n"
        "\tThread[f[{a, b, c}, {x, y, z}]]   -> {f[a, x], f[b, y], f[c, z]}\n"
        "\tThread[{a, b, c} == {x, y, z}]    -> {a == x, b == y, c == z}\n"
        "\tThread[Log[x == y], Equal]        -> Log[x] == Log[y]");
    symtab_set_docstring("Select", "Select[list, crit] selects elements of list that satisfy crit.");
    symtab_set_docstring("FreeQ", "FreeQ[expr, form] yields True if no subexpression in expr matches form, and yields False otherwise.");
    symtab_set_docstring("Function",
        "body & or Function[body]\n"
        "\trepresents a pure function with formal parameters #, #1, #2, ... and ##, ##1, ##2, ... for sequences of arguments.\n"
        "Function[x, body] or Function[{x1, x2, ...}, body]\n"
        "\trepresents a pure function with named formal parameters x or x1, x2, ....\n"
        "Function[params, body, attrs]\n"
        "\tis a pure function that is treated as having attributes attrs for purposes of evaluation.\n"
        "\tattrs can be a single attribute or a list of attributes; recognised attributes include HoldFirst, HoldRest, HoldAll, HoldAllComplete, Listable, Flat, Orderless, OneIdentity, NumericFunction, SequenceHold, and NHoldRest.\n"
        "Function[Null, body, attrs]\n"
        "\trepresents a function in which the parameters in body are given using # etc.\n"
        "\n"
        "Parameter binding is lexical: named parameters are substituted into the body before evaluation. Nested Function expressions shadow their own parameters.\n"
        "By default Function has no Hold attributes; the arguments are evaluated before substitution. Adding HoldAll (or HoldFirst / HoldRest / HoldAllComplete) in the 3-arg form holds arguments in the chosen positions.");
    symtab_set_docstring("Slot", "# or Slot[n] represents the n-th argument of a pure function.");
    symtab_set_docstring("SlotSequence", "## or SlotSequence[n] represents arguments from the n-th onward.");

    // Predicates
    symtab_set_docstring("AtomQ", "AtomQ[expr] gives True if expr is an atomic object.");
    symtab_set_docstring("Identity", "Identity[expr] gives expr (the identity operation).");
    symtab_set_docstring("Composition",
        "Composition[f1, f2, f3, ...]\n"
        "\trepresents a composition of the functions f1, f2, f3, ....\n"
        "\n"
        "Composition allows you to build up compositions of functions which can\n"
        "later be applied to specific arguments. Applied to arguments, the\n"
        "composition acts innermost-first:\n"
        "\tComposition[f, g, h][x, y]  ->  f[g[h[x, y]]].\n"
        "\n"
        "Composition has the attributes Flat and OneIdentity.\n"
        "Composition can be entered in the form f1 @* f2 @* ....\n"
        "\n"
        "Composition objects containing Identity or InverseFunction[f] are\n"
        "automatically simplified when possible:\n"
        "\tComposition[]                       ->  Identity\n"
        "\tComposition[f]                      ->  f\n"
        "\tComposition[f, Identity, g]         ->  Composition[f, g]\n"
        "\tComposition[f, InverseFunction[f]]  ->  Identity.");
    symtab_set_docstring("ComposeList",
        "ComposeList[{f1, f2, ...}, x]\n"
        "\tgenerates a list of the form {x, f1[x], f2[f1[x]], ...}.\n"
        "\n"
        "ComposeList applies its functions innermost-first and accumulates\n"
        "the intermediate results. The output list has one more element than\n"
        "the input list of functions. Function applications are evaluated\n"
        "in the normal way after construction:\n"
        "\tComposeList[{a, b, c}, x]  ->  {x, a[x], b[a[x]], c[b[a[x]]]}.\n"
        "\n"
        "ComposeList has the attribute Protected.");
    symtab_set_docstring("Accumulate",
        "Accumulate[list]\n"
        "\tgives a list of the successive accumulated totals of elements in\n"
        "\tlist. The result has the same length as list.\n"
        "\n"
        "Accumulate[list] is effectively equivalent to FoldList[Plus, list].\n"
        "Accumulate works with integers, arbitrary-precision bignums, machine\n"
        "doubles, and symbolic expressions, and threads naturally over rows\n"
        "(so for a matrix it accumulates within columns). The head of the\n"
        "input is preserved:\n"
        "\tAccumulate[{a, b, c, d}]    ->  {a, a + b, a + b + c, a + b + c + d}\n"
        "\tAccumulate[f[a, b, c, d]]   ->  f[a, a + b, a + b + c, a + b + c + d]\n"
        "\n"
        "Accumulate[list, Method -> \"CompensatedSummation\"] uses Kahan\n"
        "compensated summation to reduce numerical error when every element\n"
        "of list is a machine number. For symbolic or mixed input the option\n"
        "is ignored and the standard symbolic accumulation is returned.\n"
        "\n"
        "Accumulate has the attribute Protected.");
    symtab_set_docstring("NumberQ", "NumberQ[expr] gives True if expr is a number.");
    symtab_set_docstring("NumericQ", "NumericQ[expr] gives True if expr is a numeric quantity, and False otherwise.\nAn expression is considered a numeric quantity if it is either an explicit number or a mathematical constant such as Pi, or is a function that has attribute NumericFunction and all of whose arguments are numeric quantities.");
    symtab_set_docstring("IntegerQ", "IntegerQ[expr] gives True if expr is an integer.");
    symtab_set_docstring("EvenQ", "EvenQ[n] gives True if n is an even integer.");
    symtab_set_docstring("OddQ", "OddQ[n] gives True if n is an odd integer.");
    symtab_set_docstring("PrimeQ", "PrimeQ[n]\n\tgives True if n is a prime number.\nPrimeQ[z]\n\tfor a Gaussian integer z = a + b I, gives True if z is a Gaussian prime.");
    symtab_set_docstring("PolynomialQ", "PolynomialQ[expr, var] gives True if expr is a polynomial in var.");
    symtab_set_docstring("ListQ", "ListQ[expr] gives True if expr is a list.");
    symtab_set_docstring("VectorQ",
        "VectorQ[expr]\n"
        "\tgives True if expr is a list, none of whose elements are themselves lists, and gives False otherwise.\n"
        "VectorQ[expr, test]\n"
        "\tgives True only if test yields True when applied to each of the elements in expr.\n"
        "\n"
        "VectorQ[expr, NumberQ] tests whether expr is a vector of numbers.");
    symtab_set_docstring("MatrixQ",
        "MatrixQ[expr]\n"
        "\tgives True if expr is a list of lists that can represent a matrix, and gives False otherwise.\n"
        "MatrixQ[expr, test]\n"
        "\tgives True only if test yields True when applied to each of the matrix elements in expr.\n"
        "\n"
        "MatrixQ[expr] gives True only if expr is a list and each of its elements is a list of the same length,\n"
        "containing no elements that are themselves lists.\n"
        "MatrixQ[expr, NumberQ] tests whether expr is a numerical matrix.");
    symtab_set_docstring("MatchQ", "MatchQ[expr, form] gives True if expr matches form.");

    // Trigonometric
    symtab_set_docstring("Sin", "Sin[z] gives the sine of z.");
    symtab_set_docstring("Cos", "Cos[z] gives the cosine of z.");
    symtab_set_docstring("Tan", "Tan[z] gives the tangent of z.");
    symtab_set_docstring("Cot", "Cot[z] gives the cotangent of z.");
    symtab_set_docstring("Sec", "Sec[z] gives the secant of z.");
    symtab_set_docstring("Csc", "Csc[z] gives the cosecant of z.");
    symtab_set_docstring("ArcSin", "ArcSin[z] gives the inverse sine of z.");
    symtab_set_docstring("ArcCos", "ArcCos[z] gives the inverse cosine of z.");
    symtab_set_docstring("ArcTan", "ArcTan[z] gives the inverse tangent of z.");

    // Hyperbolic
    symtab_set_docstring("Sinh", "Sinh[z] gives the hyperbolic sine of z.");
    symtab_set_docstring("Cosh", "Cosh[z] gives the hyperbolic cosine of z.");
    symtab_set_docstring("Tanh", "Tanh[z] gives the hyperbolic tangent of z.");

    // Log/Exp
    symtab_set_docstring("Log", "Log[z] gives the natural logarithm of z. Log[b, z] gives the logarithm to base b.");
    symtab_set_docstring("Exp", "Exp[z] gives the exponential of z.");

    // Trig simplification / conversion
    symtab_set_docstring("TrigToExp",
        "TrigToExp[expr]\n\trewrites circular and hyperbolic trigonometric functions (and their\n\tinverses) in expr in terms of exponentials and logarithms.");
    symtab_set_docstring("ExpToTrig",
        "ExpToTrig[expr]\n\trewrites exponentials and logarithms in expr in terms of circular and\n\thyperbolic trigonometric functions when possible.");
    symtab_set_docstring("TrigExpand",
        "TrigExpand[expr]\n\texpands out trigonometric functions in expr.\n\tTrigExpand operates on both circular and hyperbolic functions.\n\tTrigExpand splits up sums and integer multiples that appear in arguments of\n\ttrigonometric functions, and then expands out products of trigonometric\n\tfunctions into sums of powers, using trigonometric identities when possible.\n\tTrigExpand automatically threads over lists, as well as equations,\n\tinequalities, and logic functions.");
    symtab_set_docstring("TrigFactor",
        "TrigFactor[expr]\n\tfactors trigonometric functions in expr.\n\tTrigFactor operates on both circular and hyperbolic functions.\n\tTrigFactor factors polynomials in trigonometric functions and collapses\n\tPythagorean, angle-addition, and double-angle identities where possible,\n\tbroadly acting as the inverse of TrigExpand.\n\tTrigFactor automatically threads over lists, as well as equations,\n\tinequalities, and logic functions.");
    symtab_set_docstring("TrigReduce",
        "TrigReduce[expr]\n"
        "\trewrites products and powers of trigonometric functions in expr in terms\n"
        "\tof trigonometric functions with combined arguments.\n"
        "\tTrigReduce operates on both circular and hyperbolic functions.\n"
        "\tGiven a trigonometric polynomial, TrigReduce typically yields a linear\n"
        "\texpression involving trigonometric functions with more complicated\n"
        "\targuments.\n"
        "\tTrigReduce automatically threads over lists, as well as equations,\n"
        "\tinequalities and logic functions.\n"
        "\n"
        "\tReduce trigonometric expressions:\n"
        "\t  TrigReduce[2 Cos[x]^2]            ->  1 + Cos[2 x]\n"
        "\t  TrigReduce[2 Sin[x] Cos[y]]       ->  Sin[x - y] + Sin[x + y]\n"
        "\n"
        "\tReduce hyperbolic trigonometric expressions:\n"
        "\t  TrigReduce[2 Cosh[x]^2]           ->  1 + Cosh[2 x]\n"
        "\t  TrigReduce[2 Sinh[x] Cosh[y]]     ->  Sinh[x - y] + Sinh[x + y]\n"
        "\n"
        "\tTrigonometric expressions:\n"
        "\t  TrigReduce[2 Sin[x + y] Cos[x - y]]  ->  Sin[2 x] + Sin[2 y]\n"
        "\t  TrigReduce[Tan[x] + Tan[y]]          ->  Sec[x] Sec[y] Sin[x + y]\n"
        "\n"
        "\tHyperbolic trigonometric expressions:\n"
        "\t  TrigReduce[2 Cosh[x] Cosh[y]]     ->  Cosh[x - y] + Cosh[x + y]\n"
        "\t  TrigReduce[Coth[x] + Coth[y]]     ->  Csch[x] Csch[y] Sinh[x + y]\n"
        "\n"
        "\tThreads over lists:\n"
        "\t  TrigReduce[{Tan[x] + Cot[y], Tanh[x] - Coth[y]}]\n"
        "\t    -> {Cos[x - y] Csc[y] Sec[x], -Cosh[x - y] Csch[y] Sech[x]}\n"
        "\n"
        "\tThreads over equations, inequalities, and logical operations:\n"
        "\t  TrigReduce[4 Sin[x]^4 == 1 && 2 Cos[x]^2 >= 1]\n"
        "\t    -> 1/2 (3 - 4 Cos[2 x] + Cos[4 x]) == 1 && 1 + Cos[2 x] >= 1");

    // Piecewise / Rounding
    symtab_set_docstring("Floor", "Floor[x] gives the greatest integer less than or equal to x.");
    symtab_set_docstring("Ceiling", "Ceiling[x] gives the smallest integer greater than or equal to x.");
    symtab_set_docstring("Round", "Round[x] rounds x to the nearest integer.");
    symtab_set_docstring("Chop",
        "Chop[expr]\n"
        "\treplaces approximate real numbers in expr that are close to zero\n"
        "\tby the exact integer 0.\n"
        "Chop[expr, delta]\n"
        "\treplaces numbers smaller in absolute magnitude than delta by 0.\n"
        "\n"
        "The default tolerance is 10^-10. Chop walks the entire expression\n"
        "tree, so small real-valued subterms inside arbitrary heads, lists,\n"
        "and held forms are all chopped. Exact numbers -- integers, bigints,\n"
        "rationals, and symbolic constants -- pass through untouched.\n"
        "\n"
        "For machine complex numbers Complex[re, im] whose real and imaginary\n"
        "parts are both machine reals: if only the imaginary part is below\n"
        "tolerance the whole Complex wrapper is dropped and the real part\n"
        "is returned; if only the real part is below tolerance the result\n"
        "is Complex[0., im], preserving the machine-complex shape with a\n"
        "machine zero. If both parts are below tolerance the result is the\n"
        "exact integer 0.");
    symtab_set_docstring("Clip",
        "Clip[x]\n"
        "\tgives x clipped to be between -1 and +1.\n"
        "Clip[x, {min, max}]\n"
        "\tgives x for min <= x <= max, min for x < min, and max for x > max.\n"
        "Clip[x, {min, max}, {vmin, vmax}]\n"
        "\tgives vmin for x < min and vmax for x > max.\n"
        "\n"
        "Clip threads over lists in its first argument and works at machine\n"
        "or arbitrary precision (via N). Symbolic constants such as Pi are\n"
        "numericalized only to decide which side of the interval x lies on;\n"
        "the original symbolic x is returned unchanged when min <= x <= max.\n"
        "\n"
        "Infinity and -Infinity are clipped to the upper and lower\n"
        "replacement values respectively. Clip is not defined for non-real\n"
        "complex values: Clip::ncompl is issued and the call is returned\n"
        "unevaluated. Clip[a] for an otherwise undetermined a also stays\n"
        "unevaluated so user-supplied rules can intercept it.");

    // Calculus
    symtab_set_docstring("D", "D[f, x] gives the partial derivative of f with respect to x.\nD[f, {x, n}] gives the nth partial derivative.\nD[f, x, y, ...] gives the mixed derivative.");
    symtab_set_docstring("Dt", "Dt[f] gives the total derivative of f.\nDt[f, x] gives the total derivative of f with respect to x.\nDt[f, {x, n}] gives the nth total derivative.");
    symtab_set_docstring("Derivative", "Derivative[n][f][x] represents the nth derivative of a function f evaluated at x.");

    // Control Flow
    symtab_set_docstring("Do",
        "Do[expr, n] evaluates expr n times.\n"
        "Do[expr, {i, imax}] evaluates expr with i successively taking on values 1 through imax.\n"
        "Do[expr, {i, imin, imax}] starts with i = imin.\n"
        "Do[expr, {i, imin, imax, di}] uses steps di.\n"
        "Do[expr, {i, {i1, i2, ...}}] uses the successive values i1, i2, ....\n"
        "Do[expr, {n}] evaluates expr n times with no iteration variable.\n"
        "Do[expr, iter1, iter2, ...] iterates over multiple variables, with the rightmost varying fastest.\n"
        "Do has attribute HoldAll: expr and the iterator specifications are held unevaluated until each iteration.\n"
        "Break[] inside expr exits the innermost Do loop.\n"
        "Continue[] inside expr skips the rest of expr and proceeds to the next iteration.\n"
        "Return[v] inside expr causes the enclosing function to yield v; Do itself returns Null.");
    symtab_set_docstring("For",
        "For[start, test, incr, body] executes start, then repeatedly evaluates body and incr until test fails to give True.\n"
        "For[start, test, incr] runs with an empty body, useful when incr or test carries the side-effect.\n"
        "For has attribute HoldAll: start, test, incr, and body are all held unevaluated until For drives them.\n"
        "Break[] inside body exits the loop.\n"
        "Continue[] inside body skips the rest of body and re-evaluates incr then test.\n"
        "Return[v] inside body causes the enclosing function to yield v; For itself returns Null.");
    symtab_set_docstring("While",
        "While[test, body] evaluates test, then body, repeatedly, until test first fails to give True.\n"
        "While[test] does the loop with a Null body, which is useful when test has side-effects.\n"
        "While has attribute HoldAll.\n"
        "Break[] inside body exits the loop.\n"
        "Continue[] inside body skips the rest of body and re-evaluates test.\n"
        "Return[v] inside body causes While to yield v; otherwise While returns Null.");
    symtab_set_docstring("If", "If[condition, t, f] gives t if condition evaluates to True, and f if it evaluates to False.\nIf[condition, t, f, u] gives u if condition evaluates to neither True nor False.");
    symtab_set_docstring("Which",
        "Which[test1, value1, test2, value2, ...]\n"
        "\tevaluates each test_i in turn, returning the corresponding value_i\n"
        "\tfor the first test that yields True.\n"
        "Which has attribute HoldAll, so the tests and values are not\n"
        "evaluated until Which examines them.\n"
        "If a test evaluates to neither True nor False, a Which object\n"
        "containing that test (in evaluated form) and the remaining\n"
        "elements is returned unevaluated.\n"
        "If all tests evaluate to False (or no tests are supplied), Which\n"
        "returns Null.\n"
        "Use True as the final test to supply a default value.");
    symtab_set_docstring("Switch",
        "Switch[expr, form_1, value_1, form_2, value_2, ...]\n"
        "\tevaluates expr, then compares it with each of the form_i in turn,\n"
        "\tevaluating and returning the value_i corresponding to the first\n"
        "\tmatch found.\n"
        "Only the value_i corresponding to the first form_i that matches\n"
        "expr is evaluated. Each form_i is itself evaluated only when the\n"
        "match is tried.\n"
        "If the last form_i is the pattern _, then the corresponding\n"
        "value_i is always returned if this case is reached -- it acts as\n"
        "a default branch.\n"
        "If none of the form_i match expr, the Switch is returned\n"
        "unevaluated.\n"
        "Switch has attribute HoldRest, so the form/value pairs are not\n"
        "evaluated until Switch examines them.\n"
        "Break, Return, and Throw inside the chosen value behave as they\n"
        "do in any other held context.");
    symtab_set_docstring("Piecewise",
        "Piecewise[{{val_1, cond_1}, {val_2, cond_2}, ...}]\n"
        "\trepresents a piecewise function with values val_i in the regions\n"
        "\tdefined by the conditions cond_i.\n"
        "Piecewise[{{val_1, cond_1}, ...}, val]\n"
        "\tuses the default value val if none of the cond_i apply. The\n"
        "\tdefault for val is 0.\n"
        "\n"
        "The cond_i are evaluated in turn until one yields True. If all\n"
        "preceding cond_i yield False, the corresponding val_i of the\n"
        "first True cond_i is returned. If any preceding cond_i does not\n"
        "literally yield False, the Piecewise expression is returned in\n"
        "symbolic form.\n"
        "\n"
        "Only those val_i explicitly included in the returned form are\n"
        "evaluated (Piecewise has attribute HoldAll). Pairs of the form\n"
        "{val_i, False} are dropped, and all clauses after the first\n"
        "{val_i, True} are dropped together with the default value.\n"
        "Consecutive clauses with structurally identical values are\n"
        "merged: their conditions are combined with Or.\n"
        "\n"
        "Piecewise[conds] automatically evaluates to Piecewise[conds, 0].");
    symtab_set_docstring("TrueQ", "TrueQ[expr] yields True if expr is True, and False otherwise.");
    symtab_set_docstring("Boole",
        "Boole[expr]\n"
        "\tyields 1 if expr is True and 0 if it is False.\n"
        "Boole is also known as the Iverson bracket, indicator function, and characteristic function.\n"
        "Boole is typically used to express integrals and sums over regions given by logical combinations of predicates, and as a dummy-variable encoding for categorical variables in statistics.\n"
        "Boole[expr] remains unchanged if expr is neither True nor False.\n"
        "Boole[expr] is effectively equivalent to If[expr, 1, 0].\n"
        "Boole has attributes {Listable, Protected}, so Boole[{e1, e2, ...}] automatically threads to {Boole[e1], Boole[e2], ...}.");
    symtab_set_docstring("Evaluate", "Evaluate[expr]\n\tcauses expr to be evaluated even if it appears as the argument of a function whose attributes specify that it should be held unevaluated.\nEvaluate only overrides HoldFirst, HoldRest, and HoldAll attributes when it appears directly as the head of the function argument that would otherwise be held.\nEvaluate does not override HoldAllComplete.");
    symtab_set_docstring("ReleaseHold", "ReleaseHold[expr]\n\tremoves Hold, HoldForm, HoldPattern, and HoldComplete in expr.\nReleaseHold removes only one layer of Hold etc.; it does not remove inner occurrences in nested Hold etc. functions.");
    symtab_set_docstring("HoldPattern",
        "HoldPattern[expr]\n"
        "\tis equivalent to expr for pattern matching, but maintains expr in an unevaluated form.\n"
        "HoldPattern has attributes {HoldAll, Protected}.\n"
        "The left-hand sides of rules and assignments are normally evaluated before being used for matching; wrap the LHS in HoldPattern to stop that evaluation (e.g. HoldPattern[_+_] -> 0 matches any two-term sum, whereas _+_ -> 0 would match a pattern simplified by Plus before the rule is applied).\n"
        "HoldPattern is removed by one layer of ReleaseHold.");
    symtab_set_docstring("Unevaluated",
        "Unevaluated[expr]\n"
        "\trepresents the unevaluated form of expr when it appears as the argument to a function.\n"
        "f[Unevaluated[expr]] effectively works by temporarily holding that argument, then evaluating f[expr] with the wrapper removed.\n"
        "The wrapper is NOT removed when the argument is held (e.g. when f has HoldAll, HoldFirst, or HoldRest applies) or when f has attribute HoldAllComplete.\n"
        "Sequence expressions directly inside Unevaluated are NOT flattened: Length[Unevaluated[Sequence[a, b]]] gives 2.\n"
        "Unevaluated has attributes {HoldAllComplete, Protected}, so its own argument is itself protected from evaluation, Sequence flattening, and wrapper stripping.");
    symtab_set_docstring("Hold",
        "Hold[expr]\n"
        "\tmaintains expr in an unevaluated form.\n"
        "Hold has attribute HoldAll: its arguments are not evaluated.\n"
        "Evaluate[expr] inside Hold overrides the hold and evaluates expr once.\n"
        "Sequence expressions inside Hold are flattened; use HoldComplete to prevent this.");
    symtab_set_docstring("HoldComplete",
        "HoldComplete[expr]\n"
        "\tshields expr completely from evaluation.\n"
        "HoldComplete has attribute HoldAllComplete: it prevents argument evaluation, Sequence flattening, Unevaluated stripping, and Evaluate from firing.\n"
        "Substitution (via ReplaceAll, etc.) still happens inside HoldComplete.\n"
        "HoldComplete is removed by one level of ReleaseHold.");
    symtab_set_docstring("HoldAllComplete",
        "HoldAllComplete\n"
        "\tis an attribute which specifies that all arguments to a function are not to be modified or looked at in any way in the process of evaluation.\n"
        "HoldAllComplete prevents argument evaluation, Sequence flattening inside arguments, Unevaluated wrapper stripping, and application of Evaluate.\n"
        "Evaluate cannot override HoldAllComplete.");
    symtab_set_docstring("Set", "lhs = rhs assigns rhs to lhs.");

    // In-place numeric assignment operators
    symtab_set_docstring("Increment",
        "Increment[x] or x++\n"
        "\tincreases the value of x by 1, returning the old value of x.\n"
        "\n"
        "Increment has attribute HoldFirst. In Increment[x], x can be a symbol\n"
        "or a Part expression referring to an existing value (e.g. list[[2]]++).\n"
        "Increment threads over list values because Plus is Listable.\n"
        "If x has no assigned value, Increment::rvalue is emitted and the\n"
        "expression is left unevaluated.");
    symtab_set_docstring("Decrement",
        "Decrement[x] or x--\n"
        "\tdecreases the value of x by 1, returning the old value of x.\n"
        "\n"
        "Decrement has attribute HoldFirst. In Decrement[x], x can be a symbol\n"
        "or a Part expression referring to an existing value (e.g. list[[2]]--).\n"
        "If x has no assigned value, Decrement::rvalue is emitted and the\n"
        "expression is left unevaluated.");
    symtab_set_docstring("PreIncrement",
        "PreIncrement[x] or ++x\n"
        "\tincreases the value of x by 1, returning the new value of x.\n"
        "\t++x is equivalent to x = x + 1.\n"
        "\n"
        "PreIncrement has attribute HoldFirst. In PreIncrement[x], x can be a\n"
        "symbol or a Part expression referring to an existing value.\n"
        "If x has no assigned value, PreIncrement::rvalue is emitted and the\n"
        "expression is left unevaluated.");
    symtab_set_docstring("PreDecrement",
        "PreDecrement[x] or --x\n"
        "\tdecreases the value of x by 1, returning the new value of x.\n"
        "\t--x is equivalent to x = x - 1.\n"
        "\n"
        "PreDecrement has attribute HoldFirst. In PreDecrement[x], x can be a\n"
        "symbol or a Part expression referring to an existing value.\n"
        "If x has no assigned value, PreDecrement::rvalue is emitted and the\n"
        "expression is left unevaluated.");
    symtab_set_docstring("AddTo",
        "AddTo[x, dx] or x += dx\n"
        "\tadds dx to x and returns the new value of x.\n"
        "\tx += dx is equivalent to x = x + dx.\n"
        "\n"
        "AddTo has attribute HoldFirst. The first argument x can be a symbol or\n"
        "a Part expression referring to an existing value; dx may be a number,\n"
        "a symbolic expression, or a list (combined element-wise via the Listable\n"
        "attribute of Plus). If x has no assigned value, AddTo::rvalue is emitted\n"
        "and the expression is left unevaluated.");
    symtab_set_docstring("SubtractFrom",
        "SubtractFrom[x, dx] or x -= dx\n"
        "\tsubtracts dx from x and returns the new value of x.\n"
        "\tx -= dx is equivalent to x = x - dx.\n"
        "\n"
        "SubtractFrom has attribute HoldFirst. The first argument x can be a\n"
        "symbol or a Part expression referring to an existing value; dx may be a\n"
        "number, a symbolic expression, or a list. If x has no assigned value,\n"
        "SubtractFrom::rvalue is emitted and the expression is left unevaluated.");
    symtab_set_docstring("SetDelayed", "lhs := rhs assigns rhs to lhs, evaluating it only when needed.");
    symtab_set_docstring("Default", "Default[f] gives the default value for arguments of the function f obtained with a _. pattern object.");
    symtab_set_docstring("Optional", "patt:def or Optional[patt,def] is a pattern object that represents an expression of the form patt, which, if omitted, should be replaced by the default value def.");
    symtab_set_docstring("Longest", "Longest[p] is a pattern object that matches the longest sequence consistent with the pattern p.");
    symtab_set_docstring("Shortest", "Shortest[p] is a pattern object that matches the shortest sequence consistent with the pattern p.");
    symtab_set_docstring("Repeated", "p.. or Repeated[p] is a pattern object that represents a sequence of one or more expressions, each matching p.\nRepeated[p, max] represents from 1 to max expressions matching p.\nRepeated[p, {min, max}] represents between min and max expressions matching p.\nRepeated[p, {n}] represents exactly n expressions matching p.");
    symtab_set_docstring("RepeatedNull", "p... or RepeatedNull[p] is a pattern object that represents a sequence of zero or more expressions, each matching p.\nRepeatedNull[p, max] represents from 0 to max expressions matching p.\nRepeatedNull[p, {min, max}] represents between min and max expressions matching p.\nRepeatedNull[p, {n}] represents exactly n expressions matching p.");
    symtab_set_docstring("Blank", "_ or Blank[] represents any single expression.\n_h or Blank[h] represents any single expression with head h.");
    symtab_set_docstring("BlankSequence", "__ or BlankSequence[] represents a sequence of one or more expressions.");
    symtab_set_docstring("BlankNullSequence", "___ or BlankNullSequence[] represents a sequence of zero or more expressions.");
    symtab_set_docstring("Clear", "Clear[x, y, ...] clears the values of symbols.");
    symtab_set_docstring("Flat", "Flat is an attribute that can be assigned to a symbol f to indicate that all expressions involving nested functions f should be flattened out. This property is accounted for in pattern matching.");
    symtab_set_docstring("Orderless", "Orderless is an attribute that can be assigned to a symbol f to indicate that the elements e_i in expressions of the form f[e_1, e_2, ...] should automatically be sorted into canonical order. This property is accounted for in pattern matching.");
    symtab_set_docstring("OneIdentity", "OneIdentity is an attribute that can be assigned to a symbol f to indicate that f[x], f[f[x]], etc. are all equivalent to x for the purpose of pattern matching.");
    symtab_set_docstring("Information", "Information[symbol] or ?symbol returns information on symbol.");
    symtab_set_docstring("OwnValues", "OwnValues[s] gives a list of own-value rules for s.");
    symtab_set_docstring("DownValues", "DownValues[s] gives a list of down-value rules for s.");
    symtab_set_docstring("Attributes", "Attributes[s] gives the list of attributes for s.");
    symtab_set_docstring("SetAttributes", "SetAttributes[s, attr] sets the attributes for s.");
    symtab_set_docstring("ClearAttributes",
        "ClearAttributes[s, attr] removes attr from the list of attributes of s.\n"
        "ClearAttributes[s, {attr1, attr2, ...}] removes several attributes at a time.\n"
        "ClearAttributes[{s1, s2, ...}, attrs] removes attributes from several symbols at a time.");
    symtab_set_docstring("CompoundExpression", "expr1; expr2; ... evaluates its arguments in sequence, returning the last result.");
    symtab_set_docstring("Module", "Module[{x, y, ...}, expr] specifies that x, y, ... are local variables.");
    symtab_set_docstring("Block", "Block[{x, y, ...}, expr] evaluates expr with local values for x, y, ....");
    symtab_set_docstring("With", "With[{x = x0, ...}, expr] specifies that x should be replaced by x0 throughout expr.");
    symtab_set_docstring("Return",
        "Return[expr]\n"
        "\treturns the value expr from a function.\n"
        "Return[]\n"
        "\treturns the value Null.\n"
        "Return[expr, h]\n"
        "\treturns expr from the nearest enclosing construct whose head is h.\n"
        "\n"
        "Return[expr] exits control structures within the definition of a function,\n"
        "and gives the value expr for the whole function.\n"
        "Return takes effect as soon as it is evaluated, even if it appears inside\n"
        "other functions.\n"
        "\n"
        "The recognised boundary heads are Function, Module, Block, With, Do, For,\n"
        "and While. CompoundExpression and other Hold-free heads pass Return through\n"
        "unchanged so that it can reach the boundary.");

    // Rules and Replacements
    symtab_set_docstring("Rule", "lhs -> rhs represents a rule that transforms lhs to rhs.");
    symtab_set_docstring("RuleDelayed", "lhs :> rhs represents a rule that transforms lhs to rhs, evaluating rhs only when the rule is used.");
    symtab_set_docstring("Replace", "Replace[expr, rules] applies rules to the entire expr. Replace[expr, rules, levelspec] applies rules to parts of expr specified by levelspec.");
    symtab_set_docstring("ReplaceAll", "expr /. rules or ReplaceAll[expr, rules] applies rules to transform each subpart of expr.");
    symtab_set_docstring("ReplaceRepeated", "expr //. rules or ReplaceRepeated[expr, rules] repeatedly applies rules until the expression no longer changes.");
    symtab_set_docstring("Print", "Print[expr1, expr2, ...] prints the expressions to stdout and returns Null.");

    // File I/O
    symtab_set_docstring("Get",
        "Get[\"filename\"]\n"
        "\treads expressions from a file, evaluates them in order, and returns the last result.\n"
        "Returns $Failed if the file cannot be opened.\n"
        "It is conventional to use names ending in .m for files containing Mathilda input.");
    symtab_set_docstring("Put",
        "Put[expr, \"filename\"] or expr >> \"filename\"\n"
        "\twrites expr to the file, replacing any prior contents.\n"
        "Put[expr1, expr2, ..., \"filename\"]\n"
        "\twrites a sequence of expressions to the file, each followed by a newline.\n"
        "Put[\"filename\"]\n"
        "\tcreates an empty file with the specified name (or truncates an existing one).\n"
        "expr >> \"filename\" is equivalent to Put[expr, \"filename\"]; the bare-word form\n"
        "expr >> filename is equivalent to expr >> \"filename\".\n"
        "Returns Null on success and $Failed if the file cannot be opened.");
    symtab_set_docstring("PutAppend",
        "PutAppend[expr, \"filename\"] or expr >>> \"filename\"\n"
        "\tappends expr to the end of the file, creating the file if it does not exist.\n"
        "PutAppend[expr1, expr2, ..., \"filename\"]\n"
        "\tappends a sequence of expressions, one per line.\n"
        "PutAppend works the same as Put, except that it preserves any existing\n"
        "contents of the file rather than truncating them.\n"
        "expr >>> filename is equivalent to expr >>> \"filename\".\n"
        "Returns Null on success and $Failed if the file cannot be opened.");
    symtab_set_docstring("FileExistsQ",
        "FileExistsQ[\"name\"]\n"
        "\tgives True if the file with the specified name exists, and gives False otherwise.\n"
        "In FileExistsQ[\"name\"], name is interpreted relative to your current directory.\n"
        "FileExistsQ does not search $Path.\n"
        "FileExistsQ tests for files, directories, or any other filesystem objects.");
    symtab_set_docstring("FileExtension",
        "FileExtension[\"file\"]\n"
        "\tgives the file extension for a file name.\n"
        "FileExtension[\"name.ext\"] gives \"ext\".\n"
        "FileExtension gives the extension that appears after the last . in a file name.\n"
        "If there are multiple endings to a file name, separated by ., FileExtension gives only the last one.\n"
        "FileExtension gives \"\" if there is no file extension, if the file name has the form of a directory name, or ends with a . character.\n"
        "FileExtension ignores any directory specification.\n"
        "FileExtension by default assumes pathname separators and other conventions suitable for your operating system.");
    symtab_set_docstring("FileBaseName",
        "FileBaseName[\"file\"]\n"
        "\tgives the base name for a file without its extension.\n"
        "FileBaseName[\"name.ext\"] gives \"name\".\n"
        "FileBaseName drops all directory specifications.\n"
        "FileBaseName[\"file.tar.gz\"] gives \"file.tar\" — only the final extension is split off.\n"
        "FileBaseName by default assumes pathname separators and other conventions suitable for your operating system.");
    symtab_set_docstring("FilePrint",
        "FilePrint[\"file\"]\n"
        "\tprints out the raw textual contents of file.\n"
        "FilePrint[\"file\", n]\n"
        "\tprints out the first n raw textual lines of file.\n"
        "FilePrint[\"file\", -n]\n"
        "\tprints out the last n raw textual lines of file.\n"
        "FilePrint[\"file\", m;;n]\n"
        "\tprints out lines m through n of file.\n"
        "FilePrint[\"file\", m;;n;;s]\n"
        "\tprints out lines m through n of file in steps of s; s may be negative\n"
        "\tto walk the range backwards.\n"
        "FilePrint returns Null on success and $Failed if the file cannot be opened.\n"
        "Negative indices inside the Span count from the end of the file (-1 is the last line).");

    symtab_set_docstring("FullForm", "FullForm[expr] is a wrapper that causes expr to be printed in full form.");
    symtab_set_docstring("InputForm", "InputForm[expr] is a wrapper that causes expr to be printed in input form.");
    symtab_set_docstring("TeXForm",
        "TeXForm[expr]\n"
        "\tprints as a TeX version of expr.\n"
        "TeXForm produces AMS-LaTeX-compatible TeX output.\n"
        "When an input evaluates to TeXForm[expr], TeXForm does not appear in the output.\n"
        "TeXForm translates standard mathematical functions and operations.\n"
        "Following standard mathematical conventions, single-character symbol names are given in italic font, while multiple character names are given in roman font.");
    symtab_set_docstring("HoldForm", "HoldForm[expr] prints as the expression expr, with expr maintained in an unevaluated form.");
    symtab_set_docstring("ToString",
        "ToString[expr]\n"
        "\tgives a string corresponding to the printed form of expr in InputForm.\n"
        "ToString[expr, form]\n"
        "\tgives the string corresponding to output in the specified form.\n"
        "Supported forms: InputForm (default), FullForm, TeXForm.\n"
        "Examples:\n"
        "\tToString[x^2 + y^3]              -> \"x^2 + y^3\"\n"
        "\tToString[x^2 + y^3, FullForm]    -> \"Plus[Power[x, 2], Power[y, 3]]\"\n"
        "\tToString[x^2 + y^3, TeXForm]     -> \"x^2 + y^3\" (rendered in TeX)");
    symtab_set_docstring("ToExpression",
        "ToExpression[input]\n"
        "\tparses the string input as Mathilda input and returns the\n"
        "\tresulting expression (after evaluation).\n"
        "ToExpression[input, form]\n"
        "\tuses interpretation rules corresponding to the specified form.\n"
        "\tForm may be InputForm or FullForm (both currently use the same parser).\n"
        "ToExpression[input, form, h]\n"
        "\twraps the head h around the parsed expression before evaluation.\n"
        "\tUse h = Hold to obtain the unevaluated parsed form.\n"
        "Returns $Failed if a syntax error is encountered.\n"
        "Examples:\n"
        "\tToExpression[\"1+1\"]                  -> 2\n"
        "\tToExpression[\"1+1\", InputForm, Hold] -> Hold[1 + 1]");
    symtab_set_docstring("Symbol",
        "Symbol[\"name\"]\n"
        "\trefers to a symbol with the specified name.\n"
        "\n"
        "All symbols, whether explicitly entered using Symbol or not, have head Symbol.\n"
        "x_Symbol can be used as a pattern to represent any symbol.\n"
        "The string \"name\" in Symbol[\"name\"] must be an appropriate name for a symbol.\n"
        "It can contain any letters, letter-like forms, or digits, but cannot start with\n"
        "a digit. A backtick (`) separates context prefixes; a leading backtick makes the\n"
        "name relative to the current context $Context.\n"
        "Symbol[\"name\"] creates a new symbol if none exists with the specified name.\n"
        "If Symbol[\"name\"] creates a new symbol, it does so in the context specified by\n"
        "$Context.\n"
        "Examples:\n"
        "\tSymbol[\"x\"]                       -> x\n"
        "\tHead[Symbol[\"x\"]]                 -> Symbol\n"
        "\tSymbol[\"a`x\"]                     -> a`x\n"
        "\t{f[x], f[\"x\"], f[2]} /. f[s_Symbol] :> g[s] -> {g[x], f[\"x\"], f[2]}\n"
        "\n"
        "Attributes: Protected.");
    symtab_set_docstring("ReplaceList", "ReplaceList[expr, rules] attempts to transform the entire expression expr by applying a rule or list of rules in all possible ways, and returns a list of the results obtained.\nReplaceList[expr, rules, n] gives a list of at most n results.");
    symtab_set_docstring("Cases", "Cases[{e1, e2, ...}, pattern] gives a list of the ei that match the pattern.\nCases[{e1, ...}, pattern -> rhs] gives a list of the values of rhs corresponding to the ei that match the pattern.\nCases[expr, pattern, levelspec] gives a list of all parts of expr on levels specified by levelspec that match the pattern.\nCases[expr, pattern -> rhs, levelspec] gives the values of rhs that match the pattern.\nCases[expr, pattern, levelspec, n] gives the first n parts in expr that match the pattern.\nCases[pattern] represents an operator form of Cases that can be applied to an expression.");
    symtab_set_docstring("DeleteCases", "DeleteCases[expr, pattern] removes all elements of expr that match pattern.\nDeleteCases[expr, pattern, levelspec] removes all parts of expr on levels specified by levelspec that match pattern.\nDeleteCases[expr, pattern, levelspec, n] removes the first n parts of expr that match pattern.\nDeleteCases[pattern] represents an operator form of DeleteCases that can be applied to an expression.\nThe default levelspec is {1}. With Heads -> True, the heads of expressions are also tested; deleting a head is equivalent to applying FlattenAt at that location.\nDeleteCases traverses expr in depth-first post-order (leaves before roots).");
    symtab_set_docstring("Position", "Position[expr, pattern] gives a list of the positions at which objects matching pattern appear in expr.\nPosition[expr, pattern, levelspec] finds only objects that appear on levels specified by levelspec.\nPosition[expr, pattern, levelspec, n] gives the positions of the first n objects found.\nPosition[pattern] represents an operator form of Position that can be applied to an expression.");
    symtab_set_docstring("OrderedQ", "OrderedQ[h[e1, e2, ...]] gives True if the elements are in canonical order, and False otherwise.\nOrderedQ[expr, p] uses the ordering function p to determine whether each pair of elements is in order.");
    symtab_set_docstring("Sort", "Sort[list] sorts the elements of list into canonical order.\nSort[list, p] sorts using the ordering function p.");
    symtab_set_docstring("MemberQ", "MemberQ[list, form] returns True if an element of list matches form, and False otherwise.\nMemberQ[list, form, levelspec] tests all parts of list specified by levelspec.\nMemberQ[form] represents an operator form of MemberQ that can be applied to an expression.");
    symtab_set_docstring("Count", "Count[list, pattern] gives the number of elements in list that match pattern.\nCount[expr, pattern, levelspec] gives the total number of subexpressions matching pattern that appear at the levels in expr specified by levelspec.\nCount[pattern] represents an operator form of Count that can be applied to an expression.");

    // Series
    symtab_set_docstring("Series",
        "Series[f, {x, x0, n}]\n"
        "\tgenerates a power-series expansion of f about x = x0 up to order (x - x0)^n.\n"
        "Series[f, x -> x0]\n"
        "\tgenerates the leading term of a power-series expansion of f about x = x0.\n"
        "Series[f, {x, x0, nx}, {y, y0, ny}, ...]\n"
        "\titeratively expands f, first in x, then in y, etc.\n"
        "Series handles Taylor, Laurent (negative powers), and Puiseux (fractional powers)\n"
        "\texpansions, as well as logarithmic and symbolic-exponent cases such as x^x\n"
        "\tand (1+x)^n.\n"
        "Series[f, {x, Infinity, n}] expands around x = Infinity by substituting x -> 1/u.\n"
        "The result of Series is a SeriesData object; use Normal to convert it back to\n"
        "\tan ordinary expression by dropping the O-term.\n"
        "Series is Protected and HoldAll so the expansion variable is not evaluated.");
    symtab_set_docstring("Normal",
        "Normal[expr]\n"
        "\tconverts expr to a normal expression. If expr is a SeriesData object, the\n"
        "\tO-term is dropped and the truncated polynomial (or Laurent/Puiseux sum) is\n"
        "\treturned. Other expressions pass through unchanged.");
    symtab_set_docstring("SeriesData",
        "SeriesData[x, x0, {a0, a1, ...}, nmin, nmax, den]\n"
        "\trepresents a power series in the variable x about the point x0.\n"
        "\tThe ai are the coefficients in the power series. The powers of\n"
        "\t(x - x0) that appear are nmin/den, (nmin+1)/den, ..., (nmax-1)/den,\n"
        "\tand an O[x - x0]^(nmax/den) term represents the omitted higher-order terms.\n"
        "SeriesData objects are generated by Series.\n"
        "SeriesData objects print as sums of coefficients multiplied by powers of x - x0;\n"
        "\tInputForm prints the literal SeriesData[...] form instead.\n"
        "Normal[expr] converts a SeriesData object into a normal expression, dropping the O-term.");

    // Polynomials
    symtab_set_docstring("Expand", "Expand[expr] expands out products and powers in expr.\nExpand[expr, patt] leaves unexpanded any parts of expr that are free of the pattern patt.");
    symtab_set_docstring("ExpandNumerator",
        "ExpandNumerator[expr]\n\texpands out products and powers that appear in the numerator of expr.\n"
        "ExpandNumerator works on terms that have positive integer exponents.\n"
        "ExpandNumerator applies only to the top level in expr.\n"
        "ExpandNumerator does not separate the fraction; Expand does.\n"
        "ExpandNumerator leaves the denominator unexpanded.\n"
        "ExpandNumerator automatically threads over lists, as well as equations,\n"
        "\tinequalities, and logic functions.");
    symtab_set_docstring("ExpandDenominator",
        "ExpandDenominator[expr]\n\texpands out products and powers that appear as denominators in expr.\n"
        "ExpandDenominator works only on negative integer powers.\n"
        "ExpandDenominator applies only to the top level in expr.\n"
        "ExpandDenominator leaves the numerator unexpanded.\n"
        "ExpandDenominator automatically threads over lists, as well as equations,\n"
        "\tinequalities, and logic functions.");
    symtab_set_docstring("Coefficient", "Coefficient[expr, form] gives the coefficient of form in expr.\nCoefficient[expr, form, n] gives the coefficient of form^n in expr.");
    symtab_set_docstring("CoefficientList", "CoefficientList[poly, var] gives a list of coefficients of powers of var in poly, starting with power 0.\nCoefficientList[poly, {var1, var2, ...}] gives an array of coefficients of the variables.");
    symtab_set_docstring("PolynomialGCD",
        "PolynomialGCD[poly1, poly2, ...] gives the greatest common divisor of the polynomials.\n"
        "Option Extension -> alpha computes the GCD over Q(alpha), where alpha is an\n"
        "algebraic number recognised by qa_resolve_extension (Sqrt[c], c^(1/n), or I).\n"
        "Default Extension -> None and Extension -> Automatic compute over the rationals,\n"
        "treating any algebraic numbers in the input as independent variables.");
    symtab_set_docstring("PolynomialLCM",
        "PolynomialLCM[poly1, poly2, ...] gives the least common multiple of the polynomials.\n"
        "Option Extension -> alpha computes the LCM over Q(alpha) via\n"
        "lcm(a, b) = a*b / PolynomialGCD[a, b, Extension -> alpha].\n"
        "Default Extension -> None computes over the rationals.");
    symtab_set_docstring("PolynomialExtendedGCD", "PolynomialExtendedGCD[poly1, poly2, x] gives the extended GCD of poly1 and poly2 treated as univariate polynomials in x.");
    symtab_set_docstring("PolynomialQuotient",
        "PolynomialQuotient[p, q, x] gives the quotient of p and q, treated as polynomials in x, with any remainder dropped.\n"
        "Option: Extension -> alpha (default None) divides over Q(alpha) using the Q(alpha)[x] long-division substrate; Sqrt[c], c^(1/n), and I are recognised forms for alpha.\n"
        "Extension -> None and Extension -> Automatic are accepted and currently behave as the default (no extension).");
    symtab_set_docstring("PolynomialRemainder",
        "PolynomialRemainder[p, q, x] gives the remainder from dividing p by q, treated as polynomials in x.\n"
        "Option: Extension -> alpha (default None) computes the remainder over Q(alpha); see PolynomialQuotient for the recognised alpha forms.");
    symtab_set_docstring("Collect", "Collect[expr, x] collects together terms involving the same powers of objects matching x.");
    symtab_set_docstring("PolynomialMod", "PolynomialMod[poly,m] gives the polynomial poly reduced modulo m.\nPolynomialMod[poly,{Subscript[m, 1],Subscript[m, 2],...}] reduces modulo all of the Subscript[m, i].");
    symtab_set_docstring("FactorSquareFree", "FactorSquareFree[poly] pulls out any multiple factors in a polynomial.");
    symtab_set_docstring("Factor",
        "Factor[poly] factors a polynomial over the integers.\n"
        "Factor[poly, Extension -> alpha] factors over Q(alpha), where alpha is\n"
        "Sqrt[c], c^(1/n) (rational c), or I.  Implements Trager's algebraic-\n"
        "factoring algorithm via norm + sqfr_norm + alg_factor.\n"
        "Factor[poly, Extension -> {alpha_1, ..., alpha_n}] factors over the\n"
        "compositum Q(alpha_1, ..., alpha_n).  The tower is reduced to a single\n"
        "primitive element gamma = alpha_1 + s_2 alpha_2 + ... via Trager's\n"
        "primitive-element algorithm (Phase G6).");
    symtab_set_docstring("FactorTerms",
        "FactorTerms[poly]\n\tpulls out any overall numerical factor in poly.\n"
        "FactorTerms[poly, x]\n\tpulls out any overall factor in poly that does not depend on x.\n"
        "FactorTerms[poly, {x1, x2, ...}]\n\tpulls out any overall factor in poly that does not depend on any of the xi, then progressively factors with respect to smaller subsets {x1, ..., x_{k-1}}.\n"
        "FactorTerms[poly, x] extracts the content of poly with respect to x.\n"
        "FactorTerms automatically threads over lists, equations, inequalities and logic functions.");
    symtab_set_docstring("FactorTermsList",
        "FactorTermsList[poly]\n\tgives a list in which the first element is the overall numerical factor in poly, and the second element is the polynomial with the overall factor removed.\n"
        "FactorTermsList[poly, {x1, x2, ...}]\n\tgives a list of factors of poly. The first element in the list is the overall numerical factor. The second element is a factor that does not depend on any of the xi. Subsequent elements are factors which depend on progressively more of the xi.");
    symtab_set_docstring("Variables", "Variables[poly] gives a list of all independent variables in the polynomial poly.");
    symtab_set_docstring("PolynomialQ", "PolynomialQ[expr, var] yields True if expr is a polynomial in var.");
    symtab_set_docstring("Decompose", "Decompose[poly, x] decomposes a polynomial, if possible, into a composition of simpler polynomials.");
    symtab_set_docstring("HornerForm", "HornerForm[poly] puts the polynomial poly in Horner form.\nHornerForm[poly, vars] puts poly in Horner form with respect to vars.\nHornerForm[poly1/poly2, vars1, vars2] puts the rational function in Horner form.");
    symtab_set_docstring("Resultant", "Resultant[poly1, poly2, var] computes the resultant of the polynomials poly1 and poly2 with respect to the variable var.");
    symtab_set_docstring("Discriminant", "Discriminant[poly, var] computes the discriminant of the polynomial poly with respect to the variable var.");
    symtab_set_docstring("Numerator", "Numerator[expr] gives the numerator of expr.\nNumerator picks out terms which do not have superficially negative exponents.");
    symtab_set_docstring("Denominator", "Denominator[expr] gives the denominator of expr.\nDenominator picks out terms which have superficially negative exponents.");
    symtab_set_docstring("Cancel",
        "Cancel[expr] cancels out common factors in the numerator and denominator of expr.\n"
        "Option Extension -> alpha cancels factors over Q(alpha) (e.g. simplifies\n"
        "(x^2 - 2)/(x - Sqrt[2]) to x + Sqrt[2] when Extension -> Sqrt[2]).\n"
        "Default Extension -> None treats algebraic numbers as opaque.");
    symtab_set_docstring("Together",
        "Together[expr] combines fractions over a common denominator, then cancels.\n"
        "Option Extension -> alpha runs the final cancellation over Q(alpha) so\n"
        "algebraic-coefficient simplifications fire (e.g. 1/(x-Sqrt[2]) + 1/(x+Sqrt[2])\n"
        "collapses to (2 x)/(x^2 - 2) under Extension -> Sqrt[2]).\n"
        "Default Extension -> None combines and cancels over the rationals only.");
    symtab_set_docstring("Rationalize",
        "Rationalize[x]\n"
        "\tconverts an approximate number x to a nearby rational with small denominator.\n"
        "Rationalize[x, dx]\n"
        "\tyields the rational number with smallest denominator that lies within dx of x.\n"
        "\n"
        "Rationalize[x] yields x unchanged if there is no rational number close enough to x to satisfy |p/q - x| < c/q^2, with c = 10^-4.\n"
        "Rationalize[x, dx] works with exact numbers x: the value is first numericalised, then rationalised.\n"
        "Rationalize[x, 0] forces conversion of any inexact number x to rational form, using a tolerance derived from the precision of x.\n"
        "Rationalize threads over compound expressions and Complex[re, im], so e.g. Rationalize[1.2 + 6.7 x] gives 6/5 + (67 x)/10.");
    symtab_set_docstring("Apart",
        "Apart[expr] rewrites a rational expression as a sum of terms with minimal denominators.\n"
        "Apart[expr, var] treats all variables other than var as constants.\n"
        "Option Extension -> alpha factors the denominator over Q(alpha) before\n"
        "decomposition, splitting (x^2 - 2) into (x - Sqrt[2])(x + Sqrt[2]) under\n"
        "Extension -> Sqrt[2] and producing the corresponding linear-factor\n"
        "partial fractions.  Default Extension -> None decomposes over Q.");
    symtab_set_docstring("Simplify",
        "Simplify[expr]\n\tperforms a sequence of algebraic and other transformations on expr and returns the simplest form it finds.\n"
        "Simplify[expr, assum]\n\tdoes simplification using assumptions assum.\n"
        "Simplify accepts the options ComplexityFunction (default: leaf count plus integer-digit count, matching Mathematica) and Assumptions (default: $Assumptions).\n"
        "Simplify tries Together, Cancel, Expand, Factor, FactorSquareFree, Apart, TrigExpand, TrigFactor, and a TrigToExp/ExpToTrig roundtrip, keeping the smallest result.\n"
        "Under positivity / reality assumptions Simplify also applies Log/Power identities -- Log[a b] -> Log[a] + Log[b], (a b)^c -> a^c b^c, (a^p)^q -> a^(p q), Log[a^p] -> p Log[a] and the like -- whenever the operand-domain conditions are provable from the assumption set.\n"
        "Assumptions can be equations, inequalities, domain specifications such as Element[x, Integers], or logical combinations of these. Lists of assumptions are converted to conjunctions.\n"
        "Simplify automatically threads over lists, equations, inequalities, and logic functions.");
    symtab_set_docstring("Assuming",
        "Assuming[assum, expr]\n\tevaluates expr with assum appended to $Assumptions, so that assum is included in the default assumptions used by functions such as Simplify.\n"
        "Assuming converts lists of assumptions to conjunctions.\n"
        "Assuming[assum, expr] is effectively equivalent to Block[{$Assumptions = $Assumptions && assum}, expr], so nested invocations compose and the rebinding of $Assumptions is restored on exit.");
    symtab_set_docstring("$Assumptions",
        "$Assumptions\n\tis the default setting for the Assumptions option used in Simplify and other functions that take assumptions.\n"
        "$Assumptions defaults to True (no assumptions). Functions like Assuming temporarily extend $Assumptions for the duration of their body.");
    symtab_set_docstring("Element",
        "Element[x, dom]\n\treturns True if x is provably an element of the domain dom under the current $Assumptions, False if it is provably not, and stays unevaluated otherwise.\n"
        "Supported domains: Integers, Rationals, Reals, Algebraics, Complexes, Booleans, Primes, Composites.\n"
        "Numeric and structural literals decide directly: Element[5, Integers] -> True, Element[5/2, Integers] -> False, Element[1+I, Reals] -> False, Element[2.5, Integers] -> False.\n"
        "Element consults $Assumptions for symbolic queries, so under Assuming[Element[x, Integers], ...] a query Element[x, Reals] returns True via the Integer => Real lattice.\n"
        "Element threads over lists in its first argument: Element[{x1, x2, ...}, dom] returns the list of per-element decisions.");
    symtab_set_docstring("LeafCount", "LeafCount[expr] gives the total number of indivisible subexpressions in expr.");
    symtab_set_docstring("ByteCount", "ByteCount[expr] gives the number of bytes used internally by Mathilda to store expr.");
    symtab_set_docstring("Level",
        "Level[expr, levelspec]\n"
        "\tgives a list of all subexpressions of expr on levels specified by levelspec.\n"
        "Level[expr, levelspec, f]\n"
        "\tapplies f to the sequence of subexpressions.\n\n"
        "Level uses standard level specifications:\n"
        "  n            levels 1 through n\n"
        "  Infinity     levels 1 through Infinity\n"
        "  {n}          level n only\n"
        "  {n1, n2}     levels n1 through n2\n\n"
        "Level[expr, {-1}] gives a list of all \"atomic\" objects in expr.\n"
        "A positive level n consists of all parts of expr specified by n indices.\n"
        "A negative level -n consists of all parts of expr with depth n.\n"
        "Level 0 corresponds to the whole expression.\n"
        "With the option setting Heads->True, Level includes heads of expressions and their parts.\n"
        "Level traverses expressions in depth-first order, so that the subexpressions in the final list are ordered lexicographically by their indices.");

    // Time and Date
    symtab_set_docstring("Timing", "Timing[expr] evaluates expr, and returns a list of the time in seconds used, together with the result obtained.");
    symtab_set_docstring("RepeatedTiming", "RepeatedTiming[expr] evaluates expr repeatedly and returns a list of the average time in seconds used, together with the result obtained.\nRepeatedTiming[expr, t] does repeated evaluation for at least t seconds.");
    symtab_set_docstring("TimeConstrained",
        "TimeConstrained[expr, t]\n"
        "\tevaluates expr, stopping after t seconds.\n"
        "TimeConstrained[expr, t, failexpr]\n"
        "\treturns failexpr if the time constraint is not met.\n"
        "\n"
        "TimeConstrained generates an interrupt to abort the evaluation of\n"
        "expr if the evaluation is not completed within the specified time.\n"
        "TimeConstrained evaluates failexpr only if the evaluation is\n"
        "aborted. TimeConstrained returns $Aborted if the evaluation is\n"
        "aborted and no failexpr is specified.\n"
        "\n"
        "TimeConstrained[expr, Infinity] imposes no time constraint.\n"
        "TimeConstrained may give different results on different occasions\n"
        "within a single session, for example as a result of different\n"
        "conditions of internal system caches.\n"
        "TimeConstrained takes account only of CPU time spent inside the\n"
        "main Mathilda kernel process; it does not include additional\n"
        "threads or processes.");
    symtab_set_docstring("AbsoluteTime",
        "AbsoluteTime[]\n"
        "\tgives the total number of seconds since the beginning of January 1, 1900, in your time zone.\n"
        "AbsoluteTime[date]\n"
        "\tgives the absolute time specification corresponding to the given date specification.\n"
        "\n"
        "The supported date specifications are:\n"
        "\t{y, m, d, h, m, s}\tDateList specification\n"
        "\ttime\t\t\tAbsoluteTime specification (a number, returned unchanged)\n"
        "\n"
        "DateList entries may be elided from the right: {y}, {y, m}, {y, m, d}, etc. fill the\n"
        "missing fields with {_, 1, 1, 0, 0, 0}. Day, hour, minute, and second values may be\n"
        "noninteger; the year and month must be integers. Date lists are converted to standard\n"
        "normalized form, so e.g. AbsoluteTime[{2022, 2, 31}] = AbsoluteTime[{2022, 3, 3}].\n"
        "\n"
        "AbsoluteTime[] uses whatever date and time have been set on your computer system. It\n"
        "performs no corrections for time zones, daylight saving time, or leap seconds.");

    // Comparisons
    symtab_set_docstring("SameQ", "lhs === rhs yields True if the expression lhs is identical to rhs, and yields False otherwise.");
    symtab_set_docstring("UnsameQ", "lhs =!= rhs yields True if the expression lhs is not identical to rhs, and yields False otherwise.");
    symtab_set_docstring("Equal", "lhs == rhs yields True if lhs and rhs are equal.");
    symtab_set_docstring("Unequal", "lhs != rhs yields True if lhs and rhs are not equal.");
    symtab_set_docstring("Less", "x < y yields True if x is strictly less than y.");
    symtab_set_docstring("Greater", "x > y yields True if x is strictly greater than y.");
    symtab_set_docstring("LessEqual", "x <= y yields True if x is less than or equal to y.");
    symtab_set_docstring("GreaterEqual", "x >= y yields True if x is greater than or equal to y.");

    // Primes
    symtab_set_docstring("FactorInteger", "FactorInteger[n] gives a list of the prime factors of the integer n, together with their exponents.");
    symtab_set_docstring("EulerPhi", "EulerPhi[n] gives the Euler totient function phi(n).");
    symtab_set_docstring("PrimePi", "PrimePi[x] gives the number of primes less than or equal to x.");
    symtab_set_docstring("NextPrime", "NextPrime[x] gives the next prime after x.");
    symtab_set_docstring("Distribute", "Distribute[f[x1, x2, ...]]\n\tdistributes f over Plus appearing in any of the xi.\nDistribute[expr, g]\n\tdistributes over g.\nDistribute[expr, g, f]\n\tperforms the distribution only if the head of expr is f.\nDistribute[expr, g, f, gp, fp]\n\tgives gp and fp in place of g and f respectively in the result of the distribution.");
    symtab_set_docstring("Distribute", "Distribute[f[x1, x2, ...]]\n\tdistributes f over Plus appearing in any of the xi.\nDistribute[expr, g]\n\tdistributes over g.\nDistribute[expr, g, f]\n\tperforms the distribution only if the head of expr is f.\nDistribute[expr, g, f, gp, fp]\n\tgives gp and fp in place of g and f respectively in the result of the distribution.");
    symtab_set_docstring("Inner", "Inner[f,list1,list2,g]\n\tis a generalization of Dot in which f plays the role of multiplication and g of addition.\nInner[f,list1,list2]\n\tuses Plus for g.\nInner[f,list1,list2,g,n]\n\tcontracts index n of the first tensor with the first index of the second tensor.");
    symtab_set_docstring("Outer", "Outer[f,list1,list2,...]\n\tgives the generalized outer product of the listi, forming all possible combinations of the lowest-level elements in each of them, and feeding them as arguments to f.\nOuter[f,list1,list2,...,n]\n\ttreats as separate elements only sublists at level n in the listi.\nOuter[f,list1,list2,...,n1,n2,...]\n\ttreats as separate elements only sublists at level ni in the corresponding listi.");
    symtab_set_docstring("Tuples", "Tuples[list,n]\n\tgenerates a list of all possible n-tuples of elements from list.\nTuples[{list1,list2,...}]\n\tgenerates a list of all possible tuples whose ith element is from listi.\nTuples[list,{n1,n2,...}]\n\tgenerates a list of all possible n1 x n2 x ... arrays of elements in list.");
    symtab_set_docstring("Permutations", "Permutations[list]\n\tgenerates a list of all possible permutations of the elements in list.\nPermutations[list,n]\n\tgives all permutations containing at most n elements.\nPermutations[list,{n}]\n\tgives all permutations containing exactly n elements.");
    symtab_set_docstring("Min", "Min[x1, x2, ...]\n\tyields the numerically smallest of the xi.\nMin[{x1, x2, ...}, {y1, ...}, ...]\n\tyields the smallest element of any of the lists.");
    symtab_set_docstring("Max", "Max[x1, x2, ...]\n\tyields the numerically largest of the xi.\nMax[{x1, x2, ...}, {y1, ...}, ...]\n\tyields the largest element of any of the lists.");
    symtab_set_docstring("Median", "Median[data]\n\tgives the median estimate of the elements in data.\nMedian[dist]\n\tgives the median of the distribution dist.");

    symtab_set_docstring("Quartiles", "Quartiles[data]\n\tgives the {q_1/4, q_2/4, q_3/4} quantile estimates of the elements in data.\nQuartiles[data,{{a,b},{c,d}}]\n\tuses the quantile definition specified by parameters a, b, c, d.\nQuartiles[dist]\n\tgives the {q_1/4, q_2/4, q_3/4} quantiles of the distribution dist.");

    // Random
    symtab_set_docstring("RandomInteger",
        "RandomInteger[{imin, imax}]\n\tgives a pseudorandom integer in the range {imin, ..., imax}.\n"
        "RandomInteger[imax]\n\tgives a pseudorandom integer in the range {0, ..., imax}.\n"
        "RandomInteger[]\n\tpseudorandomly gives 0 or 1.\n"
        "RandomInteger[range, n]\n\tgives a list of n pseudorandom integers.\n"
        "RandomInteger[range, {n1, n2, ...}]\n\tgives an n1 x n2 x ... array of pseudorandom integers.");
    symtab_set_docstring("RandomReal",
        "RandomReal[]\n\tgives a pseudorandom real number in the range 0 to 1.\n"
        "RandomReal[{xmin, xmax}]\n\tgives a pseudorandom real number in the range xmin to xmax.\n"
        "RandomReal[xmax]\n\tgives a pseudorandom real number in the range 0 to xmax.\n"
        "RandomReal[range, n]\n\tgives a list of n pseudorandom reals.\n"
        "RandomReal[range, {n1, n2, ...}]\n\tgives an n1 x n2 x ... array of pseudorandom reals.\n"
        "RandomReal[spec, WorkingPrecision -> n]\n\tyields reals with n digits of precision.\n"
        "\tLeading or trailing digits of the generated number can be 0.\n"
        "\tn may be MachinePrecision (the default) or a positive number of decimal digits.");
    symtab_set_docstring("SeedRandom",
        "SeedRandom[n]\n\tseeds the pseudorandom generator with the integer n.\n"
        "SeedRandom[]\n\treseeds the pseudorandom generator from system entropy.");
    symtab_set_docstring("RandomComplex",
        "RandomComplex[]\n\tgives a pseudorandom complex number with real and imaginary parts in the range 0 to 1.\n"
        "RandomComplex[{zmin, zmax}]\n\tgives a pseudorandom complex number in the rectangle with corners given by the complex numbers zmin and zmax.\n"
        "RandomComplex[zmax]\n\tgives a pseudorandom complex number in the rectangle whose corners are the origin and zmax.\n"
        "RandomComplex[range, n]\n\tgives a list of n pseudorandom complex numbers.\n"
        "RandomComplex[range, {n1, n2, ...}]\n\tgives an n1 x n2 x ... array of pseudorandom complex numbers.\n"
        "RandomComplex[spec, WorkingPrecision -> n]\n\tyields complex numbers whose real and imaginary parts have n digits of precision.\n"
        "\tLeading or trailing digits of the generated parts can be 0.\n"
        "\tn may be MachinePrecision (the default) or a positive number of decimal digits.");
    symtab_set_docstring("RandomChoice",
        "RandomChoice[{e1, e2, ...}]\n\tgives a pseudorandom choice of one of the ei.\n"
        "RandomChoice[list, n]\n\tgives a list of n pseudorandom choices.\n"
        "RandomChoice[list, {n1, n2, ...}]\n\tgives an n1 x n2 x ... array of pseudorandom choices.\n"
        "RandomChoice[{w1, w2, ...} -> {e1, e2, ...}]\n\tgives a pseudorandom choice weighted by the wi.\n"
        "RandomChoice[wlist -> elist, n]\n\tgives a list of n weighted choices.\n"
        "RandomChoice[wlist -> elist, {n1, n2, ...}]\n\tgives an n1 x n2 x ... array of weighted choices.");
    symtab_set_docstring("RandomSample",
        "RandomSample[{e1, e2, ...}, n]\n\tgives a pseudorandom sample of n of the ei, without replacement.\n"
        "RandomSample[{w1, w2, ...} -> {e1, e2, ...}, n]\n\tgives a weighted pseudorandom sample of n of the ei.\n"
        "RandomSample[{e1, e2, ...}]\n\tgives a pseudorandom permutation of the ei.\n"
        "RandomSample[list, UpTo[n]]\n\tgives a sample of n of the ei, or as many as are available.\n"
        "RandomSample never samples any element more than once.\n"
        "Use SeedRandom to seed the pseudorandom generator for reproducible results.");
}
