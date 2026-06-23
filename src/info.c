#include "info.h"
#include "symtab.h"

void info_init(void) {
    // Arithmetic
    symtab_set_docstring("Plus",
        "x + y + ... or Plus[x, y, ...] represents a sum of terms.\n"
        "Plus is Flat, Orderless, OneIdentity, Listable, and NumericFunction:\n"
        "nested Plus is auto-flattened, terms are sorted into canonical order,\n"
        "like terms are combined, and integer arguments are summed exactly via\n"
        "GMP (with int64 fast path and BigInt overflow promotion).");
    symtab_set_docstring("Times",
        "x * y * ... or Times[x, y, ...] represents a product of terms.\n"
        "Times is Flat, Orderless, OneIdentity, Listable, and NumericFunction:\n"
        "nested Times is auto-flattened, factors are sorted, like factors are\n"
        "merged into Power, and integer products use exact GMP arithmetic.\n"
        "Numeric zero collapses the product; a Plus factor is left distributed\n"
        "(use Expand to distribute).");
    symtab_set_docstring("Power",
        "x ^ y or Power[x, y] represents x to the power y.\n"
        "Power is Listable, NumericFunction, and OneIdentity. Integer exponents\n"
        "are reduced exactly (repeated squaring on GMP); Rational and Real\n"
        "exponents evaluate numerically when the base is numeric; Power[0, 0]\n"
        "stays Indeterminate; Power[x, 1/2] is canonicalised to Sqrt[x].");
    symtab_set_docstring("Subtract",
        "x - y or Subtract[x, y] represents x - y; rewritten by the evaluator\n"
        "to Plus[x, Times[-1, y]] so it inherits Plus's flattening and ordering.");
    symtab_set_docstring("Divide",
        "x / y or Divide[x, y] represents x / y; rewritten by the evaluator to\n"
        "Times[x, Power[y, -1]] so it inherits Times's flattening and ordering.");
    symtab_set_docstring("Sqrt",
        "Sqrt[z]\n"
        "\trepresents the principal square root of z.\n"
        "Sqrt is Listable. Sqrt[z] is canonicalised to Power[z, 1/2]; perfect\n"
        "integer / rational squares reduce to exact form, negative real inputs\n"
        "yield I * Sqrt[-x], and numeric inputs (Real / MPFR / Complex) are\n"
        "evaluated directly. Branch cut along the negative real axis.");
    symtab_set_docstring("Complex",
        "Complex[re, im]\n"
        "\trepresents the complex number re + im I.\n"
        "Complex is the canonical head produced by arithmetic when an Integer,\n"
        "Real, or Rational acquires an imaginary part. Pure-real inputs collapse\n"
        "to the underlying number; im == 0 unwraps to re.");
    symtab_set_docstring("Rational",
        "Rational[n, d]\n"
        "\trepresents the rational number n/d.\n"
        "When n and d are integers, Rational auto-reduces by gcd, normalises\n"
        "the sign onto the numerator, and collapses to an Integer when d == 1.\n"
        "Rationals propagate through Plus / Times exactly via GMP.");
    symtab_set_docstring("GCD",
        "GCD[n1, n2, ...]\n"
        "\tgives the greatest common divisor of the integers ni.\n"
        "Computed via GMP's binary-GCD (mpz_gcd) folded across the arguments.\n"
        "Accepts BigInt and Rational inputs (gcd(p1/q1, p2/q2) = gcd(p1,p2) /\n"
        "lcm(q1,q2)); non-integer Real or symbolic inputs leave GCD unevaluated.");
    symtab_set_docstring("LCM",
        "LCM[n1, n2, ...]\n"
        "\tgives the least common multiple of the integers ni.\n"
        "Computed via GMP's mpz_lcm folded across the arguments; sign is\n"
        "normalised non-negative. Accepts BigInt and Rational inputs.");
    symtab_set_docstring("ExtendedGCD",
        "ExtendedGCD[n1, n2, ...]\n"
        "\tgives the extended GCD {g, {r1, r2, ...}} of the integers ni,\n"
        "\twhere g == GCD[n1, ...] and g == r1 n1 + r2 n2 + ....\n"
        "Computed by folding GMP's mpz_gcdext pairwise; accepts machine and\n"
        "BigInt integers and threads over lists. Non-integer or inexact\n"
        "arguments leave ExtendedGCD unevaluated.");
    symtab_set_docstring("Divisible",
        "Divisible[n, m]\n"
        "\tyields True if n is divisible by m, and False otherwise.\n"
        "n is divisible by m when n is an integer multiple of m; this is\n"
        "effectively Mod[n, m] == 0.  Works for machine and BigInt integers,\n"
        "Gaussian integers, rationals, and exact numeric quantities (the\n"
        "quotient n/m must reduce to an integer or Gaussian integer).  Returns\n"
        "False unless n and m are manifestly divisible; symbolic, non-numeric\n"
        "arguments are left unevaluated.  Listable.");
    symtab_set_docstring("CoprimeQ",
        "CoprimeQ[n1, n2, ...]\n"
        "\tyields True if the arguments are pairwise relatively prime, and\n"
        "\tFalse otherwise.\n"
        "Integers are relatively prime when their GCD is 1.  Works for machine\n"
        "and BigInt integers.  With GaussianIntegers -> True, or when any\n"
        "argument is an exact Gaussian integer, coprimality is tested over the\n"
        "Gaussian integers Z[i].  Returns False unless the arguments are\n"
        "manifestly coprime; CoprimeQ[] is False and CoprimeQ[n] is True.\n"
        "Listable and Orderless.");
    symtab_set_docstring("PowerMod", "PowerMod[a, b, m] gives a^b mod m.\nPowerMod[a, -1, m] finds the modular inverse of a modulo m.\nPowerMod[a, 1/r, m] finds a modular r-th root of a.");
    symtab_set_docstring("PrimitiveRoot",
        "PrimitiveRoot[n]\n"
        "\tgives a primitive root of n.\n"
        "PrimitiveRoot[n, k]\n"
        "\tgives the smallest primitive root of n greater than or equal to k.\n\n"
        "A primitive root of n is a generator of the multiplicative group of integers modulo n "
        "relatively prime to n.  PrimitiveRoot returns unevaluated unless n is 2, 4, an odd prime "
        "power p^k, or twice an odd prime power 2 p^k.");
    symtab_set_docstring("PrimitiveRootList",
        "PrimitiveRootList[n]\n"
        "\tgives the sorted list of all primitive roots of n in the canonical residues {1, ..., n-1}.\n\n"
        "Returns an empty list unless n is 2, 4, an odd prime power p^k, or twice an odd prime "
        "power 2 p^k.");
    symtab_set_docstring("MultiplicativeOrder",
        "MultiplicativeOrder[k, n]\n"
        "\tgives the multiplicative order of k modulo n, the smallest positive integer m "
        "such that k^m is congruent to 1 modulo n.\n"
        "MultiplicativeOrder[k, n, {r1, r2, ...}]\n"
        "\tgives the smallest positive integer m such that k^m is congruent to one of the ri modulo n.\n\n"
        "Returns unevaluated when gcd(k, n) is not 1, when no power of k lands in the residue set, "
        "or when n is zero.  All arithmetic is exact via GMP, so k and n may be arbitrary-precision "
        "integers.");
    symtab_set_docstring("FindMinimum",
        "FindMinimum[f, {x, x0}]\n"
        "\tsearches for a local minimum of f starting from x = x0.\n"
        "FindMinimum[f, {x, x0, x1}]\n"
        "\tderivative-free 1D search bracketing the minimum from two starts (Brent).\n"
        "FindMinimum[f, {x, xstart, xmin, xmax}]\n"
        "\tbracketed 1D Brent search on [xmin, xmax] starting from xstart.\n"
        "FindMinimum[f, {{x, x0}, {y, y0}, ...}]\n"
        "\tn-D local minimum from a user-supplied start.\n"
        "FindMinimum[f, {x, y, ...}]\n"
        "\tn-D local minimum auto-starting each variable at 0.\n"
        "FindMinimum[{f, cons}, vars]\n"
        "\tlocal minimum subject to box and Inequality constraints.\n\n"
        "Methods (Method -> ...):\n"
        "\tAutomatic           picks Brent for 1D, QuasiNewton (BFGS) for n-D.\n"
        "\t\"Brent\"             derivative-free golden-section + parabolic "
        "interpolation; 1D only; honours MPFR WorkingPrecision.\n"
        "\t\"QuasiNewton\"       BFGS with cubic line search; uses analytic "
        "Gradient if given, otherwise central differences; default for n>=2; "
        "honours MPFR WorkingPrecision.\n"
        "\t\"ConjugateGradient\" Polak-Ribiere CG with line search; lower "
        "memory than BFGS for large n; gradient-based.\n"
        "\t\"Newton\"            full Hessian step via modified Cholesky "
        "factorization; falls back to a steepest-descent step when the "
        "Hessian is not positive definite or unavailable.\n\n"
        "Options:\n"
        "\tMethod              algorithm selector (see above).\n"
        "\tWorkingPrecision    MachinePrecision (double) or a positive digit "
        "count (MPFR; honoured by Brent and BFGS).\n"
        "\tMaxIterations       positive integer cap on outer iterations; "
        "default 500.\n"
        "\tAccuracyGoal        Automatic | Infinity | digits; absolute "
        "tolerance on |f| (and |x| where applicable).\n"
        "\tPrecisionGoal       Automatic | Infinity | digits; relative "
        "tolerance on step size.\n"
        "\tGradient            Automatic (finite differences) or an explicit "
        "list { dfdx1, dfdx2, ... } in the same order as vars.\n"
        "\tStepMonitor         :> body run after each accepted step, with "
        "the variables locally bound to their current values.\n"
        "\tEvaluationMonitor   :> body run on every function/gradient "
        "evaluation.\n\n"
        "FindMinimum has HoldAll and effectively uses Block to localize the "
        "variables.  Returns {fmin, {x -> xmin, ...}}.");
    symtab_set_docstring("FindMaximum",
        "FindMaximum[f, {x, x0}]\n"
        "\tsearches for a local maximum of f starting from x = x0.\n"
        "FindMaximum[f, {x, x0, x1}]\n"
        "\tderivative-free 1D search bracketing the maximum from two starts (Brent on -f).\n"
        "FindMaximum[f, {x, xstart, xmin, xmax}]\n"
        "\tbracketed 1D Brent search on [xmin, xmax] starting from xstart.\n"
        "FindMaximum[f, {{x, x0}, {y, y0}, ...}]\n"
        "\tn-D local maximum from a user-supplied start.\n"
        "FindMaximum[f, {x, y, ...}]\n"
        "\tn-D local maximum auto-starting each variable at 0.\n"
        "FindMaximum[{f, cons}, vars]\n"
        "\tlocal maximum subject to box and Inequality constraints.\n\n"
        "Methods (Method -> ...):\n"
        "\tAutomatic           picks Brent for 1D, QuasiNewton (BFGS) for n-D.\n"
        "\t\"Brent\"             derivative-free golden-section + parabolic "
        "interpolation; 1D only; honours MPFR WorkingPrecision.\n"
        "\t\"QuasiNewton\"       BFGS with cubic line search; uses analytic "
        "Gradient if given, otherwise central differences; default for n>=2; "
        "honours MPFR WorkingPrecision.\n"
        "\t\"ConjugateGradient\" Polak-Ribiere CG with line search; lower "
        "memory than BFGS for large n; gradient-based.\n"
        "\t\"Newton\"            full Hessian step via modified Cholesky "
        "factorization; falls back to a steepest-descent step when the "
        "Hessian is not negative definite or unavailable.\n\n"
        "Options:\n"
        "\tMethod              algorithm selector (see above).\n"
        "\tWorkingPrecision    MachinePrecision (double) or a positive digit "
        "count (MPFR; honoured by Brent and BFGS).\n"
        "\tMaxIterations       positive integer cap on outer iterations; "
        "default 500.\n"
        "\tAccuracyGoal        Automatic | Infinity | digits; absolute "
        "tolerance on |f| (and |x| where applicable).\n"
        "\tPrecisionGoal       Automatic | Infinity | digits; relative "
        "tolerance on step size.\n"
        "\tGradient            Automatic (finite differences) or an explicit "
        "list { dfdx1, dfdx2, ... } in the same order as vars.  The gradient "
        "is taken with respect to f, not -f.\n"
        "\tStepMonitor         :> body run after each accepted step, with "
        "the variables locally bound to their current values.\n"
        "\tEvaluationMonitor   :> body run on every function/gradient "
        "evaluation.\n\n"
        "FindMaximum has HoldAll and effectively uses Block to localize the "
        "variables.  Internally maximises by minimising -f, then negates the "
        "objective value in the result.  Returns {fmax, {x -> xmax, ...}}.");
    symtab_set_docstring("FindRoot",
        "FindRoot[f, {x, x0}]\n"
        "\tsearches for a numerical root of f starting from x = x0.\n"
        "FindRoot[lhs == rhs, {x, x0}]\n"
        "\tsearches for a numerical solution to the equation.\n"
        "FindRoot[f, {x, x0, x1}]\n"
        "\tuses a variant of the secant method with x0 and x1 as the first two approximations.\n"
        "FindRoot[f, {x, xstart, xmin, xmax}]\n"
        "\tuses Brent's method on the bracket [xmin, xmax].\n"
        "FindRoot[{f1, f2, ...}, {{x, x0}, {y, y0}, ...}]\n"
        "\tsearches for a simultaneous numerical root of the system.\n\n"
        "Options: Method ('Newton' | 'Secant' | 'Brent' | Automatic), WorkingPrecision, "
        "MaxIterations, AccuracyGoal, PrecisionGoal, DampingFactor, Jacobian, StepMonitor, "
        "EvaluationMonitor.  FindRoot has HoldAll and effectively uses Block to localize variables.");

    symtab_set_docstring("NResidue",
        "NResidue[expr, {z, z0}]\n"
        "\tnumerically finds the residue of expr near z = z0 (the coefficient of "
        "(z - z0)^-1 in the Laurent expansion) by integrating around a small circle "
        "in the complex plane.\n"
        "NResidue[{e1, e2, ...}, {z, z0}]\n"
        "\tthreads element-wise over the first argument.\n\n"
        "Works for essential singularities where the symbolic Residue (which needs a "
        "power series) cannot. Cannot distinguish a tiny spurious residual from a true "
        "zero -- Chop the result when needed; returns an incorrect value if the contour "
        "encloses another singularity or crosses a branch cut.\n\n"
        "Options: Radius (contour radius, default 1/100, or Automatic), WorkingPrecision, "
        "PrecisionGoal, MaxRecursion (max contour refinements, default 10), "
        "Method ('Trapezoidal').");
    symtab_set_docstring("ND",
        "ND[expr, x, x0]\n"
        "\tgives a numerical approximation to the derivative of expr with respect "
        "to x at the point x0.\n"
        "ND[expr, {x, n}, x0]\n"
        "\tgives a numerical approximation to the n-th derivative.\n"
        "ND[{e1, e2, ...}, x, x0]\n"
        "\tthreads element-wise over the first argument.\n\n"
        "Default Method -> EulerSum uses Richardson extrapolation of forward, "
        "direction-Scale finite differences (works for non-analytic expr; needs "
        "integer n >= 0). Method -> NIntegrate uses Cauchy's integral formula via "
        "NResidue (needs expr analytic near x0; allows fractional/complex order). "
        "ND cannot recognize small numbers that should be zero -- Chop if needed.\n\n"
        "Options: Method (EulerSum | NIntegrate), Scale (step size / contour radius / "
        "complex direction, default 1), Terms (EulerSum extrapolation terms, default 7), "
        "WorkingPrecision, PrecisionGoal, MaxRecursion.");
    symtab_set_docstring("NSeries",
        "NSeries[f, {x, x0, n}]\n"
        "\tgives a numerical approximation to the series expansion of f about "
        "x = x0, including the terms (x - x0)^-n through (x - x0)^n, as a "
        "SeriesData object.\n\n"
        "f is sampled on a circle in the complex plane centred at x0 and a "
        "discrete Fourier transform of the samples recovers the Taylor or "
        "Laurent coefficients (Cauchy's integral formula). The region of "
        "convergence is the annulus, containing the sampled circle, where f is "
        "analytic. Works for essential singularities (e.g. Sin[x + 1/x]) where "
        "the symbolic Series cannot. Returns an incorrect result if the disk "
        "centred at x0 contains a branch cut of f; for a Laurent series the "
        "SeriesData neglects higher-order poles. No effort is made to justify "
        "the precision of the coefficients, and small spurious residuals are "
        "not recognised as zero -- Chop the result when needed.\n\n"
        "Options: Radius (radius of the sampled circle, default 1), "
        "WorkingPrecision (default MachinePrecision).");
    symtab_set_docstring("NLimit",
        "NLimit[expr, z -> z0]\n"
        "\tnumerically finds the limiting value of expr as z approaches z0.\n\n"
        "A geometric sequence of sample points approaching z0 is constructed "
        "(z0 may be finite, complex, or an infinite point such as Infinity or "
        "I Infinity) and the limit is recovered by sequence acceleration. "
        "Method -> EulerSum (default) uses Richardson/Romberg extrapolation; "
        "Method -> SequenceLimit uses Wynn's epsilon algorithm. expr must be "
        "numerical when z is numerical. Small spurious residuals are not "
        "recognised as zero -- Chop if needed.\n\n"
        "Options: Method (EulerSum | SequenceLimit), WorkingPrecision (default "
        "MachinePrecision), Direction (Automatic == -1, or a complex approach "
        "vector), Scale (initial step / distance, default 1), Terms (default 7), "
        "WynnDegree (SequenceLimit iterations, default 1).");
    symtab_set_docstring("NSum",
        "NSum[f, {i, imin, imax}]\n"
        "\tgives a numerical approximation to the sum of f for i from imin to "
        "imax.\n\n"
        "NSum[f, {i, imin, imax, di}] uses step di. imax may be Infinity. "
        "NSum[f, {i, ...}, {j, ...}, ...] evaluates a multidimensional sum (an "
        "inner bound may depend on an outer index). The index is localised "
        "(HoldAll). Method -> Automatic picks Euler-Maclaurin for monotone "
        "series, the Cohen-Villegas-Zagier method for alternating series, and "
        "Wynn's epsilon (partial-sum acceleration) otherwise; large finite sums "
        "use the difference of two infinite tails. Machine or arbitrary "
        "precision via WorkingPrecision.\n\n"
        "Options: Method (Automatic | EulerMaclaurin | AlternatingSigns | "
        "WynnEpsilon), WorkingPrecision (default MachinePrecision), NSumTerms "
        "(head terms summed explicitly, default 15), NSumExtraTerms, WynnDegree, "
        "VerifyConvergence (default True; a divergent sum gives ComplexInfinity), "
        "AccuracyGoal, PrecisionGoal.");
    symtab_set_docstring("NProduct",
        "NProduct[f, {i, imin, imax}]\n"
        "\tgives a numerical approximation to the product of f for i from imin "
        "to imax.\n\n"
        "NProduct[f, {i, imin, imax, di}] uses step di. imax may be Infinity. "
        "NProduct[f, {i, ...}, {j, ...}, ...] evaluates a multidimensional "
        "product (an inner bound may depend on an outer index). The index is "
        "localised (HoldAll). Evaluated as Exp[NSum[Log[f], ...]], so the NSum "
        "engine (Euler-Maclaurin for monotone factors, Wynn's epsilon "
        "otherwise) and its convergence test carry over. Machine or arbitrary "
        "precision via WorkingPrecision.\n\n"
        "Options: Method (Automatic | EulerMaclaurin | WynnEpsilon), "
        "WorkingPrecision (default MachinePrecision), NProductFactors (leading "
        "factors taken explicitly, default 15), NProductExtraFactors, "
        "WynnDegree, VerifyConvergence (default True; a divergent product gives "
        "ComplexInfinity), AccuracyGoal, PrecisionGoal.");
    symtab_set_docstring("Product",
        "Product[f, {i, imax}]\n"
        "\tgives the product of f for i from 1 to imax.\n\n"
        "Product[f, {i, imin, imax}], Product[f, {i, imin, imax, di}] and "
        "Product[f, {i, {i1, i2, ...}}] use the standard iterator forms; "
        "multiple iterators give nested products (an inner bound may depend on "
        "an outer index). Product[f, i] gives the indefinite product "
        "(anti-quotient). The index is localised (HoldAll). Finite ranges are "
        "multiplied out directly; symbolic, indefinite and convergent infinite "
        "products are evaluated in exact closed form (n!, Pochhammer, Gamma "
        "ratios, base^k, QPochhammer, BarnesG) via a Method polyalgorithm.\n\n"
        "Options: Method (Automatic | \"Telescoping\" | \"Rational\" | "
        "\"Geometric\" | \"QProduct\"), VerifyConvergence (default True; a "
        "divergent infinite product gives Product::div), GenerateConditions, "
        "Assumptions. N[Product[...]] routes to NProduct.");
    symtab_set_docstring("Hyperfactorial",
        "Hyperfactorial[n]\n"
        "\tgives the hyperfactorial prod_{k=1}^{n} k^k.\n"
        "Exact (GMP) for a non-negative integer n; other orders stay symbolic. "
        "Listable, NumericFunction.");
    symtab_set_docstring("BarnesG",
        "BarnesG[z]\n"
        "\tgives the Barnes G-function.\n"
        "G(z+1) = Gamma[z] G(z) with G(1)=G(2)=1; for a positive integer n, "
        "G(n+1) = prod_{k=1}^{n-1} k! (exact via GMP), and G(m)=0 for "
        "non-positive integer m. Non-integer orders stay symbolic. Listable, "
        "NumericFunction.");
    symtab_set_docstring("QPochhammer",
        "QPochhammer[a, q, n]\n"
        "\tgives the q-Pochhammer symbol prod_{k=0}^{n-1} (1 - a q^k).\n"
        "QPochhammer[a, q] gives the infinite q-Pochhammer (a;q)_Inf for |q|<1. "
        "The finite form is exact/symbolic for a non-negative integer n; the "
        "infinite form evaluates for machine-real a, q. Listable, "
        "NumericFunction.");
    symtab_set_docstring("NIntegrate",
        "NIntegrate[f, {x, xmin, xmax}]\n"
        "\tgives a numerical approximation to the integral of f with respect "
        "to x from xmin to xmax.\n\n"
        "NIntegrate[f, {x, xmin, xmax}, {y, ymin, ymax}, ...] evaluates a "
        "multidimensional integral by adaptive cubature over a constant box, or "
        "iterated 1D quadrature when an inner bound depends on an outer "
        "variable. The variable is localised (HoldAll). "
        "xmin/xmax may be Infinity, -Infinity, or complex (a straight-line "
        "contour); extra nodes {x, x0, x1, ..., xk} give a piecewise-linear "
        "contour or mark interior singularities. Method -> Automatic chooses "
        "globally-adaptive Gauss-Kronrod for smooth finite integrands, "
        "double-exponential (tanh-sinh / sinh-sinh / exp-sinh) for endpoint "
        "singularities and infinite ranges and high precision, a Levin/zeros "
        "scheme for oscillatory integrands, an exponential endpoint map plus "
        "integration-between-the-zeros for an oscillatory endpoint singularity, "
        "and Monte-Carlo for high dimensions and region (Boole) integrands. "
        "Machine or arbitrary precision via WorkingPrecision.\n\n"
        "Options: Method (Automatic | GlobalAdaptive | GaussKronrodRule | "
        "DoubleExponential | TrapezoidalRule | LevinRule | "
        "OscillatorySingularity | MonteCarlo | "
        "QuasiMonteCarlo | AdaptiveMonteCarlo | PrincipalValue), "
        "WorkingPrecision (default MachinePrecision), PrecisionGoal, "
        "AccuracyGoal, MaxRecursion, MinRecursion, MaxPoints, Exclusions.");
    symtab_set_docstring("NRoots",
        "NRoots[lhs == rhs, var]\n"
        "\tyields a disjunction of equations var==r1 || var==r2 || ... giving "
        "numerical approximations to the roots of the polynomial equation in var. "
        "Roots of multiplicity k appear as k identical equations; a single root "
        "yields a bare equation. Real and complex coefficients are handled at "
        "machine and arbitrary precision. Method -> Automatic uses the Aberth-"
        "Ehrlich simultaneous iteration; \"CompanionMatrix\" uses companion-matrix "
        "eigenvalues (real QR directly, complex via a real 2n embedding); "
        "\"JenkinsTraub\" uses the three-stage Jenkins-Traub algorithm.\n\n"
        "Options: Method (Automatic | \"Aberth\" | \"CompanionMatrix\" | "
        "\"JenkinsTraub\"), PrecisionGoal (Automatic = machine; a digit count "
        "selects arbitrary precision), MaxIterations, StepMonitor.");
    symtab_set_docstring("NSolve",
        "NSolve[expr, vars]\n"
        "\tgives numerical approximations to the solutions of the equation or "
        "system expr for the variables vars, as a list of replacement-rule "
        "lists. NSolve[expr, vars, Reals] restricts to real solutions; the "
        "default domain is the complexes. vars may be a single variable or a "
        "list; NSolve[{e1, e2, ...}, vars] is the conjunction e1 && e2 && .... "
        "A working precision may be given as a trailing positional argument or "
        "via WorkingPrecision. Results: {} no solutions, {{x->s,...},...} the "
        "solutions (univariate roots are repeated by multiplicity), {{}} the "
        "universal solution. A univariate polynomial equation is solved with "
        "NRoots; square zero-dimensional polynomial systems use a Groebner-basis "
        "multiplication-matrix eigenvalue method (Method -> \"Symbolic\" uses "
        "lexicographic elimination); other equations fall back to Solve or "
        "FindRoot seeding. Integer, real, and complex coefficients are handled "
        "at machine and arbitrary precision.\n\n"
        "Options: MaxRoots, Method (Automatic | \"EndomorphismMatrix\" | "
        "\"Homotopy\" | \"Symbolic\"), WorkingPrecision, VerifySolutions, "
        "RandomSeeding.");
    symtab_set_docstring("Factorial",
        "n! or Factorial[n]\n"
        "\tgives the factorial of n.\n"
        "For non-negative integers, n! is computed exactly via GMP's mpz_fac_ui.\n"
        "For half-integers (n = m/2 with m odd) it reduces to Sqrt[Pi] times a\n"
        "rational from the Gamma functional equation. Negative integers give\n"
        "ComplexInfinity. Other inputs stay unevaluated.");
    symtab_set_docstring("Gamma",
        "Gamma[z]\n"
        "\tis the Euler gamma function Gamma(z).\n"
        "Gamma[a, z]\n"
        "\tis the upper incomplete gamma function Gamma(a, z).\n"
        "Gamma[a, z0, z1]\n"
        "\tis the generalized incomplete gamma Gamma(a, z0) - Gamma(a, z1).\n"
        "Integer and half-integer arguments reduce to exact values ((z-1)!,\n"
        "and rational multiples of Sqrt[Pi]); non-positive integers give\n"
        "ComplexInfinity. Machine and arbitrary-precision (MPFR) real inputs\n"
        "evaluate numerically, as do machine-precision complex inputs. Listable.");
    symtab_set_docstring("Beta",
        "Beta[a, b]\n"
        "\tis the Euler beta function B(a, b) = Gamma(a) Gamma(b) / Gamma(a+b).\n"
        "Beta[z, a, b]\n"
        "\tis the incomplete beta function Integral_0^z t^(a-1) (1-t)^(b-1) dt.\n"
        "Beta[z0, z1, a, b]\n"
        "\tis the generalized incomplete beta Beta[z1, a, b] - Beta[z0, a, b].\n"
        "Exact for rational arguments (a positive integer gives a rational via\n"
        "Pochhammer); non-positive integer poles give ComplexInfinity. Machine\n"
        "and arbitrary-precision (MPFR) real and complex inputs evaluate\n"
        "numerically. The incomplete form reduces through Hypergeometric2F1.\n"
        "Listable.");
    symtab_set_docstring("AiryAi",
        "AiryAi[z]\n"
        "\tgives the Airy function Ai(z), the solution of y'' = z y that decays\n"
        "\tas z -> +Infinity.\n"
        "AiryAi[0] = 1/(3^(2/3) Gamma[2/3]), AiryAi[+-Infinity] = 0. An entire\n"
        "function of z. Real and complex inputs evaluate numerically at machine\n"
        "or arbitrary (MPFR) precision; D[AiryAi[z], z] = AiryAiPrime[z]. Listable.");
    symtab_set_docstring("AiryAiPrime",
        "AiryAiPrime[z]\n"
        "\tgives the derivative Ai'(z) of the Airy function AiryAi.\n"
        "AiryAiPrime[0] = -1/(3^(1/3) Gamma[1/3]), AiryAiPrime[+Infinity] = 0. Real\n"
        "and complex inputs evaluate numerically at machine or arbitrary (MPFR)\n"
        "precision; D[AiryAiPrime[z], z] = z AiryAi[z]. Listable.");
    symtab_set_docstring("AiryBi",
        "AiryBi[z]\n"
        "\tgives the Airy function Bi(z), the solution of y'' = z y that grows\n"
        "\texponentially as z -> +Infinity.\n"
        "AiryBi[0] = 1/(3^(1/6) Gamma[2/3]), AiryBi[Infinity] = Infinity,\n"
        "AiryBi[-Infinity] = 0. An entire function of z. Real and complex inputs\n"
        "evaluate numerically at machine or arbitrary (MPFR) precision;\n"
        "D[AiryBi[z], z] = AiryBiPrime[z]. Listable.");
    symtab_set_docstring("AiryBiPrime",
        "AiryBiPrime[z]\n"
        "\tgives the derivative Bi'(z) of the Airy function AiryBi.\n"
        "AiryBiPrime[0] = 3^(1/6)/Gamma[1/3], AiryBiPrime[+Infinity] = Infinity. Real\n"
        "and complex inputs evaluate numerically at machine or arbitrary (MPFR)\n"
        "precision; D[AiryBiPrime[z], z] = z AiryBi[z]. Listable.");
    symtab_set_docstring("BesselJ",
        "BesselJ[n, z]\n"
        "\tgives the Bessel function of the first kind J_n(z), a solution of\n"
        "\tz^2 y'' + z y' + (z^2 - n^2) y = 0 regular at the origin.\n"
        "J_0(0) = 1, J_n(0) = 0 for integer n != 0. Has a branch cut along the\n"
        "negative real z axis for non-integer n. Real and complex order and\n"
        "argument evaluate numerically at machine or arbitrary (MPFR) precision;\n"
        "D[BesselJ[n, z], z] = (BesselJ[n-1, z] - BesselJ[n+1, z])/2. Listable.");

    symtab_set_docstring("BesselK",
        "BesselK[n, z]\n"
        "\tgives the modified Bessel function of the second kind K_n(z), a\n"
        "\tsolution of z^2 y'' + z y' - (z^2 + n^2) y = 0.\n"
        "K_n(z) is even in n (K_{-n} = K_n) and decays like e^{-z} as z -> Inf.\n"
        "K_0(0) = Infinity, K_n(0) = ComplexInfinity. Has a branch cut along the\n"
        "negative real z axis. Real and complex order and argument evaluate\n"
        "numerically at machine or arbitrary (MPFR) precision;\n"
        "D[BesselK[n, z], z] = -(BesselK[n-1, z] + BesselK[n+1, z])/2. Listable.");

    symtab_set_docstring("BesselI",
        "BesselI[n, z]\n"
        "\tgives the modified Bessel function of the first kind I_n(z), the\n"
        "\tsolution of z^2 y'' + z y' - (z^2 + n^2) y = 0 regular at the origin.\n"
        "I_0(0) = 1, I_n(0) = 0 for integer n != 0; I_n grows like e^z as\n"
        "z -> Inf and is even in n (I_{-n} = I_n). Has a branch cut along the\n"
        "negative real z axis for non-integer n. Real and complex order and\n"
        "argument evaluate numerically at machine or arbitrary (MPFR) precision;\n"
        "D[BesselI[n, z], z] = (BesselI[n-1, z] + BesselI[n+1, z])/2. Listable.");

    symtab_set_docstring("BesselY",
        "BesselY[n, z]\n"
        "\tgives the Bessel function of the second kind Y_n(z), the solution of\n"
        "\tz^2 y'' + z y' + (z^2 - n^2) y = 0 singular at the origin.\n"
        "Y_0(0) = -Infinity, Y_n(0) = ComplexInfinity for integer n != 0; Y_n has\n"
        "a logarithmic branch point at 0 and a branch cut along the negative real\n"
        "z axis, with Y_{-n} = (-1)^n Y_n for integer n. Real and complex order\n"
        "and argument evaluate numerically at machine or arbitrary (MPFR)\n"
        "precision; D[BesselY[n, z], z] = (BesselY[n-1, z] - BesselY[n+1, z])/2.\n"
        "Listable.");
    symtab_set_docstring("Erf",
        "Erf[z]\n"
        "\tgives the error function erf(z) = (2/Sqrt[Pi]) Integral_0^z e^(-t^2) dt.\n"
        "Erf[z0, z1]\n"
        "\tgives the generalized error function erf(z1) - erf(z0).\n"
        "Erf[0] = 0, Erf[Infinity] = 1, Erf[-Infinity] = -1. An entire function,\n"
        "odd in z. Real and complex inputs evaluate numerically at machine or\n"
        "arbitrary (MPFR) precision; D[Erf[z], z] = (2/Sqrt[Pi]) E^(-z^2). Listable.");
    symtab_set_docstring("Erfc",
        "Erfc[z]\n"
        "\tgives the complementary error function erfc(z) = 1 - erf(z).\n"
        "Erfc[0] = 1, Erfc[Infinity] = 0, Erfc[-Infinity] = 2. An entire\n"
        "function. Real inputs evaluate via libm/MPFR erfc (cancellation-free);\n"
        "complex inputs via 1 - erf(z) at machine or arbitrary (MPFR) precision.\n"
        "D[Erfc[z], z] = -(2/Sqrt[Pi]) E^(-z^2). Listable.");
    symtab_set_docstring("Erfi",
        "Erfi[z]\n"
        "\tgives the imaginary error function erfi(z) = -I Erf[I z]\n"
        "\t= (2/Sqrt[Pi]) Integral_0^z e^(t^2) dt.\n"
        "Erfi[0] = 0, Erfi[Infinity] = Infinity, Erfi[I Infinity] = I. An entire\n"
        "function, odd in z. Real and complex inputs evaluate numerically at\n"
        "machine or arbitrary (MPFR) precision; D[Erfi[z], z] = (2/Sqrt[Pi]) E^(z^2).\n"
        "Listable.");
    symtab_set_docstring("ExpIntegralEi",
        "ExpIntegralEi[z]\n"
        "\tgives the exponential integral Ei(z), the principal value of\n"
        "\t-Integral_{-z}^Infinity e^-t/t dt, with a branch cut on (-Infinity, 0).\n"
        "ExpIntegralEi[0] = -Infinity, ExpIntegralEi[Infinity] = Infinity,\n"
        "ExpIntegralEi[-Infinity] = 0, ExpIntegralEi[+-I Infinity] = +-I Pi. Real and\n"
        "complex inputs evaluate numerically at machine or arbitrary (MPFR) precision;\n"
        "D[ExpIntegralEi[z], z] = E^z/z. Listable.");
    symtab_set_docstring("LogIntegral",
        "LogIntegral[z]\n"
        "\tgives the logarithmic integral li(z), the principal value of\n"
        "\tIntegral_0^z dt/ln t, equal to ExpIntegralEi[Log[z]], with a branch cut\n"
        "\ton (-Infinity, 1). LogIntegral[0] = 0, LogIntegral[1] = -Infinity,\n"
        "LogIntegral[Infinity] = Infinity. Real and complex inputs evaluate\n"
        "numerically at machine or arbitrary (MPFR) precision; D[LogIntegral[z], z] =\n"
        "1/Log[z]. Listable.");
    symtab_set_docstring("InverseErf",
        "InverseErf[s]\n"
        "\tgives the inverse error function: the z solving s = Erf[z].\n"
        "InverseErf[z0, s]\n"
        "\tgives the inverse of the generalized error function Erf[z0, z].\n"
        "InverseErf[0] = 0, InverseErf[1] = Infinity, InverseErf[-1] = -Infinity.\n"
        "Odd in s. Numerical values are given only for real s in [-1, 1], at\n"
        "machine or arbitrary (MPFR) precision; D[InverseErf[z], z] =\n"
        "(Sqrt[Pi]/2) E^(InverseErf[z]^2). Listable.");
    symtab_set_docstring("InverseErfc",
        "InverseErfc[s]\n"
        "\tgives the inverse complementary error function: the z solving s = Erfc[z].\n"
        "InverseErfc[0] = Infinity, InverseErfc[1] = 0, InverseErfc[2] = -Infinity.\n"
        "Numerical values are given only for real s in [0, 2], at machine or\n"
        "arbitrary (MPFR) precision; D[InverseErfc[z], z] =\n"
        "-(Sqrt[Pi]/2) E^(InverseErfc[z]^2). Listable.");
    symtab_set_docstring("PolyGamma",
        "PolyGamma[z]\n"
        "\tgives the digamma function psi(z) (rewritten as PolyGamma[0, z]).\n"
        "PolyGamma[n, z]\n"
        "\tgives the n-th derivative of the digamma function, psi^(n)(z).\n"
        "Positive-integer arguments reduce to exact values: psi(m) to a rational\n"
        "minus EulerGamma, and psi^(n)(m) for odd n to a rational plus a rational\n"
        "multiple of Pi^(n+1); even orders stay symbolic. Non-positive integer\n"
        "arguments give ComplexInfinity. Inexact real and complex arguments evaluate\n"
        "numerically at machine or arbitrary (MPFR) precision. PolyGamma[-1, z] gives\n"
        "LogGamma[z]. Listable.");
    symtab_set_docstring("HarmonicNumber",
        "HarmonicNumber[n]\n"
        "\tgives the n-th harmonic number H_n = Sum_{i=1}^n 1/i.\n"
        "HarmonicNumber[n, r]\n"
        "\tgives the order-r harmonic number H_n^(r) = Sum_{i=1}^n 1/i^r.\n"
        "Non-negative integer n expands to the exact finite sum (a rational for\n"
        "integer r, an explicit sum for symbolic r); HarmonicNumber[Infinity, r] is\n"
        "Zeta[r]; a non-positive integer order r gives the Faulhaber polynomial in n.\n"
        "Inexact arguments evaluate numerically at machine or arbitrary (MPFR)\n"
        "precision, including complex order, via Zeta[r] - Zeta[r, n+1] (and the\n"
        "digamma form for r = 1). Listable.");
    symtab_set_docstring("LogGamma",
        "LogGamma[z]\n"
        "\tgives the log-gamma function log(Gamma(z)), analytic except for a branch\n"
        "cut on the negative reals. Exact at integer and half-integer z (with the\n"
        "negative-axis branch term), divergent (Infinity) at non-positive integers,\n"
        "and evaluated numerically for real and complex z at machine or arbitrary\n"
        "(MPFR) precision. D[LogGamma[z], z] is PolyGamma[0, z]. Listable.");
    symtab_set_docstring("Pochhammer",
        "Pochhammer[a, n]\n"
        "\tgives the Pochhammer symbol (a)_n = a (a+1) ... (a+n-1) = Gamma[a+n]/Gamma[a].\n"
        "Exact integer n expands to the product of n linear factors: a polynomial\n"
        "product for symbolic a, an exact Integer/Rational for numeric a (negative n\n"
        "gives the reciprocal product). Other numeric arguments evaluate via the\n"
        "Gamma ratio, reducing exact half-integers to rational multiples of Sqrt[Pi]\n"
        "and tracking machine or arbitrary (MPFR) precision; machine-precision\n"
        "complex arguments evaluate too. Other arguments stay symbolic. Listable.");
    symtab_set_docstring("Zeta",
        "Zeta[s]\n"
        "\tis the Riemann zeta function zeta(s) = Sum_{k>=1} k^-s.\n"
        "Zeta[s, a]\n"
        "\tis the Hurwitz zeta function zeta(s, a) = Sum_{k>=0} (k+a)^-s.\n"
        "Even positive integers give rational multiples of Pi^(2n), negative\n"
        "integers give rationals, Zeta[0] is -1/2, and Zeta[1] is ComplexInfinity;\n"
        "odd positive integers stay symbolic. Hurwitz zeta at a positive integer a\n"
        "reduces to Zeta[s] minus a finite power sum, and Zeta[s, 1/2] is\n"
        "(2^s - 1) Zeta[s]. Real, complex, machine and\n"
        "arbitrary-precision (MPFR) numeric arguments evaluate numerically via\n"
        "mpfr_zeta (real Riemann) or an Euler-Maclaurin kernel. Listable.");
    symtab_set_docstring("HurwitzZeta",
        "HurwitzZeta[s, a]\n"
        "\tis the Hurwitz zeta function zeta(s, a) = Sum_{k>=0} (k + a)^-s.\n"
        "Identical to Zeta[s, a] for Re(a) > 0, but built on the principal-branch\n"
        "power (k + a)^-s, so it differs from Zeta for non-positive real a and has\n"
        "poles at a = 0, -1, -2, ... . HurwitzZeta[s, 1] is Zeta[s], "
        "HurwitzZeta[s, 1/2]\n"
        "is (2^s - 1) Zeta[s], and a positive integer a reduces to Zeta[s] minus a\n"
        "finite power sum. A non-positive integer a gives ComplexInfinity for "
        "positive\n"
        "integer s and the Bernoulli-polynomial value for non-positive integer s.\n"
        "Real, complex, machine and arbitrary-precision (MPFR) arguments evaluate\n"
        "numerically via an Euler-Maclaurin kernel. Listable.");
    symtab_set_docstring("BernoulliB",
        "BernoulliB[n]\n"
        "\tgives the Bernoulli number B_n.\n"
        "BernoulliB[n, x]\n"
        "\tgives the Bernoulli polynomial B_n(x).\n"
        "Non-negative integer n gives the exact rational B_n (odd n > 1 give 0,\n"
        "B_0 = 1, B_1 = -1/2); an inexact integer-valued n evaluates it at machine\n"
        "or arbitrary (MPFR) precision. BernoulliB[n, x] expands the degree-n\n"
        "polynomial with exact rational coefficients, staying symbolic in x or\n"
        "evaluating numerically when x is inexact. Listable.");
    symtab_set_docstring("EulerE",
        "EulerE[n]\n"
        "\tgives the Euler number E_n.\n"
        "EulerE[n, x]\n"
        "\tgives the Euler polynomial E_n(x).\n"
        "Non-negative integer n gives the exact integer E_n (odd n give 0,\n"
        "E_0 = 1, E_2 = -1, E_4 = 5); an inexact integer-valued n evaluates it at\n"
        "machine or arbitrary (MPFR) precision. EulerE[n, x] expands the degree-n\n"
        "polynomial with exact rational coefficients, staying symbolic in x or\n"
        "evaluating numerically when x is inexact; EulerE[n, 1/2] folds to\n"
        "2^-n EulerE[n]. Listable.");
    symtab_set_docstring("PolyLog",
        "PolyLog[n, z]\n"
        "\tgives the polylogarithm Li_n(z) = Sum_{k>=1} z^k/k^n.\n"
        "PolyLog[n, p, z]\n"
        "\tgives the Nielsen generalized polylogarithm S_{n,p}(z) (accepted but\n"
        "left unevaluated).\n"
        "Special arguments reduce in closed form: PolyLog[n, 0] = 0,\n"
        "PolyLog[1, z] = -Log[1-z], PolyLog[0, z] = z/(1-z), negative integer\n"
        "orders give rational functions, PolyLog[n, 1] = Zeta[n] and\n"
        "PolyLog[n, -1] = (2^(1-n)-1) Zeta[n] for integer n >= 2, with exact forms\n"
        "for PolyLog[2, 1/2] and PolyLog[3, 1/2]. Inexact real or complex arguments\n"
        "evaluate numerically at machine or arbitrary (MPFR) precision via a power\n"
        "series or the Jonquiere/zeta expansion. There is a branch cut from 1 to\n"
        "Infinity. Listable.");
    symtab_set_docstring("LerchPhi",
        "LerchPhi[z, s, a]\n"
        "\tis the Lerch transcendent Phi(z, s, a) = Sum_{k>=0} z^k/(k + a)^s.\n"
        "It generalizes Zeta, HurwitzZeta and PolyLog: LerchPhi[1, s, a] is\n"
        "Zeta[s, a] and z LerchPhi[z, s, 1] is PolyLog[s, z]. Exact reductions\n"
        "cover z = 0 (a^-s), s = 0 (1/(1-z)), z = +-1, positive integer a (a\n"
        "PolyLog form) and negative integer s (a rational function of z). The\n"
        "options DoublyInfinite -> True (sum k from -Infinity to Infinity) and\n"
        "IncludeSingularTerm -> True (keep the k + a = 0 term) are supported.\n"
        "Inexact arguments with |z| < 1 evaluate numerically at machine or\n"
        "arbitrary (MPFR) precision; |z| > 1 stays symbolic. Listable.");
    symtab_set_docstring("ProductLog",
        "ProductLog[z]\n"
        "\tgives the principal solution w of z == w e^w (the Lambert W function).\n"
        "ProductLog[k, z] gives the k-th solution (k any integer, k == 0 the\n"
        "principal branch); branches are ordered by imaginary part. ProductLog[z]\n"
        "is real for z >= -1/e and has a branch cut along (-Infinity, -1/e].\n"
        "Exact values include ProductLog[0] = 0, ProductLog[-1/E] = -1,\n"
        "ProductLog[E] = 1, ProductLog[-Pi/2] = I Pi/2 and ProductLog[k, 0] =\n"
        "-Infinity for k != 0. Inexact real or complex arguments evaluate\n"
        "numerically at machine or arbitrary (MPFR) precision. Satisfies\n"
        "D[ProductLog[z], z] = ProductLog[z]/(z (1 + ProductLog[z])). Listable.");
    symtab_set_docstring("StieltjesGamma",
        "StieltjesGamma[n]\n"
        "\tgives the n-th Stieltjes constant gamma_n, the Laurent coefficients of\n"
        "Zeta about s = 1. StieltjesGamma[0] is EulerGamma; higher constants are\n"
        "inert (they stay symbolic) and appear in Series expansions of Zeta. Listable.");
    symtab_set_docstring("EulerGamma",
        "EulerGamma\n"
        "\tis Euler's constant gamma, with numerical value ~= 0.5772156649.\n"
        "EulerGamma is the Euler-Mascheroni constant, the limit of\n"
        "HarmonicNumber[n] - Log[n] as n -> Infinity. It is a mathematical\n"
        "constant: it has attributes Constant and Protected, NumericQ[EulerGamma]\n"
        "is True, and D[EulerGamma, x] is 0. N[EulerGamma, prec] evaluates it to\n"
        "any precision.");
    symtab_set_docstring("Pi",
        "Pi\n"
        "\tis pi, with numerical value ~= 3.14159.\n"
        "Pi is a mathematical constant: it has attributes Constant and Protected,\n"
        "NumericQ[Pi] is True, and D[Pi, x] is 0. N[Pi, prec] evaluates it to any\n"
        "precision.");
    symtab_set_docstring("E",
        "E\n"
        "\tis the exponential constant e (base of natural logarithms), with\n"
        "\tnumerical value ~= 2.71828.\n"
        "E is a mathematical constant: it has attributes Constant and Protected,\n"
        "NumericQ[E] is True, and D[E, x] is 0. N[E, prec] evaluates it to any\n"
        "precision.");
    symtab_set_docstring("Catalan",
        "Catalan\n"
        "\tis Catalan's constant, with numerical value ~= 0.915966.\n"
        "Catalan is the sum over k >= 0 of (-1)^k (2 k + 1)^-2. It is a\n"
        "mathematical constant: it has attributes Constant and Protected,\n"
        "NumericQ[Catalan] is True, and D[Catalan, x] is 0. N[Catalan, prec]\n"
        "evaluates it to any precision.");
    symtab_set_docstring("Degree",
        "Degree\n"
        "\tgives the number of radians in one degree, with numerical value\n"
        "\tPi/180 (~= 0.0174533).\n"
        "Multiply by Degree to convert degrees to radians, so 30 Degree is 30\n"
        "degrees. It is a mathematical constant: it has attributes Constant and\n"
        "Protected, NumericQ[Degree] is True, and D[Degree, x] is 0. N[Degree,\n"
        "prec] evaluates it to any precision.");
    symtab_set_docstring("GoldenRatio",
        "GoldenRatio\n"
        "\tis the golden ratio phi = (1 + Sqrt[5])/2, with numerical value\n"
        "\t~= 1.61803.\n"
        "GoldenRatio is the positive root of x^2 == x + 1. It is a mathematical\n"
        "constant: it has attributes Constant and Protected, NumericQ[GoldenRatio]\n"
        "is True, and D[GoldenRatio, x] is 0. N[GoldenRatio, prec] evaluates it to\n"
        "any precision.");
    symtab_set_docstring("GoldenAngle",
        "GoldenAngle\n"
        "\tis the golden angle (3 - Sqrt[5]) Pi = 2 Pi / GoldenRatio^2, with\n"
        "\tnumerical value ~= 2.39996 radians (~= 137.5 degrees).\n"
        "GoldenAngle is a mathematical constant: it has attributes Constant and\n"
        "Protected, NumericQ[GoldenAngle] is True, and D[GoldenAngle, x] is 0.\n"
        "N[GoldenAngle, prec] evaluates it to any precision.");
    symtab_set_docstring("Glaisher",
        "Glaisher\n"
        "\tis the Glaisher-Kinkelin constant A, with numerical value ~= 1.28243.\n"
        "Glaisher's constant satisfies Log[A] == 1/12 - Zeta'[-1], where Zeta is\n"
        "the Riemann zeta function. It is a mathematical constant: it has\n"
        "attributes Constant and Protected, NumericQ[Glaisher] is True, and\n"
        "D[Glaisher, x] is 0. N[Glaisher, prec] evaluates it to any precision.");
    symtab_set_docstring("Khinchin",
        "Khinchin\n"
        "\tis Khinchin's constant K (also Khintchine's constant), with numerical\n"
        "\tvalue ~= 2.68545.\n"
        "Khinchin's constant is the limiting geometric mean of the partial\n"
        "quotients in the continued-fraction expansion of almost every real\n"
        "number, given by the product over s >= 1 of (1 + 1/(s (s + 2)))^Log2[s].\n"
        "It is a mathematical constant: it has attributes Constant and Protected,\n"
        "NumericQ[Khinchin] is True, and D[Khinchin, x] is 0. N[Khinchin, prec]\n"
        "evaluates it to any precision.");
    symtab_set_docstring("Factorial2", "Factorial2[n] (also typeset n!!) gives the double factorial of n.\nFor non-negative integer n: n!! = n * (n-2) * (n-4) * ... down to 2 (n even) or 1 (n odd).\nSpecial values: 0!! = 1, (-1)!! = 1.\nNegative even integers and negative odd integers below -1 give ComplexInfinity.\nFactorial2 stays unevaluated on symbolic arguments.");
    symtab_set_docstring("Fibonacci",
        "Fibonacci[n]\n"
        "\tgives the nth Fibonacci number F_n.\n"
        "Fibonacci[n, x]\n"
        "\tgives the nth Fibonacci polynomial F_n(x).\n"
        "Exact integer orders are computed via GMP fast doubling (numbers) or the\n"
        "recurrence F_k = x F_{k-1} + F_{k-2} (polynomials); negative orders use\n"
        "F_{-n} = (-1)^(n+1) F_n. Inexact or complex orders evaluate the\n"
        "generalized closed form numerically. Listable; symbolic orders stay\n"
        "unevaluated.");
    symtab_set_docstring("LucasL",
        "LucasL[n]\n"
        "\tgives the nth Lucas number L_n.\n"
        "LucasL[n, x]\n"
        "\tgives the nth Lucas polynomial L_n(x).\n"
        "Exact integer orders are computed via GMP fast doubling (numbers, using\n"
        "L_m = 2 F_{m+1} - F_m) or the recurrence L_k = x L_{k-1} + L_{k-2} with\n"
        "L_0 = 2, L_1 = x (polynomials); negative orders use L_{-n} = (-1)^n L_n.\n"
        "Inexact or complex orders evaluate the generalized closed form\n"
        "L_n = phi^n + Cos[Pi n] phi^-n (phi = GoldenRatio) numerically.\n"
        "Listable; symbolic orders stay unevaluated.");
    symtab_set_docstring("Binomial",
        "Binomial[n, m]\n"
        "\tgives the binomial coefficient C(n, m) = n! / (m! (n-m)!).\n"
        "For non-negative integer arguments, computed exactly via GMP's\n"
        "mpz_bin_uiui. Generalised forms (negative or symbolic n, half-integer\n"
        "m) reduce through the Gamma functional equation; non-decidable forms\n"
        "stay unevaluated.");

    symtab_set_docstring("HypergeometricPFQ",
        "HypergeometricPFQ[{a1, ...}, {b1, ...}, z]\n"
        "\tis the generalized hypergeometric function pFq(a;b;z), the series\n"
        "Sum (prod_i Pochhammer[a_i, k] / prod_j Pochhammer[b_j, k]) z^k / k!.\n"
        "Common upper/lower parameters cancel; a non-positive integer upper\n"
        "parameter terminates the series to a polynomial. Evaluates to machine,\n"
        "arbitrary-precision (MPFR), and complex numbers by direct summation in\n"
        "the convergent regime (p<=q for all z; p==q+1 for |z|<1).");

    symtab_set_docstring("Hypergeometric0F1",
        "Hypergeometric0F1[b, z]\n"
        "\tis the confluent hypergeometric 0F1, equal to "
        "HypergeometricPFQ[{}, {b}, z].");

    symtab_set_docstring("Hypergeometric1F1",
        "Hypergeometric1F1[a, b, z]\n"
        "\tis Kummer's confluent hypergeometric 1F1, equal to "
        "HypergeometricPFQ[{a}, {b}, z].");

    symtab_set_docstring("Hypergeometric2F1",
        "Hypergeometric2F1[a, b, c, z]\n"
        "\tis the Gauss hypergeometric 2F1, equal to "
        "HypergeometricPFQ[{a, b}, {c}, z].");

    // Structural Manipulation
    symtab_set_docstring("List",
        "{e1, e2, ...} or List[e1, e2, ...]\n"
        "\trepresents an ordered list of the elements ei.\n"
        "List is the fundamental container head: vectors are lists, matrices are\n"
        "lists of lists, and the structural operators (Part, Map, Take, Drop,\n"
        "Length, ...) act on List. Elements are evaluated normally and kept in the\n"
        "given order (List has no Orderless attribute). The parser writes the {...}\n"
        "syntax to List, and the printer renders List[...] back as {...}.");
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
    symtab_set_docstring("Take",
        "Take[list, n]\n"
        "\tgives the first n elements of list.\n"
        "Take[list, -n]\n"
        "\tgives the last n elements.\n"
        "Take[list, {m, n}]\n"
        "\tgives elements m through n.\n"
        "Take[list, {m, n, s}]\n"
        "\tgives elements m through n in steps of s.\n"
        "Take[list, {m}]\n"
        "\tgives the single element at position m (wrapped in the head of list).\n"
        "Take[list, spec1, spec2, ...]\n"
        "\ttakes elements at successive levels, e.g. a sub-block of a matrix.\n\n"
        "Negative indices count from the end; UpTo[n], All, and None are also accepted\n"
        "as specifications. Indices are 1-based; out-of-range requests leave the\n"
        "expression unevaluated. Take operates on any expression, not just List.");
    symtab_set_docstring("Drop",
        "Drop[list, n]\n"
        "\tgives list with its first n elements dropped.\n"
        "Drop[list, -n]\n"
        "\tdrops the last n elements.\n"
        "Drop[list, {m, n}]\n"
        "\tdrops elements m through n.\n"
        "Drop[list, {m, n, s}]\n"
        "\tdrops elements m through n in steps of s.\n"
        "Drop[list, {m}]\n"
        "\tdrops the single element at position m.\n"
        "Drop[list, spec1, spec2, ...]\n"
        "\tdrops elements at successive levels.\n\n"
        "Negative indices count from the end; UpTo[n], All, and None are also accepted.\n"
        "Indices are 1-based; out-of-range requests leave the expression unevaluated.");


    // Linear Algebra
    symtab_set_docstring("Dot",
        "a . b . c or Dot[a, b, c]\n"
        "\tcontracts the last index of each argument with the first index of\n"
        "\tthe next: matrix-matrix, matrix-vector, vector-vector, and general\n"
        "\ttensor inner products.\n"
        "Numeric machine-precision Real / Complex matrix-matrix dot dispatches\n"
        "to BLAS dgemm / zgemm when available; exact and symbolic inputs use\n"
        "the elementwise sum-of-products.");
    symtab_set_docstring("Det",
        "Det[m]\n"
        "\tgives the determinant of the square matrix m.\n"
        "Exact integer / rational / symbolic inputs use Bareiss-style\n"
        "fraction-free Gaussian elimination; machine-precision Real / Complex\n"
        "inputs dispatch to LAPACK LU (dgetrf / zgetrf) and accumulate the\n"
        "pivot-signed product of diagonal entries; arbitrary-precision MPFR\n"
        "inputs run a Doolittle LU at the input precision.");
    symtab_set_docstring("Cross",
        "Cross[a, b]\n"
        "\tgives the vector cross product of two length-3 vectors.\n"
        "Cross[a1, a2, ..., a(n-1)]\n"
        "\tgives the generalized (n-1)-fold cross product in n dimensions,\n"
        "\ti.e. the unique vector orthogonal to all inputs whose components\n"
        "\tare the signed cofactor minors of the matrix [a1; a2; ...; en].");
    symtab_set_docstring("Norm",
        "Norm[expr]\n"
        "\tgives the 2-norm of a number, vector, or matrix (Frobenius norm for\n"
        "\tmatrices).\n"
        "Norm[expr, p]\n"
        "\tgives the p-norm: Abs[expr] for scalars; (Sum |xi|^p)^(1/p) for\n"
        "\tvectors with 1 <= p < Infinity; Max[Abs[expr]] for p == Infinity;\n"
        "\tinduced operator norm for matrices when p is 1, 2, or Infinity.");
    symtab_set_docstring("Normalize",
        "Normalize[v]\n"
        "\tgives the normalized form of a vector v (effectively v / Norm[v]).\n"
        "Normalize[z]\n"
        "\tgives the normalized form of a scalar (incl. complex) z, namely z / Abs[z].\n"
        "Normalize[expr, f]\n"
        "\tnormalizes with respect to the norm function f, i.e. expr / f[expr].\n"
        "Zero input is returned unchanged.");
    symtab_set_docstring("Tr",
        "Tr[m]\n"
        "\tgives the trace of the matrix m, i.e. the sum of its diagonal\n"
        "\tentries (for a rectangular m, sums entries m[[i, i]] up to\n"
        "\tMin[Dimensions[m]]).\n"
        "Tr[m, f]\n"
        "\tcombines the diagonal entries with f instead of Plus.\n"
        "Tr[m, f, n]\n"
        "\twalks down to level n, summing the multi-index diagonal of a\n"
        "\trank-n tensor.");
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
    symtab_set_docstring("LatticeReduce",
        "LatticeReduce[m]\n"
        "\tgives an LLL-reduced basis for the lattice spanned by the rows\n"
        "\t(vectors) of m.  The entries of m may be integers, Gaussian\n"
        "\tintegers, rationals, or Gaussian rationals.  Reduction is exact\n"
        "\t(GMP rational arithmetic, so it is correct for both machine-size\n"
        "\tand arbitrary-precision entries) and preserves the lattice, its\n"
        "\tdeterminant, and every linear relation among the rows.  The rows\n"
        "\tmust be linearly independent.");
    symtab_set_docstring("FindIntegerNullVector",
        "FindIntegerNullVector[{x1, ..., xn}]\n"
        "\tfinds integers {a1, ..., an}, not all zero, with "
        "a1 x1 + ... + an xn == 0 (PSLQ / integer-relation detection).\n"
        "FindIntegerNullVector[{x1, ..., xn}, d]\n"
        "\trestricts the search to relations of norm <= d.\n"
        "The xi may be real or complex, exact or inexact; for complex xi the "
        "ai are Gaussian integers.  Exact relations are validated with "
        "PossibleZeroQ; for inexact xi the relation holds to the precision of "
        "the input.  When no relation is found the call is returned "
        "unevaluated.\n"
        "Options:\n"
        "\tWorkingPrecision    Automatic, or a digit count for the search.\n"
        "\tZeroTest            Automatic, or a function applied to the residual.");
    symtab_set_docstring("IdentityMatrix", "IdentityMatrix[n] gives the n x n identity matrix.\nIdentityMatrix[{m, n}] gives the m x n identity matrix.");
    symtab_set_docstring("Fit", "Fit[data, {f1, ..., fn}, vars] finds a least-squares fit a1 f1 + ... + an fn to data for functions of the variables vars (a symbol or list of symbols).\nFit[{m, v}] gives the coefficient vector minimizing ||m.a - v|| for design matrix m and response vector v.\nData may be a list of values {v1, ...} (coordinates 1, 2, ...), univariate pairs {{x, v}, ...}, or multivariate rows {{x, y, ..., v}, ...}.\nOptions: WorkingPrecision (Automatic | n | Infinity), FitRegularization ({\"Tikhonov\"|\"L2\"|\"RidgeRegression\"|\"LASSO\"|\"L1\", lambda}), NormFunction.");
    symtab_set_docstring("DesignMatrix", "DesignMatrix[data, {f1, ..., fn}, vars] gives the design matrix with entries f_i evaluated at the data coordinates.\nData shapes match Fit. The WorkingPrecision option converts entries to machine or n-digit reals; otherwise they are exact.");
    symtab_set_docstring("DiagonalMatrix", "DiagonalMatrix[list] gives a matrix with the elements of list on the leading diagonal, and zero elsewhere.\nDiagonalMatrix[list, k] gives a matrix with the elements of list on the k-th diagonal.\nDiagonalMatrix[list, k, n] pads with zeros to create an n x n matrix.");
    symtab_set_docstring("HilbertMatrix", "HilbertMatrix[n] gives the n x n Hilbert matrix with entries 1/(i + j - 1).\nHilbertMatrix[{m, n}] gives the m x n Hilbert matrix.\nEntries are exact Rationals unless the WorkingPrecision option requests MachinePrecision or a digit count.");
    symtab_set_docstring("HankelMatrix", "HankelMatrix[n] gives the n x n Hankel matrix with first row and column the integers 1..n.\nHankelMatrix[{c1, ..., cm}] gives the m x m Hankel matrix with first column the given list.\nHankelMatrix[{c1, ..., cm}, {r1, ..., rn}] gives the m x n Hankel matrix with first column the first list and last row the second.\nA Hankel matrix is constant along its antidiagonals; entries are copied verbatim.");
    symtab_set_docstring("ToeplitzMatrix", "ToeplitzMatrix[n] gives the n x n Toeplitz matrix with first row and column the integers 1..n.\nToeplitzMatrix[{c1, ..., cn}] gives the n x n symmetric Toeplitz matrix with first column the given list.\nToeplitzMatrix[{c1, ..., cm}, {r1, ..., rn}] gives the m x n Toeplitz matrix with first column the first list and first row the second.\nA Toeplitz matrix is constant along its diagonals; entries are copied verbatim.");
    symtab_set_docstring("VandermondeMatrix", "VandermondeMatrix[{x1, ..., xn}] gives the n x n Vandermonde matrix with entry (i, j) equal to xi^(j-1).\nVandermondeMatrix[{x1, ..., xn}, k] gives the n x k Vandermonde matrix.\nThe nodes need not be numerical or distinct; columns are successive powers, so the first column is all ones.");
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
    symtab_set_docstring("QRDecomposition",
        "QRDecomposition[m]\n"
        "\tgives the QR decomposition of m as a list {q, r}, where q is\n"
        "\trow-orthonormal (row-unitary in the complex case) and r is\n"
        "\tupper triangular.  The original matrix satisfies\n"
        "\tm == ConjugateTranspose[q] . r.\n"
        "\n"
        "QRDecomposition computes the \"thin\" QR factorisation: when m\n"
        "has rank r, both q and r have r rows.  For an n x p input,\n"
        "q has dimensions r x n and r has dimensions r x p, so q's\n"
        "rows live in the column space of m and r encodes the original\n"
        "columns in that basis.\n"
        "\n"
        "QRDecomposition[m, Pivoting -> True]\n"
        "\tgives a list {q, r, p} where p is a p x p permutation matrix\n"
        "\tsuch that m . p == ConjugateTranspose[q] . r.  With pivoting\n"
        "\tthe diagonal of r appears in order of decreasing magnitude.\n"
        "\n"
        "QRDecomposition works on every input family supported by the\n"
        "rest of the linear-algebra builtins:\n"
        "\t- exact integer / rational matrices (output stays exact,\n"
        "\t  with Sqrt[...] in the column norms)\n"
        "\t- complex matrices (q's rows are unitary in the Hermitian\n"
        "\t  inner product)\n"
        "\t- machine-precision Real matrices (output is Real at machine\n"
        "\t  precision, matching the inexact-in / inexact-out contract)\n"
        "\t- arbitrary-precision MPFR matrices (output at the input\n"
        "\t  precision)\n"
        "\t- free-symbolic matrices (output in closed symbolic form)\n"
        "\n"
        "The algorithm is Modified Gram-Schmidt on the columns of m,\n"
        "applied through the evaluator so symbolic, exact, and inexact\n"
        "inputs share one code path.  Rank-deficient inputs (columns in\n"
        "the span of earlier columns) produce a shorter q / r without\n"
        "any error.\n"
        "\n"
        "A non-rank-2 or empty matrix emits QRDecomposition::matrix and\n"
        "the call is left unevaluated.  Unknown option keys or values\n"
        "emit QRDecomposition::opts and the call is left unevaluated.\n"
        "TargetStructure -> \"Structured\" is reserved for a future\n"
        "release and currently leaves the call unevaluated.");
    symtab_set_docstring("LUDecomposition",
        "LUDecomposition[m]\n"
        "\tgives the LU decomposition of a square matrix m as a list\n"
        "\t{lu, p, c}.  The first element lu is the combined Doolittle\n"
        "\tfactor matrix: its strictly-lower triangle is L (with an\n"
        "\timplicit unit diagonal) and its upper triangle is U.  The\n"
        "\tsecond element p is a 1-indexed row-permutation vector such\n"
        "\tthat m[[p]] == l . u where l = LowerTriangularize[lu, -1] +\n"
        "\tIdentityMatrix[n] and u = UpperTriangularize[lu].  The third\n"
        "\telement c is an L-infinity condition-number estimate for\n"
        "\tapproximate numerical matrices, or 0 for exact / symbolic m.\n"
        "\n"
        "LUDecomposition works on every input family supported by the\n"
        "rest of the linear-algebra builtins:\n"
        "\t- exact integer / rational matrices (output stays exact)\n"
        "\t- complex matrices\n"
        "\t- machine-precision Real matrices (LAPACK dgetrf / zgetrf\n"
        "\t  with dgecon / zgecon for the condition estimate)\n"
        "\t- arbitrary-precision MPFR matrices (Householder-free\n"
        "\t  Doolittle at the input precision; condition estimate via\n"
        "\t  the explicit inverse)\n"
        "\t- free-symbolic matrices (output in closed symbolic form)\n"
        "\n"
        "The algorithm is Doolittle's elimination with partial row\n"
        "pivoting.  Numerical inputs use largest-|pivot| selection;\n"
        "symbolic / exact inputs advance to the next non-zero pivot\n"
        "only when the natural choice is provably zero.\n"
        "\n"
        "A singular m emits LUDecomposition::sing and the factorisation\n"
        "completes with a zero on U's diagonal at the singular step.\n"
        "\n"
        "A non-square or empty matrix emits LUDecomposition::matsq and\n"
        "the call is left unevaluated.");
    symtab_set_docstring("SingularValueDecomposition",
        "SingularValueDecomposition[m]\n"
        "\tgives the singular value decomposition of a matrix m as a\n"
        "\tlist {u, sigma, v}, where sigma is a diagonal matrix and\n"
        "\tm == u . sigma . ConjugateTranspose[v].  u and v have\n"
        "\torthonormal columns.\n"
        "\n"
        "SingularValueDecomposition[m, k]\n"
        "\tgives the SVD associated with the k largest singular values\n"
        "\t(or |k| smallest when k is negative).\n"
        "SingularValueDecomposition[m, UpTo[k]]\n"
        "\tgives the SVD for as many of the k largest singular values\n"
        "\tas are available (up to MatrixRank[m]).\n"
        "\n"
        "SingularValueDecomposition[{m, a}]\n"
        "\tgives the generalized singular value decomposition of m\n"
        "\twith respect to a as {{u, ua}, {sigma, sigma_a}, v} such\n"
        "\tthat m == u . sigma . ConjugateTranspose[v] and\n"
        "\ta == ua . sigma_a . ConjugateTranspose[v].  Uses LAPACK\n"
        "\tdggsvd3 / zggsvd3 for real / complex machine-precision\n"
        "\tinputs; exact-numeric input is numericalised to 53 bits and\n"
        "\trouted through the same path.  High-precision MPFR input is\n"
        "\tcurrently downgraded to machine precision and emits the\n"
        "\t::gmpdwn warning.  Free-symbolic input emits ::nogsymb and\n"
        "\tleaves the call unevaluated.\n"
        "\n"
        "SingularValueDecomposition works on every input family supported\n"
        "by the rest of the linear-algebra builtins:\n"
        "\t- exact integer / rational matrices (output stays exact;\n"
        "\t  singular values are Sqrt[...] forms when irrational)\n"
        "\t- complex matrices (u and v are unitary in the Hermitian\n"
        "\t  inner product)\n"
        "\t- machine-precision Real matrices (LAPACK dgesdd / zgesdd,\n"
        "\t  or dggsvd3 / zggsvd3 for the generalized form)\n"
        "\t- arbitrary-precision MPFR matrices (one-sided Jacobi SVD\n"
        "\t  at the input precision, real and complex)\n"
        "\t- free-symbolic matrices (eigendecomposition of m^H . m;\n"
        "\t  the call is left unevaluated when no closed form exists)\n"
        "\n"
        "Options:\n"
        "\tTolerance -> t       :  zero out singular values below t\n"
        "\tTargetStructure ->\n"
        "\t  \"Dense\"            :  u, sigma, v all dense (default)\n"
        "\t  \"Structured\"       :  sigma returned as DiagonalMatrix[{..}]\n"
        "\n"
        "A non-rank-2 or empty matrix emits\n"
        "SingularValueDecomposition::matrix and the call is left\n"
        "unevaluated.  Generalized SVD with mismatched column counts\n"
        "emits ::matdims.  An out-of-range k or UpTo[k] emits ::sval.\n"
        "Unknown option keys / values emit ::opts and leave the call\n"
        "unevaluated.");
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
    symtab_set_docstring("FullForm",
        "FullForm[expr]\n"
        "\tprints expr as its raw internal tree (heads written before arguments\n"
        "\tin functional form, no operator or infix sugar).\n"
        "FullForm is a wrapper recognised by Print/Out; when an input evaluates\n"
        "to FullForm[expr] the wrapper is consumed by the printer and does not\n"
        "appear in the output.");
    symtab_set_docstring("Head",
        "Head[expr]\n"
        "\tgives the head of expr.\n"
        "Head[expr, h]\n"
        "\twraps the result with h, i.e. returns h[Head[expr]].\n"
        "\n"
        "For atoms, Head returns Integer, Real, BigInt, Rational, Complex,\n"
        "Symbol, or String; for a compound expression f[...], Head returns f.");
    symtab_set_docstring("Length",
        "Length[expr]\n"
        "\tgives the number of top-level elements in expr (the arity of its\n"
        "\thead).  Length of any atom is 0.");
    symtab_set_docstring("Dimensions",
        "Dimensions[expr]\n"
        "\tgives a list of the dimensions of expr.\n"
        "Dimensions[expr, n]\n"
        "\tgives the dimensions of expr down to at most level n.\n"
        "\n"
        "expr is treated as a full array only at levels where every sub-piece\n"
        "shares the same head and length; the walk halts at the first ragged\n"
        "level. Dimensions always returns a List, including the empty List {}\n"
        "for atomic expressions.");
    symtab_set_docstring("First", "First[expr] gives the first element of expr.");
    symtab_set_docstring("Last", "Last[expr] gives the last element of expr.");
    symtab_set_docstring("Most", "Most[expr] gives all but the last element of expr.");
    symtab_set_docstring("Rest", "Rest[expr] gives all but the first element of expr.");
    symtab_set_docstring("Append", "Append[expr, elem] adds elem to the end of expr.");
    symtab_set_docstring("Prepend", "Prepend[expr, elem] adds elem to the beginning of expr.");
    symtab_set_docstring("Insert", "Insert[expr, elem, n] inserts elem at position n in expr.");
    symtab_set_docstring("Delete", "Delete[expr, n] deletes the element at position n in expr.");
    symtab_set_docstring("Reverse", "Reverse[expr] reverses the order of elements in expr.");
    symtab_set_docstring("Rescale",
        "Rescale[x, {min, max}]\n"
        "\tgives x rescaled to run from 0 to 1 over the range min to max, "
        "equivalent to (x - min)/(max - min).\n"
        "Rescale[x, {min, max}, {ymin, ymax}]\n"
        "\tgives x rescaled to run from ymin to ymax over the range min to max.\n"
        "Rescale[list]\n"
        "\trescales each element of list to run from 0 to 1 over the range "
        "Min[list] to Max[list].\n"
        "Rescale threads over a list first argument and works with exact, "
        "real, complex, and symbolic quantities.");
    symtab_set_docstring("PadLeft",
        "PadLeft[list, n]\n"
        "\tmakes a list of length n by padding list with zeros on the left.\n"
        "PadLeft[list, n, x]\n"
        "\tpads by repeating the element x.\n"
        "PadLeft[list, n, {x1, x2, ...}]\n"
        "\tpads by cyclically repeating the elements xi.\n"
        "PadLeft[list, n, padding, m]\n"
        "\tleaves a margin of m elements of padding on the right.\n"
        "PadLeft[list, {n1, n2, ...}]\n"
        "\tmakes a nested list with length ni at level i.\n"
        "PadLeft[list]\n"
        "\tpads a ragged array list with zeros to make it full.\n"
        "A negative length pads on the right; a negative margin truncates "
        "trailing elements. The head of list need not be List.");
    symtab_set_docstring("PadRight",
        "PadRight[list, n]\n"
        "\tmakes a list of length n by padding list with zeros on the right.\n"
        "PadRight[list, n, x]\n"
        "\tpads by repeating the element x.\n"
        "PadRight[list, n, {x1, x2, ...}]\n"
        "\tpads by cyclically repeating the elements xi.\n"
        "PadRight[list, n, padding, m]\n"
        "\tleaves a margin of m elements of padding on the left.\n"
        "PadRight[list, {n1, n2, ...}]\n"
        "\tmakes a nested list with length ni at level i.\n"
        "PadRight[list]\n"
        "\tpads a ragged array list with zeros to make it full.\n"
        "A negative length pads on the left; a negative margin truncates "
        "leading elements. The head of list need not be List.");
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
    symtab_set_docstring("Flatten",
        "Flatten[list]\n"
        "\tflattens out nested lists, collapsing every level into a flat list\n"
        "\twith the same head as the top level.\n"
        "Flatten[list, n]\n"
        "\tflattens only the top n levels.\n"
        "Flatten[list, n, h]\n"
        "\tflattens only sublists whose head matches h, leaving other heads\n"
        "\tin place.");
    symtab_set_docstring("Partition",
        "Partition[list, n]\n"
        "\tpartitions list into non-overlapping sublists of length n; trailing\n"
        "\telements that do not fill a block are discarded.\n"
        "Partition[list, n, d]\n"
        "\tuses offset d between successive sublists; d = 1 gives a moving\n"
        "\twindow, d = n gives non-overlapping blocks.");

    // List Operations
    symtab_set_docstring("Table", "Table[expr, n]\n\tgenerates a list of n copies of expr.\nTable[expr, {i, imax}]\n\tgenerates a list of the values of expr with i running from 1 to imax.");
    symtab_set_docstring("Range", "Range[n]\n\tgenerates the list {1, 2, 3, ..., n}.\nRange[n, m]\n\tgenerates the list {n, n + 1, ..., m - 1, m}.\nRange[n, m, d]\n\tuses step d.");
    symtab_set_docstring("Array",
        "Array[f, n]\n"
        "\tgenerates a list {f[1], f[2], ..., f[n]}.\n"
        "Array[f, n, r]\n"
        "\tgenerates a list of length n starting from index r.\n"
        "Array[f, {n1, n2, ...}]\n"
        "\tgenerates an n1 x n2 x ... nested-list array with elements\n"
        "\tf[i1, i2, ...].");
    symtab_set_docstring("Union",
        "Union[list]\n"
        "\tgives the sorted list of distinct elements in list.\n"
        "Union[l1, l2, ...]\n"
        "\tgives the sorted list of distinct elements appearing in any of the\n"
        "\tinput lists (set union).\n"
        "Comparison is by canonical structural equality.");
    symtab_set_docstring("Tally", "Tally[list] counts the number of occurrences of each distinct element in list.");
    symtab_set_docstring("Commonest", "Commonest[list] gives a list of the elements that are the most common in list.\nCommonest[list, n] gives a list of the n most common elements in list.");
    symtab_set_docstring("DeleteDuplicates",
        "DeleteDuplicates[list]\n"
        "\treturns list with duplicate elements removed, keeping the first\n"
        "\toccurrence of each element and preserving the original order.\n"
        "DeleteDuplicates[list, test]\n"
        "\ttreats two elements as duplicates when test[a, b] yields True.");
    symtab_set_docstring("Split",
        "Split[list]\n"
        "\tsplits list into runs of consecutive identical elements, returning\n"
        "\ta list of these runs.\n"
        "Split[list, test]\n"
        "\tgroups runs of consecutive elements ei, ej for which test[ei, ej]\n"
        "\tyields True.");

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
    symtab_set_docstring("Map",
        "f /@ expr or Map[f, expr]\n"
        "\tapplies f to each element at level 1 of expr, preserving expr's head.\n"
        "Map[f, expr, levelspec]\n"
        "\tapplies f at the parts of expr selected by levelspec (e.g. {2} for\n"
        "\tlevel 2 only, Infinity for every level).");
    symtab_set_docstring("Apply",
        "f @@ expr or Apply[f, expr]\n"
        "\treplaces the head of expr with f.\n"
        "Apply[f, expr, levelspec]\n"
        "\tperforms the head replacement at the parts of expr specified by\n"
        "\tlevelspec; the default levelspec is {0} (top level only).");
    symtab_set_docstring("MapAll",
        "f //@ expr or MapAll[f, expr]\n"
        "\tapplies f to every subexpression in expr (equivalent to\n"
        "\tMap[f, expr, {0, Infinity}]).  Atomic leaves are wrapped too.");
    symtab_set_docstring("Through",
        "Through[p[f1, f2, ...][x1, x2, ...]]\n"
        "\tdistributes the trailing argument list across the inner functions,\n"
        "\tgiving p[f1[x1, x2, ...], f2[x1, x2, ...], ...].\n"
        "Through[expr, h]\n"
        "\tdistributes only when the outer head equals h.");
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
        "\t{m, n, s} elements m through n in steps of s");
    symtab_set_docstring("Select",
        "Select[list, crit]\n"
        "\tselects elements e of list for which crit[e] yields True, preserving\n"
        "\tthe head of list.\n"
        "Select[list, crit, n]\n"
        "\tstops after the first n matching elements.\n"
        "Select[crit]\n"
        "\tis the operator form: Select[crit][list] == Select[list, crit].");
    symtab_set_docstring("FreeQ",
        "FreeQ[expr, form]\n"
        "\tyields True if no subexpression of expr matches form, False otherwise.\n"
        "FreeQ[expr, form, levelspec]\n"
        "\trestricts the search to parts of expr at the levels specified by\n"
        "\tlevelspec.\n"
        "FreeQ[form]\n"
        "\tis the operator form: FreeQ[form][expr] == FreeQ[expr, form].");
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
    symtab_set_docstring("AtomQ",
        "AtomQ[expr]\n"
        "\tgives True if expr is an atomic object (Integer, Real, BigInt,\n"
        "\tRational, Complex, Symbol, or String), and False if expr is a\n"
        "\tcompound expression of the form head[...].");
    symtab_set_docstring("Identity", "Identity[expr] gives expr unchanged (the identity function).");
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
    symtab_set_docstring("Differences",
        "Differences[list]\n"
        "\tgives the successive differences of the elements of list.\n"
        "Differences[list, n] gives the n-th differences (length l - n).\n"
        "Differences[list, n, s] takes differences of elements step s apart\n"
        "(length l - n |s|).\n"
        "Differences[list, {n1, n2, ...}] gives the successive n_k-th differences\n"
        "at level k of a nested list; for a matrix m, Differences[m, n] (= "
        "Differences[m, {n, 0}]) differences successive rows.\n"
        "FoldList[Plus, x, Differences[list]] inverts Differences.\n"
        "Differences has the attribute Protected.");
    symtab_set_docstring("Ratios",
        "Ratios[list]\n"
        "\tgives the successive ratios list[[k+1]]/list[[k]] of the elements\n"
        "\tof list (length l - 1).\n"
        "Ratios[list, n] gives the n-th iterated ratios (length l - n); n must\n"
        "be a non-negative integer (n = 0 returns list unchanged).\n"
        "Ratios[list, {n1, n2, ...}] gives the successive n_k-th ratios at\n"
        "level k of a nested list; for a matrix m, Ratios[m, n] (= "
        "Ratios[m, {n, 0}]) takes ratios of successive rows.\n"
        "FoldList[Times, x, Ratios[list]] inverts Ratios.\n"
        "Ratios has the attribute Protected.");
    symtab_set_docstring("NumberQ",
        "NumberQ[expr]\n"
        "\tgives True if expr is an explicit number (Integer, BigInt, Rational,\n"
        "\tReal, MPFR, or Complex), and False otherwise.  Symbolic constants\n"
        "\tsuch as Pi give False; use NumericQ for those.");
    symtab_set_docstring("MachineNumberQ",
        "MachineNumberQ[expr] gives True if expr is a machine-precision real or complex number, and False otherwise.");
    symtab_set_docstring("NumericQ", "NumericQ[expr] gives True if expr is a numeric quantity, and False otherwise.\nAn expression is considered a numeric quantity if it is either an explicit number or a mathematical constant such as Pi, or is a function that has attribute NumericFunction and all of whose arguments are numeric quantities.");
    symtab_set_docstring("Positive",
        "Positive[x]\n"
        "\tgives True if x is a positive real number, and False if x is a\n"
        "manifestly negative real number, a non-real complex number, or zero.\n"
        "For non-numeric x the expression is left unevaluated. Positive is\n"
        "Listable, so it threads over lists element by element.");
    symtab_set_docstring("Negative",
        "Negative[x]\n"
        "\tgives True if x is a negative real number, and False if x is a\n"
        "manifestly non-negative real number (including zero) or a non-real\n"
        "complex number. For non-numeric x the expression is left unevaluated.\n"
        "Negative is Listable, so it threads over lists element by element.");
    symtab_set_docstring("NonNegative",
        "NonNegative[x]\n"
        "\tgives True if x is a real number that is positive or zero, and False\n"
        "if x is a manifestly negative real number or a non-real complex number.\n"
        "For non-numeric x the expression is left unevaluated. NonNegative is\n"
        "Listable, so it threads over lists element by element.");
    symtab_set_docstring("NonPositive",
        "NonPositive[x]\n"
        "\tgives True if x is a real number that is negative or zero, and False\n"
        "if x is a manifestly positive real number or a non-real complex number.\n"
        "For non-numeric x the expression is left unevaluated. NonPositive is\n"
        "Listable, so it threads over lists element by element.");
    symtab_set_docstring("IntegerQ",
        "IntegerQ[expr]\n"
        "\tgives True if expr is an Integer or BigInt, False otherwise.\n"
        "Returns False on rationals with denominator > 1, reals, and symbolic\n"
        "expressions (even those that are integer-valued, e.g. 2 Pi / Pi).");
    symtab_set_docstring("EvenQ", "EvenQ[n] gives True if n is an even integer (Integer or BigInt), False otherwise.");
    symtab_set_docstring("OddQ", "OddQ[n] gives True if n is an odd integer (Integer or BigInt), False otherwise.");
    symtab_set_docstring("PrimeQ",
        "PrimeQ[n]\n"
        "\tgives True if n is a prime integer, False otherwise.\n"
        "PrimeQ[z]\n"
        "\tfor a Gaussian integer z = a + b I, gives True if z is a Gaussian prime.\n"
        "PrimeQ[n, GaussianIntegers -> True]\n"
        "\ttests primality of n in Z[i] rather than in Z.\n"
        "Primality is tested with GMP's mpz_probab_prime_p using 25 Miller-Rabin "
        "rounds on top of a Baillie-PSW pre-screen, so composite false positives "
        "have probability below 4^-25 (definite for n < 2^64).");
    symtab_set_docstring("SquareFreeQ",
        "SquareFreeQ[expr]\n"
        "\tgives True if expr is a square-free polynomial or number, and False otherwise.\n"
        "SquareFreeQ[expr, vars]\n"
        "\tgives True if expr is square-free with respect to the variables vars.\n"
        "Option GaussianIntegers -> True | False | Automatic switches to Gaussian integers.");
    symtab_set_docstring("IrreduciblePolynomialQ",
        "IrreduciblePolynomialQ[poly]\n"
        "\tgives True if poly is an irreducible polynomial over the rationals.\n"
        "Option Extension -> alpha | {alpha_i} tests irreducibility over the field "
        "extension generated by the algebraic numbers alpha_i.\n"
        "Option Extension -> Automatic extends Q by every algebraic-number coefficient "
        "in poly; Extension -> All tests absolute irreducibility over the complex numbers.\n"
        "Option GaussianIntegers -> True tests irreducibility over the Gaussian rationals.");
    symtab_set_docstring("PossibleZeroQ",
        "PossibleZeroQ[expr] gives True if symbolic and numerical methods suggest "
        "that expr has value zero, and False otherwise.\n"
        "The general problem of deciding whether an expression is zero is "
        "undecidable; PossibleZeroQ is a quick but not always accurate test.");
    symtab_set_docstring("PolynomialQ",
        "PolynomialQ[expr, var]\n"
        "\tgives True if expr is a polynomial in var, False otherwise.\n"
        "PolynomialQ[expr, {v1, v2, ...}]\n"
        "\tgives True if expr is a polynomial in all of the vi simultaneously.\n"
        "Checks that expr expands to a sum of products of non-negative integer\n"
        "powers of the vars with var-free coefficients.");
    symtab_set_docstring("ListQ", "ListQ[expr] gives True if expr is a list (head List), False otherwise.");
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
    symtab_set_docstring("MatchQ",
        "MatchQ[expr, form]\n"
        "\tgives True if expr matches the pattern form, False otherwise.\n"
        "MatchQ[form]\n"
        "\tis the operator form: MatchQ[form][expr] == MatchQ[expr, form].\n"
        "Pattern matching honours sequence variables (__, ___), PatternTest,\n"
        "Condition, attribute-driven flattening / ordering, and the surrounding\n"
        "$Assumptions / DownValues environment.");
    symtab_set_docstring("HermitianMatrixQ",
        "HermitianMatrixQ[m]\n"
        "\tgives True if m is explicitly Hermitian (m == ConjugateTranspose[m]),\n"
        "\tand False otherwise.\n"
        "\n"
        "Options:\n"
        "\tSameTest  -> Automatic   function used to test equality of entries.\n"
        "\tTolerance -> Automatic   numeric tolerance for approximate matrices.\n"
        "\n"
        "With SameTest -> f, entries m[i,j] and Conjugate[m[j,i]] are taken to be\n"
        "equal when f[m[i,j], Conjugate[m[j,i]]] gives True.  With Tolerance -> t,\n"
        "entries are accepted when Abs[m[i,j] - Conjugate[m[j,i]]] <= t.\n"
        "Diagonal entries must satisfy the same test (i.e. be purely real for\n"
        "numeric matrices).");
    symtab_set_docstring("SymmetricMatrixQ",
        "SymmetricMatrixQ[m]\n"
        "\tgives True if m is explicitly symmetric (m == Transpose[m]),\n"
        "\tand False otherwise.\n"
        "\n"
        "Options:\n"
        "\tSameTest  -> Automatic   function used to test equality of entries.\n"
        "\tTolerance -> Automatic   numeric tolerance for approximate matrices.\n"
        "\n"
        "With SameTest -> f, entries m[i,j] and m[j,i] are taken to be equal\n"
        "when f[m[i,j], m[j,i]] gives True.  With Tolerance -> t, entries are\n"
        "accepted when Abs[m[i,j] - m[j,i]] <= t.  SymmetricMatrixQ uses the\n"
        "definition m^T == m for both real- and complex-valued matrices, so a\n"
        "complex symmetric matrix need not be Hermitian.");
    symtab_set_docstring("SquareMatrixQ",
        "SquareMatrixQ[m]\n"
        "\tgives True if m is a square matrix (Dimensions[m] == {n, n}),\n"
        "\tand False otherwise.\n"
        "\n"
        "Works for symbolic as well as numerical matrices.  Returns False on\n"
        "non-list, ragged, rectangular, empty, or higher-rank tensor inputs.");
    symtab_set_docstring("DiagonalMatrixQ",
        "DiagonalMatrixQ[m]\n"
        "\tgives True if m is diagonal, and False otherwise.\n"
        "DiagonalMatrixQ[m, k]\n"
        "\tgives True if m has nonzero elements only on the k-th diagonal,\n"
        "\tand False otherwise.  Positive k refers to superdiagonals above\n"
        "\tthe main diagonal; negative k refers to subdiagonals below it.\n"
        "\tWorks for rectangular as well as square matrices.\n"
        "\n"
        "Option:\n"
        "\tTolerance -> Automatic   numeric tolerance for approximate matrices.\n"
        "\n"
        "With Tolerance -> t, off-diagonal entries e are taken to be zero\n"
        "when Abs[e] <= t evaluates to True.  Without a tolerance the test\n"
        "is structural: only literal numeric zeros (Integer 0, Real 0.0,\n"
        "BigInt 0) count as zero.  Returns False on non-matrix, ragged,\n"
        "empty (i.e. {}), or higher-rank tensor inputs; an n-by-0 matrix\n"
        "(e.g. {{}, {}}) is vacuously diagonal.");
    symtab_set_docstring("UpperTriangularMatrixQ",
        "UpperTriangularMatrixQ[m]\n"
        "\tgives True if m is upper triangular, and False otherwise.\n"
        "UpperTriangularMatrixQ[m, k]\n"
        "\tgives True if m is upper triangular starting from the k-th\n"
        "\tdiagonal (every entry m[i,j] with j - i < k is zero), and\n"
        "\tFalse otherwise.  Positive k refers to superdiagonals above\n"
        "\tthe main diagonal; negative k refers to subdiagonals below it.\n"
        "\tWorks for rectangular as well as square matrices.\n"
        "\n"
        "Option:\n"
        "\tTolerance -> Automatic   numeric tolerance for approximate matrices.\n"
        "\n"
        "With Tolerance -> t, sub-diagonal entries e are taken to be zero\n"
        "when Abs[e] <= t evaluates to True.  Without a tolerance the test\n"
        "is structural: only literal numeric zeros (Integer 0, Real 0.0,\n"
        "BigInt 0) count as zero.  Returns False on non-matrix, ragged,\n"
        "empty (i.e. {}), or higher-rank tensor inputs; an n-by-0 matrix\n"
        "(e.g. {{}, {}}) is vacuously upper triangular.");
    symtab_set_docstring("PositiveDefiniteMatrixQ",
        "PositiveDefiniteMatrixQ[m]\n"
        "\tgives True if m is explicitly positive definite, and False\n"
        "\totherwise.\n"
        "\n"
        "A matrix m is positive definite if Re[Conjugate[x] . m . x] > 0\n"
        "for every nonzero vector x.  Equivalently, the Hermitian part\n"
        "(m + ConjugateTranspose[m]) / 2 has only positive eigenvalues.\n"
        "The test is performed by attempting a Cholesky factorisation of\n"
        "the Hermitian part; on numeric matrices this is dispatched to\n"
        "BLAS/LAPACK's dpotrf (real) or zpotrf (complex) when available.\n"
        "Returns False on non-numeric, non-square, ragged, empty, or\n"
        "higher-rank tensor inputs.");
    symtab_set_docstring("NegativeDefiniteMatrixQ",
        "NegativeDefiniteMatrixQ[m]\n"
        "\tgives True if m is explicitly negative definite, and False\n"
        "\totherwise.\n"
        "\n"
        "A matrix m is negative definite if Re[Conjugate[x] . m . x] < 0\n"
        "for every nonzero vector x.  Equivalently, -m is positive\n"
        "definite, i.e. the negated Hermitian part has only positive\n"
        "eigenvalues.  The test is performed by attempting a Cholesky\n"
        "factorisation of -(m + ConjugateTranspose[m]) / 2; on numeric\n"
        "matrices this is dispatched to BLAS/LAPACK's dpotrf (real) or\n"
        "zpotrf (complex) when available.  Returns False on non-numeric,\n"
        "non-square, ragged, empty, or higher-rank tensor inputs.");

    // Trigonometric (all Listable, NumericFunction; exact reductions at
    // rational multiples of Pi via the standard table; numeric inputs route
    // to libm / MPFR; symbolic arguments stay unevaluated).
    symtab_set_docstring("Sin",
        "Sin[z]\n"
        "\tgives the sine of z (argument in radians).\n"
        "Sin is Listable. Numeric inputs are evaluated via libm (Real) or MPFR\n"
        "(arbitrary precision); rational multiples of Pi reduce to exact values.");
    symtab_set_docstring("Cos",
        "Cos[z]\n"
        "\tgives the cosine of z (argument in radians).\n"
        "Cos is Listable. Numeric inputs route to libm / MPFR; rational\n"
        "multiples of Pi reduce to exact values.");
    symtab_set_docstring("Tan",
        "Tan[z]\n"
        "\tgives the tangent of z. Equivalent to Sin[z] / Cos[z].\n"
        "Tan is Listable. Singularities at z = Pi/2 + k Pi yield ComplexInfinity.");
    symtab_set_docstring("Cot",
        "Cot[z]\n"
        "\tgives the cotangent of z. Equivalent to Cos[z] / Sin[z].\n"
        "Cot is Listable. Singularities at z = k Pi yield ComplexInfinity.");
    symtab_set_docstring("Sec",
        "Sec[z]\n"
        "\tgives the secant of z (= 1 / Cos[z]).\n"
        "Sec is Listable. Singularities at z = Pi/2 + k Pi yield ComplexInfinity.");
    symtab_set_docstring("Csc",
        "Csc[z]\n"
        "\tgives the cosecant of z (= 1 / Sin[z]).\n"
        "Csc is Listable. Singularities at z = k Pi yield ComplexInfinity.");
    symtab_set_docstring("ArcSin",
        "ArcSin[z]\n"
        "\tgives the principal inverse sine of z, in [-Pi/2, Pi/2] for real z\n"
        "\tin [-1, 1].\n"
        "ArcSin is Listable. Branch cuts run along the real axis with |z| > 1.");
    symtab_set_docstring("ArcCos",
        "ArcCos[z]\n"
        "\tgives the principal inverse cosine of z, in [0, Pi] for real z\n"
        "\tin [-1, 1].\n"
        "ArcCos is Listable. Branch cuts run along the real axis with |z| > 1.");
    symtab_set_docstring("ArcTan",
        "ArcTan[z]\n"
        "\tgives the principal inverse tangent of z, in (-Pi/2, Pi/2).\n"
        "ArcTan[y, x]\n"
        "\tgives the argument of the complex number x + I y, in (-Pi, Pi]\n"
        "\t(two-argument atan2 form).\n"
        "ArcTan is Listable.");

    // Hyperbolic (Listable, NumericFunction; numeric inputs route to libm /
    // MPFR; integer / Pi-multiple arguments reduce when applicable).
    symtab_set_docstring("Sinh",
        "Sinh[z]\n"
        "\tgives the hyperbolic sine of z, (Exp[z] - Exp[-z]) / 2.\n"
        "Sinh is Listable.");
    symtab_set_docstring("Cosh",
        "Cosh[z]\n"
        "\tgives the hyperbolic cosine of z, (Exp[z] + Exp[-z]) / 2.\n"
        "Cosh is Listable.");
    symtab_set_docstring("Tanh",
        "Tanh[z]\n"
        "\tgives the hyperbolic tangent of z, Sinh[z] / Cosh[z].\n"
        "Tanh is Listable.");

    // Log/Exp
    symtab_set_docstring("Log",
        "Log[z]\n"
        "\tgives the principal natural logarithm of z, with branch cut along\n"
        "\tthe negative real axis.\n"
        "Log[b, z]\n"
        "\tgives the logarithm to base b, i.e. Log[z] / Log[b].\n"
        "Log is Listable. Log[1] = 0, Log[E] = 1, Log[E^n] = n for symbolic n.\n"
        "Numeric inputs route to libm / MPFR; negative reals yield I Pi + Log[|z|].");
    symtab_set_docstring("Exp",
        "Exp[z]\n"
        "\tgives the exponential E^z.\n"
        "Exp is Listable. Exp[0] = 1, Exp[Log[x]] -> x, Exp[I Pi] = -1.\n"
        "Numeric inputs route to libm / MPFR.");

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
        "\trewrites products and powers of trigonometric functions in expr in\n"
        "\tterms of trigonometric functions with combined arguments.\n"
        "TrigReduce operates on both circular and hyperbolic functions; given a\n"
        "trigonometric polynomial it typically yields a linear expression\n"
        "involving trigonometric functions with more complicated arguments\n"
        "(broadly the inverse of TrigExpand).\n"
        "TrigReduce automatically threads over lists, equations, inequalities,\n"
        "and logic functions.");

    // Piecewise / Rounding
    symtab_set_docstring("Floor",
        "Floor[x]\n"
        "\tgives the greatest integer less than or equal to x.\n"
        "Floor[x, a]\n"
        "\tgives the greatest multiple of a less than or equal to x.\n"
        "Floor is Listable. Exact (Integer / BigInt / Rational) inputs return\n"
        "exact integers; Real / MPFR inputs are rounded toward -Infinity at\n"
        "the input precision; symbolic inputs stay unevaluated.");
    symtab_set_docstring("Ceiling",
        "Ceiling[x]\n"
        "\tgives the smallest integer greater than or equal to x.\n"
        "Ceiling[x, a]\n"
        "\tgives the smallest multiple of a greater than or equal to x.\n"
        "Ceiling is Listable. Exact inputs return exact integers; Real / MPFR\n"
        "inputs are rounded toward +Infinity at the input precision.");
    symtab_set_docstring("ContinuedFraction",
        "ContinuedFraction[x, n]\n"
        "\tgives a list of the first n terms in the continued-fraction\n"
        "\trepresentation of x.\n"
        "ContinuedFraction[x]\n"
        "\tgives all terms determinable from the precision of x.\n"
        "The list {a1, a2, a3, ...} corresponds to a1 + 1/(a2 + 1/(a3 + ...)).\n"
        "Exact rationals give a finite (canonical, last term >= 2) expansion.\n"
        "For Sqrt[d] with d a non-square integer the no-count form returns\n"
        "{a1, ..., {b1, ...}}, the bracketed block repeating cyclically. Inexact\n"
        "Real / MPFR inputs yield terms only as far as the precision determines\n"
        "them. ContinuedFraction is Listable.");
    symtab_set_docstring("FromContinuedFraction",
        "FromContinuedFraction[{a1, a2, ..., an}]\n"
        "\treconstructs a1 + 1/(a2 + 1/(a3 + ... + 1/an)). The terms may be\n"
        "\tsymbolic; the result is the convergent in nested (un-expanded) form.\n"
        "FromContinuedFraction[{a1, ..., am, {b1, ..., bk}}]\n"
        "\tgives the exact quadratic irrational whose continued-fraction terms\n"
        "\tbegin with the ai then cycle through the bi forever; all ai and bi\n"
        "\tmust be integers. FromContinuedFraction[{}] is 0. It is the inverse\n"
        "\tof ContinuedFraction.");
    symtab_set_docstring("Round",
        "Round[x]\n"
        "\trounds x to the nearest integer, breaking ties to the nearest even\n"
        "\tinteger (banker's rounding).\n"
        "Round[x, a]\n"
        "\trounds x to the nearest multiple of a.\n"
        "Round is Listable. Exact inputs return exact integers; Real / MPFR\n"
        "inputs round at the input precision.");
    symtab_set_docstring("UnitStep",
        "UnitStep[x]\n"
        "\tgives 0 for x < 0 and 1 for x >= 0 (the value at 0 is 1).\n"
        "UnitStep[x1, x2, ...]\n"
        "\tgives 1 only when none of the xi are negative, otherwise 0.\n"
        "UnitStep[] is 1. The result is always exact. Exact symbolic real\n"
        "arguments are resolved by numerical certification; non-real or\n"
        "unresolved arguments are left unevaluated. UnitStep is Listable and\n"
        "Orderless.");
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
    symtab_set_docstring("D", "D[f, x] gives the partial derivative of f with respect to x.\nD[f, {x, n}] gives the nth partial derivative.\nD[f, x, y, ...] gives the mixed derivative.\nD[f, x, NonConstants -> {y, ...}] treats the listed symbols as implicit functions of x.\nDistributes over Equal: D[a == b, x] gives D[a, x] == D[b, x].");
    symtab_set_docstring("Dt", "Dt[f] gives the total derivative of f.\nDt[f, x] gives the total derivative of f with respect to x.\nDt[f, {x, n}] gives the nth total derivative.");
    symtab_set_docstring("Derivative",
        "f' represents the derivative of a function f of one argument.\n"
        "Derivative[n1, n2, ...][f] is the general form, representing a function\n"
        "obtained from f by differentiating n1 times with respect to the first\n"
        "argument, n2 times with respect to the second argument, and so on.\n"
        "\n"
        "f' is equivalent to Derivative[1][f]; f'' evaluates to Derivative[2][f].\n"
        "Derivative is a functional operator acting on functions to give derivative\n"
        "functions. Derivative is generated when D is applied to functions whose\n"
        "derivatives the system does not know.\n"
        "\n"
        "Mathilda attempts to convert Derivative[n1,...,nm][f] to a pure function.\n"
        "When f is a symbol carrying DownValues, the evaluator rewrites the head\n"
        "as Function[{t1,...,tm}, f[t1,...,tm]] with the rule expanded into the\n"
        "body, then differentiates that pure function. If no DownValue matches,\n"
        "the original Derivative form is returned.\n"
        "\n"
        "Attributes: Protected, ReadProtected.");

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
    symtab_set_docstring("If",
        "If[cond, t]\n"
        "\tgives t if cond evaluates to True; gives Null otherwise.\n"
        "If[cond, t, f]\n"
        "\tgives t if cond is True, f if False, and is left unevaluated\n"
        "\tif cond is neither.\n"
        "If[cond, t, f, u]\n"
        "\talso supplies u as the result when cond is neither True nor False.\n"
        "If has attribute HoldRest: only the branch chosen by cond is evaluated.");
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
    symtab_set_docstring("InterpolatingFunction",
        "InterpolatingFunction[domain, table]\n"
        "\trepresents an approximate function whose values are found by\n"
        "\tinterpolation. domain is {{x1min, x1max}, ...} with one interval\n"
        "\tper dimension; table is a list of {coord, value} data points on a\n"
        "\tfull tensor grid (coord is a scalar for 1-D, an {x1, ..., xm} list\n"
        "\tfor m-D).\n"
        "InterpolatingFunction[...][x1, ..., xm]\n"
        "\tgives the interpolated value using tensor-product piecewise-\n"
        "\tpolynomial (default order 3) interpolation. Arguments outside the\n"
        "\tdomain are extrapolated with a warning.\n"
        "Derivative[d1, ..., dm][InterpolatingFunction[...]]\n"
        "\tgives an InterpolatingFunction for the mixed partial derivative.\n"
        "In standard output only the domain is shown; the rest is <>.");
    symtab_set_docstring("Interpolation",
        "Interpolation[data]\n"
        "\tconstructs an InterpolatingFunction that interpolates data, given\n"
        "\tas {f1, f2, ...} (values at x = 1, 2, ...), {{x1, f1}, ...} (values\n"
        "\tat given abscissae), or {{{x1, y1, ...}, f1}, ...} (an m-D tensor\n"
        "\tgrid).\n"
        "Interpolation[data, x]\n"
        "\tbuilds the interpolating function and evaluates it at x (a number,\n"
        "\tor a coordinate list in m-D).\n"
        "Interpolation[{{{x1,...}, f1, df1, ddf1, ...}, ...}]\n"
        "\treproduces supplied derivatives at the nodes (df = gradient, ddf =\n"
        "\tHessian, ...) by tensor-product Hermite interpolation.\n"
        "Interpolation[data, InterpolationOrder -> n]\n"
        "\tuses piecewise-polynomial pieces of degree n (default 3; 0 gives a\n"
        "\tpiecewise-constant and 1 a piecewise-linear interpolant).\n"
        "Interpolation[data, Method -> m]\n"
        "\tselects \"Spline\" (natural/cyclic cubic spline) or \"Hermite\"\n"
        "\t(piecewise cubic Hermite with estimated slopes).\n"
        "Interpolation[data, PeriodicInterpolation -> True]\n"
        "\tbuilds a periodic interpolant (period = the data span; the data must\n"
        "\trepeat its first sample at the last). A per-dimension {True, False}\n"
        "\tlist selects periodicity per axis.\n"
        "Vector- or array-valued samples (f_i a list) are interpolated\n"
        "component-wise and return an array of the same shape.\n"
        "Works at machine or arbitrary (MPFR) precision, matching the data.");
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
    symtab_set_docstring("ConditionalExpression",
        "ConditionalExpression[expr, cond]\n"
        "\tis a symbolic construct that represents expr when cond is True.\n"
        "ConditionalExpression[expr, True] evaluates to expr.\n"
        "ConditionalExpression[expr, False] evaluates to Undefined.\n"
        "Nested forms collapse: ConditionalExpression[ConditionalExpression[e, c1], c2] evaluates to ConditionalExpression[e, c1 && c2].\n"
        "ConditionalExpression has attribute Protected.");
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
    symtab_set_docstring("Set",
        "lhs = rhs or Set[lhs, rhs]\n"
        "\tevaluates rhs once and assigns the result to lhs.  When lhs is a\n"
        "\tsymbol, the assignment is stored as an OwnValue; when lhs has the\n"
        "\tform f[args...] it is stored as a DownValue on f.  Set has attribute\n"
        "\tHoldFirst so lhs is not evaluated before the assignment.");

    // Options machinery
    symtab_set_docstring("Options",
        "Options[f] gives the list of default option rules for the symbol f.\n"
        "\tOptions[expr] gives the options explicitly set in an expression such\n"
        "\tas a graphics object.  Options[obj, name] gives the setting for the\n"
        "\tnamed option; Options[obj, {names}] gives a list of settings.  Assign\n"
        "\tto Options[f] to redefine all default options at once.");
    symtab_set_docstring("SetOptions",
        "SetOptions[s, name -> value, ...] sets default options for the symbol\n"
        "\ts and returns the new Options[s].  It can change Protected (but not\n"
        "\tLocked) symbols, and only changes existing options -- an unknown name\n"
        "\traises SetOptions::optnf.  Use AppendTo[Options[s], ...] to add one.");
    symtab_set_docstring("OptionValue",
        "OptionValue[name] gives the value of an option named name in the\n"
        "\toptions matched by OptionsPattern[] in the enclosing rule.\n"
        "\tOptionValue[f, name] uses options associated with the head f;\n"
        "\tOptionValue[f, opts, name] resolves from the explicit rules opts then\n"
        "\tthe defaults from f; a trailing Hold wraps the result in Hold.");

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
    symtab_set_docstring("SetDelayed",
        "lhs := rhs or SetDelayed[lhs, rhs]\n"
        "\tassigns rhs to lhs as a delayed rule: rhs is held and evaluated\n"
        "\teach time the rule fires (with bindings from lhs substituted in),\n"
        "\tnot at assignment time.  SetDelayed has attribute HoldAll.");
    symtab_set_docstring("Default",
        "Default[f]\n"
        "\tgives the default value supplied for a missing optional argument of\n"
        "\tf when the pattern _. (Optional[Blank[]]) appears in a rule.\n"
        "Default[f, i]\n"
        "\tgives the default value for the i-th argument position of f.");
    symtab_set_docstring("Optional",
        "patt:def or Optional[patt, def]\n"
        "\tis a pattern object that matches patt if it is present; if patt is\n"
        "\tomitted from the argument sequence, def is used in its place.\n"
        "patt_. (sugar for Optional[patt_, Default[f]]) draws the default value\n"
        "from Default[f] at the call site.");
    symtab_set_docstring("Longest", "Longest[p] is a pattern object that matches the longest sequence consistent with the pattern p.");
    symtab_set_docstring("Shortest", "Shortest[p] is a pattern object that matches the shortest sequence consistent with the pattern p.");
    symtab_set_docstring("Repeated", "p.. or Repeated[p] is a pattern object that represents a sequence of one or more expressions, each matching p.\nRepeated[p, max] represents from 1 to max expressions matching p.\nRepeated[p, {min, max}] represents between min and max expressions matching p.\nRepeated[p, {n}] represents exactly n expressions matching p.");
    symtab_set_docstring("RepeatedNull", "p... or RepeatedNull[p] is a pattern object that represents a sequence of zero or more expressions, each matching p.\nRepeatedNull[p, max] represents from 0 to max expressions matching p.\nRepeatedNull[p, {min, max}] represents between min and max expressions matching p.\nRepeatedNull[p, {n}] represents exactly n expressions matching p.");
    symtab_set_docstring("Blank", "_ or Blank[] represents any single expression.\n_h or Blank[h] represents any single expression with head h.");
    symtab_set_docstring("BlankSequence", "__ or BlankSequence[] represents a sequence of one or more expressions.");
    symtab_set_docstring("BlankNullSequence", "___ or BlankNullSequence[] represents a sequence of zero or more expressions.");
    symtab_set_docstring("Clear",
        "Clear[s1, s2, ...]\n"
        "\tclears all OwnValues and DownValues attached to the named symbols,\n"
        "\tleaving attributes and the symbol itself intact.\n"
        "Clear has attribute HoldAll; Protected symbols are skipped with a\n"
        "diagnostic.");
    symtab_set_docstring("Unset",
        "Unset[lhs] or lhs =.\n"
        "\tremoves any rule whose left-hand side is lhs, up to renaming of\n"
        "\tpattern variables. A bare symbol clears its value; a function form\n"
        "\tclears the matching definition on the head symbol.\n"
        "Unset has attribute HoldFirst; Protected symbols are not affected.");
    symtab_set_docstring("ClearAll",
        "ClearAll[s1, s2, ...]\n"
        "\tclears all values, definitions, attributes and messages for the\n"
        "\tnamed symbols. ClearAll[{s1, s2, ...}] accepts a list of specs.\n"
        "ClearAll has attribute HoldAll; symbols with attribute Locked or\n"
        "Protected are not affected.");
    symtab_set_docstring("Remove",
        "Remove[s1, s2, ...]\n"
        "\tremoves the named symbols completely, deleting their definitions\n"
        "\tfrom the symbol table. Remove[{s1, s2, ...}] accepts a list of specs.\n"
        "Remove has attribute HoldAll; symbols with attribute Locked or\n"
        "Protected are not affected.");
    symtab_set_docstring("Protect",
        "Protect[s1, s2, ...]\n"
        "\tsets the attribute Protected for the named symbols and returns the\n"
        "\tlist of their names. Protect[{s1, s2, ...}] accepts a list of specs.\n"
        "Protect has attribute HoldAll; Locked symbols are not affected.");
    symtab_set_docstring("Unprotect",
        "Unprotect[s1, s2, ...]\n"
        "\tremoves the attribute Protected from the named symbols and returns\n"
        "\tthe list of their names. Unprotect[{s1, ...}] accepts a list of specs.\n"
        "Unprotect has attribute HoldAll; Locked symbols are not affected.");
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
    symtab_set_docstring("Rule",
        "lhs -> rhs or Rule[lhs, rhs]\n"
        "\trepresents an immediate rewrite rule: rhs is evaluated when the\n"
        "\trule object is constructed, then matched against lhs at use.");
    symtab_set_docstring("RuleDelayed",
        "lhs :> rhs or RuleDelayed[lhs, rhs]\n"
        "\trepresents a delayed rewrite rule: rhs is held and evaluated only\n"
        "\teach time the rule fires, after the pattern bindings on lhs are\n"
        "\tsubstituted into rhs.");
    symtab_set_docstring("Replace",
        "Replace[expr, rules]\n"
        "\ttries to match expr at the top level against rules and returns the\n"
        "\trewritten form; if no rule matches, returns expr unchanged.\n"
        "Replace[expr, rules, levelspec]\n"
        "\tapplies rules only at the parts of expr specified by levelspec.\n"
        "Matching tries each rule in order and uses the first that succeeds;\n"
        "rules may be a single rule or a list of rules.");
    symtab_set_docstring("ReplaceAll",
        "expr /. rules or ReplaceAll[expr, rules]\n"
        "\ttraverses expr top-down and applies the first matching rule at each\n"
        "\tsubexpression. A matched subexpression is replaced and NOT recursed\n"
        "\tinto further -- ReplaceAll is a single pass, not a fixed point.");
    symtab_set_docstring("ReplaceRepeated",
        "expr //. rules or ReplaceRepeated[expr, rules]\n"
        "\trepeatedly applies ReplaceAll[expr, rules] until the result stops\n"
        "\tchanging, then returns the fixed point.  Useful for chained rewrite\n"
        "\tsystems; subject to the same recursion-limit guard as evaluator\n"
        "\tfixed-point iteration.");
    symtab_set_docstring("Print",
        "Print[expr1, expr2, ...]\n"
        "\tprints each argument to stdout, concatenated without separator and\n"
        "\tfollowed by a newline, and returns Null.  Arguments are formatted in\n"
        "\tthe default output form (matching the REPL's Out display).");

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

    symtab_set_docstring("InputForm",
        "InputForm[expr]\n"
        "\tprints expr in a form suitable to be re-read by the parser, using\n"
        "\toperator syntax (a + b, not Plus[a, b]) and explicit string quotes.\n"
        "Like FullForm, InputForm is a printer wrapper: it is consumed during\n"
        "output and does not appear in the printed result.");
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
        "\tgives the printed form of expr (as InputForm) as a String.\n"
        "ToString[expr, form]\n"
        "\tuses the specified output form.\n"
        "Supported forms: InputForm (default), FullForm, TeXForm.");
    symtab_set_docstring("ToExpression",
        "ToExpression[input]\n"
        "\tparses input (a String) as Mathilda input and returns the resulting\n"
        "\texpression after evaluation.\n"
        "ToExpression[input, form]\n"
        "\tuses interpretation rules for the specified form. form may be\n"
        "\tInputForm or FullForm (both currently use the same parser).\n"
        "ToExpression[input, form, h]\n"
        "\twraps the head h around the parsed expression before evaluation;\n"
        "\tuse h = Hold to obtain the unevaluated parsed form.\n"
        "Returns $Failed if a syntax error is encountered.");
    symtab_set_docstring("Symbol",
        "Symbol[\"name\"]\n"
        "\trefers to the symbol with the specified name, creating it in\n"
        "\t$Context if none yet exists.\n"
        "\n"
        "All symbols, whether explicitly entered using Symbol or not, have head\n"
        "Symbol; x_Symbol matches any symbol. The name string may contain\n"
        "letters, letter-like forms, or digits but must not start with a digit.\n"
        "A backtick (`) separates context prefixes; a leading backtick makes the\n"
        "name relative to the current context $Context.\n"
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
    symtab_set_docstring("SeriesCoefficient",
        "SeriesCoefficient[f, {x, x0, k}]\n"
        "\tgives the coefficient of (x - x0)^k in the power-series expansion of f\n"
        "\tabout x = x0. Works for a concrete integer index k and a finite expansion\n"
        "point, for any f that Series can expand. HoldAll, Protected.");
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
    symtab_set_docstring("PowerExpand",
        "PowerExpand[expr]\n\texpands (a b)^c to a^c b^c and (a^b)^c to a^(b c), and expands\n"
        "\tLog and Arg of products and powers.\n"
        "PowerExpand[expr, {x1, x2, ...}]\n\texpands only with respect to the listed variables.\n"
        "The transformations are correct in general only when c is an integer or a and b\n"
        "\tare positive reals; PowerExpand otherwise disregards branch cuts.\n"
        "With the Assumptions option, Assumptions->True gives a universally-correct result\n"
        "\tand Assumptions->assum a result valid under assum.\n"
        "PowerExpand threads over lists, equations, inequalities, and logic functions.");
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
    symtab_set_docstring("Coefficient",
        "Coefficient[expr, form]\n"
        "\tgives the coefficient of form^1 in expr.  form is matched\n"
        "\tstructurally against the bases of products in the expanded form\n"
        "\tof expr.\n"
        "Coefficient[expr, form, n]\n"
        "\tgives the coefficient of form^n.  n may be a non-negative integer\n"
        "\tor (for Laurent / Puiseux expressions) a rational.");
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
    symtab_set_docstring("Collect",
        "Collect[expr, x]\n"
        "\texpands expr and gathers terms with the same power of x, returning\n"
        "\ta sum of the form Sum[c_k x^k] with each c_k free of x.\n"
        "Collect[expr, {x1, x2, ...}]\n"
        "\tcollects with respect to each xi in turn (nested grouping).\n"
        "Collect[expr, x, f]\n"
        "\tapplies f to each coefficient before re-assembling the sum, useful\n"
        "\tfor f = Simplify or f = Factor.");
    symtab_set_docstring("PolynomialMod",
        "PolynomialMod[poly, m]\n"
        "\treduces poly modulo m.  If m is an integer, each coefficient of\n"
        "\tpoly is reduced to a canonical residue in {0, ..., m-1}.  If m is a\n"
        "\tpolynomial, poly is reduced modulo m as polynomials over the\n"
        "\trationals (in contrast to PolynomialRemainder, the leading\n"
        "\tcoefficient of m is not normalised).\n"
        "PolynomialMod[poly, {m1, m2, ...}]\n"
        "\treduces modulo each mi in turn.");
    symtab_set_docstring("FactorSquareFree",
        "FactorSquareFree[poly]\n"
        "\twrites poly as a product of pairwise-coprime square-free factors,\n"
        "\tcollecting repeated factors into powers.\n"
        "Computed via the Yun / Musser square-free decomposition using\n"
        "polynomial GCDs of poly with its derivative; cheaper than full Factor\n"
        "and sufficient when only multiplicities are needed.");
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
    symtab_set_docstring("Variables",
        "Variables[poly]\n"
        "\tgives the sorted list of independent variables that appear as bases\n"
        "\tof non-numeric subexpressions in poly.\n"
        "Walks the expression tree and collects symbols and compound forms\n"
        "(e.g. Sin[x], a[i]) that occur outside numeric arithmetic; duplicates\n"
        "are removed via canonical order.");
    symtab_set_docstring("Decompose",
        "Decompose[poly, x]\n"
        "\tdecomposes the univariate polynomial poly into the deepest possible\n"
        "\tcomposition {p1, p2, ..., pk} such that poly == p1[p2[...[pk[x]]...]],\n"
        "\twith each pi a polynomial of degree >= 2 in x.\n"
        "\tReturns {poly} if no nontrivial decomposition exists.");
    symtab_set_docstring("HornerForm",
        "HornerForm[poly]\n"
        "\trewrites the univariate polynomial poly in nested (Horner) form,\n"
        "\twhich evaluates in n multiplications and n additions instead of\n"
        "\tthe naive 2n.\n"
        "HornerForm[poly, var]\n"
        "\tuses var as the recursion variable for multivariate poly.\n"
        "HornerForm[poly1 / poly2, vars1, vars2]\n"
        "\tputs a rational function in Horner form, nested with respect to\n"
        "\tvars1 in the numerator and vars2 in the denominator.");
    symtab_set_docstring("MinimalPolynomial",
        "MinimalPolynomial[s, x]\n"
        "\tgives the lowest-degree polynomial in x with integer coefficients,\n"
        "\tpositive leading coefficient and content 1, having the algebraic\n"
        "\tnumber s as a root.  s may be built from rationals, radicals, the\n"
        "\timaginary unit, roots of unity, and Root[] objects.\n"
        "MinimalPolynomial[s]\n"
        "\tgives the minimal polynomial as a pure function.\n"
        "MinimalPolynomial[s, x, Extension -> a]\n"
        "\tgives the characteristic polynomial of s in Q(a) over Q(a).\n"
        "\tComputed by resultant elimination of the radicals; threads over\n"
        "\tlists.");
    symtab_set_docstring("Resultant",
        "Resultant[p, q, var]\n"
        "\tgives the resultant of p and q as polynomials in var: the unique\n"
        "\tinteger / polynomial scalar that vanishes iff p and q share a\n"
        "\troot in var.  Computed via a Sylvester-matrix determinant or, in\n"
        "\tthe exact path, a subresultant pseudo-remainder sequence.");
    symtab_set_docstring("Discriminant",
        "Discriminant[poly, var]\n"
        "\tgives the discriminant of poly with respect to var: up to sign and\n"
        "\tleading-coefficient scaling, Resultant[poly, D[poly, var], var] /\n"
        "\tlc[poly, var].  Vanishes iff poly has a repeated root in var.");
    symtab_set_docstring("Subresultants",
        "Subresultants[poly1, poly2, var]\n"
        "\tgives the list of principal subresultant coefficients of poly1 and\n"
        "\tpoly2 with respect to var.  The list has length\n"
        "\tMin[Exponent[poly1, var], Exponent[poly2, var]] + 1, its first\n"
        "\telement is Resultant[poly1, poly2, var], and the first k entries\n"
        "\tvanish exactly when the polynomials share k roots (with\n"
        "\tmultiplicity).  Computed by a subresultant polynomial-remainder\n"
        "\tsequence, or a Sylvester-minor determinant for algebraic\n"
        "\tcoefficients.");
    symtab_set_docstring("SubresultantPolynomials",
        "SubresultantPolynomials[poly1, poly2, var]\n"
        "\tgives the list of subresultant polynomials {S_0, ..., S_m} of\n"
        "\tpoly1 and poly2 with respect to var, where m = Exponent[poly2,\n"
        "\tvar].  The list has length m + 1, its first element is\n"
        "\tResultant[poly1, poly2, var], and the coefficient of var^j in S_j\n"
        "\tis the j-th principal subresultant coefficient.  Requires\n"
        "\tExponent[poly1, var] >= Exponent[poly2, var] and exact\n"
        "\tcoefficients.  Computed by a subresultant polynomial-remainder\n"
        "\tsequence, with a determinant-polynomial fallback for algebraic\n"
        "\tcoefficients.");
    symtab_set_docstring("GroebnerBasis",
        "GroebnerBasis[{p1, p2, ...}, {x1, x2, ...}]\n"
        "\tgives a list of polynomials that form a Gröbner basis for the\n"
        "\tideal generated by the pi over Q[x1, x2, ...].\n"
        "GroebnerBasis[{p1, ...}, {x1, ...}, {y1, ...}]\n"
        "\teliminates the yi and returns a Gröbner basis of the elimination\n"
        "\tideal in {x1, ...}.\n"
        "Free symbols in pi that are not in vars (and not a known constant\n"
        "\tlike Pi, E, EulerGamma, ...) are auto-promoted to coefficient-ring\n"
        "\tparameters and survive in the basis polynomials.\n"
        "Options: MonomialOrder -> Lexicographic (default) |\n"
        "\tDegreeReverseLexicographic | EliminationOrder | a weight matrix\n"
        "\t{{...}, ...} of integers (one column per variable, defining a\n"
        "\tcustom term order); CoefficientDomain -> Rationals (default);\n"
        "\tMethod -> \"Buchberger\" (default) | \"GroebnerWalk\"; Sort\n"
        "\t-> True reverses the main-variable list; ParameterVariables -> p\n"
        "\tor {p1, ...} marks parameters explicitly (the main-variable list\n"
        "\tis then optional and is auto-derived from the polys).  Other\n"
        "\tsettings emit GroebnerBasis::nimpl and fall back; Modulus emits\n"
        "\tGroebnerBasis::modnotimpl and falls back to the rational basis.\n"
        "Polynomial equations (Equal[a, b]) are accepted in place of\n"
        "\tpolynomials; each is rewritten as a - b before computation.\n"
        "TimeConstrained[GroebnerBasis[...], t] aborts cleanly via the\n"
        "\tcooperative deadline hook.  Lex order can be exponentially slow\n"
        "\ton systems with 3+ variables; switch to MonomialOrder ->\n"
        "\tDegreeReverseLexicographic for hard inputs.");
    symtab_set_docstring("Eliminate",
        "Eliminate[eqns, vars]\n"
        "\teliminates vars between a list/conjunction of simultaneous\n"
        "\tequations lhs == rhs, returning a balanced Equal[] or an\n"
        "\tAnd[] of Equal[]s in the remaining variables (True if the\n"
        "\telimination ideal is empty, False if the system is\n"
        "\tinconsistent). Works on polynomial equations over Q via a\n"
        "\tlexicographic Gröbner basis with elimination block. A\n"
        "\tprincipal-branch inverse-function pre-pass peels single Sin/\n"
        "\tCos/Tan/Sinh/Cosh/Tanh/Exp/Log wrappers and emits\n"
        "\tEliminate::ifun; non-polynomial systems otherwise return\n"
        "\tunevaluated with Eliminate::nlin.");
    symtab_set_docstring("SolveAlways",
        "SolveAlways[eqns, vars]\n"
        "\tfinds values of parameters appearing in eqns but not in vars\n"
        "\tsuch that eqns hold for every value of vars.\n"
        "Equations may be Equal[lhs, rhs] or a List/And of such.\n"
        "Reduction: each lhs - rhs is treated as a polynomial in vars\n"
        "\tvia CoefficientList; every coefficient must vanish, and the\n"
        "\tresulting system is passed to Solve with the remaining\n"
        "\tsymbols (the parameters) as unknowns.  Returns {} when there\n"
        "\tare no parameters to solve for.");
    symtab_set_docstring("ToRadicals",
        "ToRadicals[expr]\n"
        "\tattempts to express all Root objects in expr in terms of radicals.\n"
        "\n"
        "ToRadicals can always give expressions in terms of radicals when the\n"
        "\thighest degree of the polynomial that appears in any Root object is\n"
        "\tfour.  Binomial Root objects of the form Root[Function[a #^n + b], k]\n"
        "\tare also reduced to radicals for any degree n.  Other Root objects\n"
        "\tof degree five or higher are returned unchanged.\n"
        "If Root objects in expr contain parameters, ToRadicals[expr] may yield\n"
        "\ta result that is not equal to expr for all values of the parameters.\n"
        "ToRadicals automatically threads over lists, equations, inequalities,\n"
        "\tand logic functions.");
    symtab_set_docstring("Numerator",
        "Numerator[expr]\n"
        "\tgives the numerator of expr regarded as a rational expression.\n"
        "\tPicks out factors of expr that do not carry a superficially negative\n"
        "\texponent; constants and symbols pass through as-is.");
    symtab_set_docstring("Denominator",
        "Denominator[expr]\n"
        "\tgives the denominator of expr regarded as a rational expression.\n"
        "\tCollects factors of expr that carry a superficially negative\n"
        "\texponent, inverted; returns 1 when no such factors exist.");
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
        "\n"
        "Options:\n"
        "  Assumptions (default $Assumptions) -- facts assumed while simplifying.\n"
        "  ComplexityFunction (default: leaf count plus integer-digit count, matching Mathematica) -- ranks candidate forms; the lowest-scoring form is returned.\n"
        "  TransformationFunctions (default Automatic) -- the functions applied to try to transform parts of expr. Automatic uses the built-in collection; {f1, f2, ...} uses only the fi; {Automatic, f1, ...} uses the built-in functions together with the fi.\n"
        "\n"
        "The built-in collection tries Together, Cancel, Expand, Factor, FactorSquareFree, Apart, TrigExpand, TrigFactor, and a TrigToExp/ExpToTrig roundtrip, keeping the smallest result.\n"
        "Under positivity / reality assumptions Simplify also applies Log/Power identities -- Log[a b] -> Log[a] + Log[b], (a b)^c -> a^c b^c, (a^p)^q -> a^(p q), Log[a^p] -> p Log[a] and the like -- whenever the operand-domain conditions are provable from the assumption set.\n"
        "Assumptions can be equations, inequalities, domain specifications such as Element[x, Integers], or logical combinations of these. Lists of assumptions are converted to conjunctions.\n"
        "Simplify automatically threads over lists, equations, inequalities, and logic functions.");
    symtab_set_docstring("TransformationFunctions",
        "TransformationFunctions\n\tis an option for Simplify giving the list of functions to apply to try to transform parts of an expression.\n"
        "TransformationFunctions -> Automatic uses the built-in collection of transformation functions.\n"
        "TransformationFunctions -> {f1, f2, ...} uses only the functions fi.\n"
        "TransformationFunctions -> {Automatic, f1, ...} uses the built-in transformation functions together with the fi.\n"
        "Each function is applied to the whole expression and to its subexpressions; the lowest-complexity result (per ComplexityFunction) is kept.");
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
        "Element[{x1, ..., xN}, dom] and Element[x1 | ... | xN, dom] are shorthand for the conjunction Element[x1, dom] && ... && Element[xN, dom]: True/False if every component decides, otherwise unevaluated and treated as a joint per-variable fact by Simplify.");
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
    symtab_set_docstring("SameQ",
        "lhs === rhs or SameQ[lhs, rhs]\n"
        "\tyields True if lhs and rhs are structurally identical (head-by-head,\n"
        "\targument-by-argument), and False otherwise.  Numerically equal but\n"
        "\tdistinct heads (e.g. 1 and 1.) are NOT considered same.");
    symtab_set_docstring("UnsameQ",
        "lhs =!= rhs or UnsameQ[lhs, rhs]\n"
        "\tis the negation of SameQ: True iff lhs and rhs are not structurally\n"
        "\tidentical.");
    symtab_set_docstring("Equal",
        "lhs == rhs or Equal[lhs, rhs]\n"
        "\ttests mathematical equality. Numeric arguments decide directly\n"
        "\t(Integer / Rational exact comparison; Real / MPFR comparison with\n"
        "\tprecision tolerance); structurally identical symbolic forms decide\n"
        "\tTrue; otherwise the call stays unevaluated as a symbolic equation.\n"
        "Equal threads over Lists pairwise; chained Equal becomes Inequality.");
    symtab_set_docstring("Unequal",
        "lhs != rhs or Unequal[lhs, rhs]\n"
        "\tis the negation of Equal: True if lhs and rhs can be decided unequal,\n"
        "\tFalse if they can be decided equal, otherwise unevaluated.");
    symtab_set_docstring("Less",
        "x < y or Less[x, y]\n"
        "\tyields True if x is strictly less than y on numeric inputs, False\n"
        "\tif strictly greater or equal, otherwise unevaluated.\n"
        "Chained forms (x < y < z) become Inequality, decided pairwise.");
    symtab_set_docstring("Greater",
        "x > y or Greater[x, y]\n"
        "\tyields True if x is strictly greater than y on numeric inputs,\n"
        "\tFalse if strictly less or equal, otherwise unevaluated.\n"
        "Chained forms become Inequality.");
    symtab_set_docstring("LessEqual",
        "x <= y or LessEqual[x, y]\n"
        "\tyields True if x is less than or equal to y on numeric inputs,\n"
        "\tFalse if strictly greater, otherwise unevaluated.");
    symtab_set_docstring("GreaterEqual",
        "x >= y or GreaterEqual[x, y]\n"
        "\tyields True if x is greater than or equal to y on numeric inputs,\n"
        "\tFalse if strictly less, otherwise unevaluated.");
    symtab_set_docstring("Inequality",
        "Inequality[v0, op0, v1, op1, v2, ...] is the canonical form for a "
        "chained comparison such as a < b <= c. It returns True if every "
        "adjacent pair holds, False if any pair fails, and otherwise the "
        "residual chain with decidable pairs dropped.");

    symtab_set_docstring("IntegerDigits",
        "IntegerDigits[n] gives a list of the decimal digits in the integer n.\n"
        "IntegerDigits[n, b] gives a list of the base b digits in the integer n.\n"
        "IntegerDigits[n, b, len] pads the list on the left with zeros to give a list of length len; if n has more than len base-b digits, the last len least-significant digits are returned.");

    symtab_set_docstring("IntegerLength",
        "IntegerLength[n] gives the number of decimal digits in the integer n.\n"
        "IntegerLength[n, b] gives the number of base b digits in n.\n"
        "IntegerLength ignores the sign of n; IntegerLength[0] is 0.");

    symtab_set_docstring("IntegerExponent",
        "IntegerExponent[n, b] gives the highest power of b that divides n.\n"
        "IntegerExponent[n] is equivalent to IntegerExponent[n, 10] and "
        "gives the number of trailing zeros in the decimal digits of n.\n"
        "IntegerExponent ignores the sign of n; IntegerExponent[0, b] is "
        "Infinity.");

    symtab_set_docstring("DigitCount",
        "DigitCount[n] gives a list of the counts of digits 1, 2, ..., 9, 0 "
        "in the base-10 representation of n.\n"
        "DigitCount[n, b] gives a list of the counts of digits 1, 2, ..., "
        "b-1, 0 in the base-b representation of n.\n"
        "DigitCount[n, b, d] gives the number of d digits in the base-b "
        "representation of n.\n"
        "The sign of n is discarded; DigitCount[0] is a list of zeros.");

    symtab_set_docstring("IntegerString",
        "IntegerString[n] gives a string consisting of the decimal digits "
        "in the integer n.\n"
        "IntegerString[n, b] gives a string consisting of the base-b digits "
        "in n; digit values 10 to 35 use the letters a-z.\n"
        "IntegerString[n, b, len] pads the string on the left with zero "
        "digits to give a string of length len; if len is less than the "
        "number of digits in n, the len least-significant digits are "
        "returned.\n"
        "The maximum allowed base is 36; the sign of n is discarded.");

    symtab_set_docstring("RealDigits",
        "RealDigits[x] gives a list {digits, exp} of the digits in the "
        "approximate real number x together with the exponent such that "
        "the first digit is the coefficient of 10^(exp - 1).\n"
        "RealDigits[x, b] gives base-b digits.\n"
        "RealDigits[x, b, len] gives len digits.\n"
        "RealDigits[x, b, len, n] gives len digits starting from the "
        "coefficient of b^n.\n"
        "For rationals with non-terminating expansions the digit list "
        "ends in a nested list of the recurring block.  For inexact "
        "(machine or MPFR) reals, digits beyond the available precision "
        "are returned as Indeterminate.  The sign of x is discarded.");

    symtab_set_docstring("MantissaExponent",
        "MantissaExponent[x] gives a list {m, e} containing the mantissa "
        "and exponent of the real number x, such that x = m * 10^e and "
        "1/10 <= |m| < 1 (or m = 0 when x = 0).\n"
        "MantissaExponent[x, b] gives the base-b mantissa and exponent; "
        "the mantissa then lies in 1/b <= |m| < 1.\n"
        "Works for exact (Integer, Rational) and approximate (Real, MPFR) "
        "numeric inputs.  For exact inputs the mantissa is an exact "
        "Rational; for inexact inputs the mantissa carries the same "
        "precision as x.  Currently only integer bases >= 2 are supported.");

    symtab_set_docstring("RealExponent",
        "RealExponent[x] gives Log[10, |x|] -- the base-10 real exponent "
        "of x.\n"
        "RealExponent[x, b] gives Log[b, |x|] in the specified base b.\n"
        "Accepts Integer, BigInt, Rational, Real, and (with USE_MPFR) MPFR "
        "inputs, plus symbolic numeric values such as Pi, E, or Pi^Pi.  "
        "Result is a machine Real unless an MPFR input lifts it to MPFR "
        "at that precision.  Exact zero gives -Infinity; machine 0. gives "
        "Log[b, $MinMachineNumber] (~ -307.65 in base 10); MPFR 0 with "
        "precision p digits gives -p / Log10[b].  Threads over lists.");

    symtab_set_docstring("FromDigits",
        "FromDigits[list] constructs an integer from a list of decimal "
        "digits, most-significant first.\n"
        "FromDigits[list, b] takes the digits to be given in base b.\n"
        "FromDigits[\"string\"] constructs an integer from a string of "
        "digits, where letters a-z and A-Z denote digit values 10-35.\n"
        "FromDigits[\"string\", b] takes the digits in the string to be "
        "given in base b.\n"
        "Digits in list and characters in the string need not be less than "
        "the base; they are carried through Horner's method.  Symbolic "
        "digits or base expand to the polynomial sum d[0] b^(n-1) + ... + "
        "d[n-1].");

    // Primes
    symtab_set_docstring("FactorInteger", "FactorInteger[n] gives a list of the prime factors of the integer n, together with their exponents.");
    symtab_set_docstring("EulerPhi", "EulerPhi[n] gives the Euler totient function phi(n).");
    symtab_set_docstring("PrimePi", "PrimePi[x] gives the number of primes less than or equal to x.");
    symtab_set_docstring("NextPrime", "NextPrime[x] gives the next prime after x.");
    symtab_set_docstring("Distribute",
        "Distribute[f[x1, x2, ...]]\n"
        "\tdistributes f over Plus appearing in any of the xi, building the sum\n"
        "\tof f applied to every Cartesian-product selection of summands.\n"
        "Distribute[expr, g]\n"
        "\tdistributes over the head g instead of Plus.\n"
        "Distribute[expr, g, f]\n"
        "\tperforms the distribution only if the head of expr is f.\n"
        "Distribute[expr, g, f, gp, fp]\n"
        "\tgives gp and fp in place of g and f respectively in the result.");
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

    // System floating-point constants (registered in core.c).
    symtab_set_docstring("$MachinePrecision",
        "$MachinePrecision\n"
        "\tgives the number of decimal digits of precision used for\n"
        "\tmachine-precision numbers.\n"
        "\n"
        "Derived from the platform's DBL_MANT_DIG -- typically 53*Log[10,2]\n"
        "(~ 15.9546) on IEEE 754 systems.");
    symtab_set_docstring("$MachineEpsilon",
        "$MachineEpsilon\n"
        "\tgives the difference between 1.0 and the next-nearest number\n"
        "\trepresentable as a machine-precision number.\n"
        "\n"
        "Equals the platform's DBL_EPSILON; measures the granularity of\n"
        "machine-precision numbers.");
    symtab_set_docstring("$MinMachineNumber",
        "$MinMachineNumber\n"
        "\tgives the smallest positive machine-precision number that can\n"
        "\tbe represented in normalized form on this computer system.\n"
        "\n"
        "Equals the platform's DBL_MIN.");
    symtab_set_docstring("$MaxMachineNumber",
        "$MaxMachineNumber\n"
        "\tgives the largest machine-precision number that can be used on\n"
        "\tthis computer system.\n"
        "\n"
        "Equals the platform's DBL_MAX.");
    symtab_set_docstring("$MaxNumber",
        "$MaxNumber\n"
        "\tgives the maximum arbitrary-precision number that can be\n"
        "\trepresented on this computer system.\n"
        "\n"
        "With USE_MPFR builds, this is the largest finite value at machine\n"
        "precision under MPFR's current exponent range; otherwise it equals\n"
        "$MaxMachineNumber.");
    symtab_set_docstring("$MinNumber",
        "$MinNumber\n"
        "\tgives the minimum positive arbitrary-precision number that can\n"
        "\tbe represented on this computer system.\n"
        "\n"
        "With USE_MPFR builds, this is the smallest positive value at\n"
        "machine precision under MPFR's current exponent range; otherwise\n"
        "it equals $MinMachineNumber.");

    // Graphics options
    symtab_set_docstring("AspectRatio",
        "AspectRatio\n"
        "\tis an option for Graphics and Plot that specifies the ratio of\n"
        "\theight to width of the rendered plot.\n"
        "\n"
        "AspectRatio -> Automatic sets the ratio from the actual coordinate\n"
        "values (true geometry); AspectRatio -> Full stretches the graphics to\n"
        "fill the enclosing region; AspectRatio -> a uses the explicit\n"
        "height-to-width ratio a. Plot defaults to 1/GoldenRatio.");
    symtab_set_docstring("ImageSize",
        "ImageSize\n"
        "\tis an option for Graphics and Plot that specifies the overall size\n"
        "\tof the image to display.\n"
        "\n"
        "ImageSize -> w sets the width to w pixels, with the height following\n"
        "from AspectRatio; ImageSize -> {w, h} fixes both the width and height\n"
        "in pixels. The default is a width of 800.");
    symtab_set_docstring("Frame",
        "Frame\n"
        "\tis an option for Graphics and Plot that specifies whether to draw a\n"
        "\tframe with ticks and labels around the plot.\n"
        "\n"
        "Frame -> True boxes all four edges; Frame -> False (or None) draws no\n"
        "frame; Frame -> {{left, right}, {bottom, top}} toggles each edge with\n"
        "True or False. In Plot a frame takes the place of the default Axes.");
    symtab_set_docstring("FrameTicks",
        "FrameTicks\n"
        "\tis an option for Graphics and Plot that specifies the tick marks on\n"
        "\tthe edges of a frame.\n"
        "\n"
        "FrameTicks -> Automatic (the default) draws major and minor ticks with\n"
        "labels on every drawn frame edge; FrameTicks -> None keeps the frame\n"
        "box but draws no ticks; the {{left, right}, {bottom, top}} form selects\n"
        "Automatic or None per edge.");
    symtab_set_docstring("FrameStyle",
        "FrameStyle\n"
        "\tis an option for Graphics and Plot that specifies the style of the\n"
        "\tframe box, its ticks and its labels.\n"
        "\n"
        "FrameStyle -> RGBColor[...] or GrayLevel[...] sets the colour. The\n"
        "default is a neutral gray.");
}
