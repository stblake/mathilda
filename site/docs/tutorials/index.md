# Tutorials

A guided path through Mathilda, from your first REPL session to writing your own
pattern-based rules and doing symbolic calculus. Every example is worked end to
end and **verified against the current Mathilda build**.

Work through them in order if you're new — each one builds on the last.

<div class="grid cards" markdown>

-   :material-rocket-launch: __[1. Getting started](01-getting-started.md)__

    Build Mathilda, launch the REPL, understand `In[]`/`Out[]`, learn the
    surface syntax, and get help on any function with `?Name`.

-   :material-file-tree: __[2. Expressions & evaluation](02-expressions-and-evaluation.md)__

    Everything is an expression. Meet `FullForm`, `Head`, the attribute system,
    the fixed-point evaluator, and how `Hold` suspends evaluation.

-   :material-vector-polyline: __[3. Pattern matching & rules](03-pattern-matching-and-rules.md)__

    Blanks and named patterns, conditions and tests, transformation rules
    (`->`, `:>`), replacement (`/.`, `//.`), and defining your own functions.

-   :material-decimal: __[4. Arithmetic](04-arithmetic.md)__

    Exact integers and rationals, fast machine-precision reals, and
    arbitrary-precision arithmetic (`N`, `Precision`); the basic operators,
    digit and radix manipulation, and combinatorial functions.

-   :material-key-variant: __[5. Number theory](05-number-theory.md)__

    `GCD`, `ExtendedGCD`, modular arithmetic and `PowerMod`, primes
    (`PrimeQ`, `FactorInteger`, `NextPrime`), `EulerPhi`, and continued
    fractions — up to RSA-style worked examples.

-   :material-function-variant: __[6. Algebra](06-algebra.md)__

    Expand and factor polynomials, dissect and divide them, reshape rational
    expressions with `Together`/`Apart`, simplify, and put the polynomial
    toolkit (`Resultant`, `GroebnerBasis`) to work on real problems.

-   :material-equal: __[7. Solutions of equations](07-solutions-of-equations.md)__

    Solve polynomial, transcendental, and simultaneous equations with `Solve`;
    `Root` objects and `ToRadicals`; eliminate variables with `Eliminate`; and
    tackle geometry and optimisation problems.

-   :material-math-integral: __[8. Calculus](08-calculus.md)__

    Differentiate and integrate, expand power series, take limits, evaluate
    symbolic sums, and find roots and extrema numerically.

-   :material-calculator-variant: __[9. Numerical calculus](09-numerical-calculus.md)__

    When there is no closed form: numerical integration, differentiation,
    summation, products, limits, series, and residues — `NIntegrate`, `ND`,
    `NSum`, `NProduct`, `NLimit`, `NSeries`, `NResidue`.

-   :material-sigma-lower: __[10. Special functions](10-special-functions.md)__

    The higher transcendental functions: `Gamma`, `Zeta`, `PolyGamma`, `Erf`,
    `PolyLog`, the Bernoulli and Euler numbers, and the hypergeometric family —
    with their exact reductions and numerical values.

</div>

!!! tip "Following along"
    Start the REPL with `./Mathilda` and type each `In[...]` line yourself
    (without the prompt). Press Return to evaluate. End a line with `\` to
    continue a long expression onto the next line.
