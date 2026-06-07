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

-   :material-decimal: __[4. Machine & arbitrary precision](04-machine-and-arbitrary-precision-arithmetic.md)__

    Exact integers and rationals, fast machine-precision reals, and
    arbitrary-precision arithmetic. Meet `N`, `Precision`, and when to use each.

-   :material-function-variant: __[5. Algebra](05-algebra.md)__

    Expand and factor polynomials, dissect and divide them, reshape rational
    expressions with `Together`/`Apart`, simplify, and solve equations.

-   :material-math-integral: __[6. Calculus](06-calculus.md)__

    Differentiate and integrate, expand power series, take limits, evaluate
    symbolic sums, and find roots and extrema numerically.

</div>

!!! tip "Following along"
    Start the REPL with `./Mathilda` and type each `In[...]` line yourself
    (without the prompt). Press Return to evaluate. End a line with `\` to
    continue a long expression onto the next line.
