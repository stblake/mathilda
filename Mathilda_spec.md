# Mathilda Specification

Mathilda is a tiny, AI Agent-generated, symbolic computer algebra system (CAS) inspired by the core architecture of Mathematica. Mathilda is written in the C programming language. It supports a recursive expression model, pattern matching, rewriting rules, and a small library of built-in functions.

The name, Mathilda is inspired by David Stoutemyer's PICOMATH-80 tiny computer algebra system.

This file is the entry point. Detailed material lives under [`docs/spec/`](docs/spec/).

## Foundations

| Topic | File |
|-------|------|
| Expression structure (atoms, function nodes, types) | [`docs/spec/expression-structure.md`](docs/spec/expression-structure.md) |
| Operators and precedence | [`docs/spec/operators.md`](docs/spec/operators.md) |
| Known limitations | [`docs/spec/limitations.md`](docs/spec/limitations.md) |
| Technical implementation details | [`docs/spec/architecture.md`](docs/spec/architecture.md) |
| Performance notes (poly, deriv, limit) | [`docs/spec/performance.md`](docs/spec/performance.md) |

For a deeper architectural treatment (subsystem boundaries, evaluation pipeline, memory ownership rules, extension recipes), see the top-level [`SPEC.md`](SPEC.md).

## Built-in function reference

Each category lives in [`docs/spec/builtins/`](docs/spec/builtins/):

| Category | File |
|----------|------|
| Structural manipulation (`Part`, `Head`, `First`, `Insert`, `Span`, ...) | [`builtins/structural-manipulation.md`](docs/spec/builtins/structural-manipulation.md) |
| Expression information (`AtomQ`, `Length`, `Depth`, `LeafCount`, ...) | [`builtins/expression-information.md`](docs/spec/builtins/expression-information.md) |
| Time and date (`Timing`, `RepeatedTiming`, ...) | [`builtins/time-and-date.md`](docs/spec/builtins/time-and-date.md) |
| Linear algebra (`Dot`, `Det`, `Inverse`, `Cross`, `Norm`, `Eigenvalues`, `Eigenvectors`, ...) | [`builtins/linear-algebra.md`](docs/spec/builtins/linear-algebra.md) |
| Statistics (`Mean`, `Variance`, `Median`, ...) | [`builtins/statistics.md`](docs/spec/builtins/statistics.md) |
| Random number generation | [`builtins/random-number-generation.md`](docs/spec/builtins/random-number-generation.md) |
| String operations | [`builtins/string-operations.md`](docs/spec/builtins/string-operations.md) |
| Arithmetic (`Plus`, `Times`, `Power`, `Mod`, `Factorial`, `Binomial`, ...) | [`builtins/arithmetic.md`](docs/spec/builtins/arithmetic.md) |
| Number theory (`GCD`, `LCM`, `PowerMod`, `PrimeQ`, `FactorInteger`, `EulerPhi`, ...) | [`builtins/number-theory.md`](docs/spec/builtins/number-theory.md) |
| Algebra (`Factor`, `Expand`, `Together`, `Apart`, `GroebnerBasis`, ...) | [`builtins/algebra.md`](docs/spec/builtins/algebra.md) |
| Solutions of equations (`Solve`, `SolveAlways`, `Root`, `ToRadicals`, `Eliminate`, ...) | [`builtins/solutions-of-equations.md`](docs/spec/builtins/solutions-of-equations.md) |
| Comparisons (`Equal`, `Less`, `Greater`, `SameQ`, `Inequality`, ...) | [`builtins/comparisons.md`](docs/spec/builtins/comparisons.md) |
| Calculus (`D`, `Integrate`, `Limit`, ...) | [`builtins/calculus.md`](docs/spec/builtins/calculus.md) |
| Simplification (`Simplify`, `SimplifyCount`, `Assuming`, `$Assumptions`, `Element`, ...) | [`builtins/simplification.md`](docs/spec/builtins/simplification.md) |
| Power series (`Series`, `SeriesData`, `Normal`, ...) | [`builtins/power-series.md`](docs/spec/builtins/power-series.md) |
| Elementary functions (trig, hyperbolic, log, exp) | [`builtins/elementary-functions.md`](docs/spec/builtins/elementary-functions.md) |
| Mathematical constants (`Pi`, `E`, `Degree`, `EulerGamma`, `Catalan`, `GoldenRatio`, `Glaisher`, `Khinchin`, ...) | [`builtins/mathematical-constants.md`](docs/spec/builtins/mathematical-constants.md) |
| Special functions (`Gamma`, `Pochhammer`, `HypergeometricPFQ`, `Hypergeometric0F1`/`1F1`/`2F1`, ...) | [`builtins/special-functions.md`](docs/spec/builtins/special-functions.md) |
| Numerical calculus (`ND`, `NIntegrate`, `NSum`, `NProduct`, `NLimit`, `NSeries`, `NResidue`) | [`builtins/numerical-calculus.md`](docs/spec/builtins/numerical-calculus.md) |
| Lists and iteration (`Table`, `Range`, `Map`, `Do`, ...) | [`builtins/lists-and-iteration.md`](docs/spec/builtins/lists-and-iteration.md) |
| Data structures (`Association`/`<\|...\|>`, `Keys`, `Values`, `Lookup`, `Counts`, `GroupBy`, `Merge`, ...) | [`builtins/data-structures.md`](docs/spec/builtins/data-structures.md) |
| Functional programming (`Function`, `Apply`, `Select`, ...) | [`builtins/functional-programming.md`](docs/spec/builtins/functional-programming.md) |
| Control flow (`If`, `Which`, `Switch`, `For`, `While`) | [`builtins/control-flow.md`](docs/spec/builtins/control-flow.md) |
| Assignment and rules (`Set`, `SetDelayed`, `Rule`, `RuleDelayed`) | [`builtins/assignment-and-rules.md`](docs/spec/builtins/assignment-and-rules.md) |
| Scoping constructs (`Module`, `Block`, `With`) | [`builtins/scoping-constructs.md`](docs/spec/builtins/scoping-constructs.md) |
| Pattern matching (`MatchQ`, `Cases`, `DeleteCases`, `Position`, `Count`, ...) | [`builtins/pattern-matching.md`](docs/spec/builtins/pattern-matching.md) |
| File I/O (`Get`, `Put`, `PutAppend`, `>>`, `>>>`) | [`builtins/file-io.md`](docs/spec/builtins/file-io.md) |
| Graphics (`Graphics`, `Show`, `Plot`, `Point`, `Line`, `Rectangle`, `Circle`, `Disk`, `Polygon`, `Text`, ...) | [`builtins/graphics.md`](docs/spec/builtins/graphics.md) |
| FLINT context (`` FLINT`PolynomialGCD ``, `` FLINT`Factor ``, `` FLINT`Det ``, `` FLINT`Zeta ``, ...) — direct access to the FLINT-backed kernels | [`builtins/flint.md`](docs/spec/builtins/flint.md) |

## Changelog

Detailed feature-addition and bug-fix notes, organized by week (Mon – Sun, keyed by Monday's date), in [`docs/spec/changelog/`](docs/spec/changelog/):

| Week (Mon – Sun) | File |
|------------------|------|
| 2026-04-20 → 2026-04-26 | [`changelog/2026-04-20.md`](docs/spec/changelog/2026-04-20.md) |
| 2026-04-27 → 2026-05-03 | [`changelog/2026-04-27.md`](docs/spec/changelog/2026-04-27.md) |
| 2026-05-04 → 2026-05-10 | [`changelog/2026-05-04.md`](docs/spec/changelog/2026-05-04.md) |
| 2026-05-11 → 2026-05-17 | [`changelog/2026-05-11.md`](docs/spec/changelog/2026-05-11.md) |
| 2026-05-18 → 2026-05-24 | [`changelog/2026-05-18.md`](docs/spec/changelog/2026-05-18.md) |
| 2026-05-25 → 2026-05-31 | [`changelog/2026-05-25.md`](docs/spec/changelog/2026-05-25.md) |
| 2026-06-01 → 2026-06-07 | [`changelog/2026-06-01.md`](docs/spec/changelog/2026-06-01.md) |
| 2026-06-08 → 2026-06-14 | [`changelog/2026-06-08.md`](docs/spec/changelog/2026-06-08.md) |
| 2026-06-15 → 2026-06-21 | [`changelog/2026-06-15.md`](docs/spec/changelog/2026-06-15.md) |
| 2026-06-22 → 2026-06-28 | [`changelog/2026-06-22.md`](docs/spec/changelog/2026-06-22.md) |
| 2026-06-29 → 2026-07-05 | [`changelog/2026-06-29.md`](docs/spec/changelog/2026-06-29.md) |

New entries land in the file for the current week (use the Monday-date of that week as the filename, format `YYYY-MM-DD.md`). When a change touches a built-in's documented behavior, the corresponding `docs/spec/builtins/*.md` file is updated as well; the changelog records the rationale and timing.
