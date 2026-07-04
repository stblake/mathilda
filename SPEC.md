# Mathilda System Architecture

## 1. Overview

Mathilda is a symbolic computer algebra system written in C99, modelled on the
core architecture and evaluation semantics of Mathematica (the Wolfram
Language). It implements a recursive expression model, structural pattern
matching with backtracking, rewriting rules, and an extensive library of
built-in mathematical functions. The name pays homage to David Stoutemyer's
PICOMATH-80.

**Key characteristics:**
- ~22 kLoC of C99 across ~44 source modules.
- Arbitrary-precision integers via GMP.
- Interactive REPL with GNU Readline support.
- Mathematica-compatible surface syntax and evaluation semantics.
- Licensed under GPLv3.

**External dependencies:**
- **GMP** — arbitrary-precision integers.
- **GNU Readline** — interactive line editing and history.
- **GMP-ECM** — advanced integer factorization (vendored in `src/external/ecm/`,
  do not modify).
- **Raylib** (optional, `USE_GRAPHICS`) — windowing/rendering backend for
  `Graphics[]`/`Show[]`/`Plot[]`. Autodetected via `pkg-config`; the build
  degrades gracefully to a text placeholder when absent.

**Companion documents:**
- [`Mathilda_spec.md`](Mathilda_spec.md) — built-in function reference index.
- [`docs/spec/`](docs/spec/) — per-category built-in docs and changelog.
- [`docs/extending.md`](docs/extending.md) — how to add builtins, modules,
  patterns, internal rules, and operators.
- [`CLAUDE.md`](CLAUDE.md) — contributor workflow.

---

## 2. Project Layout

```
src/
  expr.{c,h}        Expression representation (tagged union AST)
  parse.{c,h}       Pratt parser
  eval.{c,h}        Core evaluator (fixed-point loop)
  symtab.{c,h}      Symbol table: OwnValues, DownValues, builtins, attributes
  match.{c,h}       Pattern matcher with backtracking
  replace.{c,h}     Rule engine (Replace, ReplaceAll, ReplaceRepeated)
  attr.{c,h}        Attribute bitflags
  print.{c,h}       Expression formatting
  repl.c            REPL entry point + main()
  core.{c,h}        Core builtins and module initialization hub
  plus.c / times.c / power.c / arithmetic.c   Arithmetic heads
  list.c / part.c / sort.c                    List + structural ops
  trig.c / hyperbolic.c / logexp.c / complex.c
  comparisons.c / boolean.c / cond.c / iter.c
  funcprog.c / purefunc.c / patterns.c
  rat.c / parfrac.c / expand.c
  modular.c / facint.c / piecewise.c / stats.c
  load.c / info.c / datetime.c
  numbertheory/     GCD/LCM/ExtendedGCD/PowerMod/Factorial/Binomial/
                    PrimitiveRoot/MultiplicativeOrder/Divisors/DivisorSigma/...
                    (one builtin per file; shared helpers in nt_util.c /
                    nt_gaussian.c, declared in numbertheory_internal.h)
  poly/             Polynomial subsystem (univariate, multivariate, factoring,
                    algebraic-number factoring, polynomial solving)
  linalg/           Dense linear algebra; eigen kernels split by algorithm
  graph/            Graph subsystem: Graph[] construction/validation, queries
                    (VertexList/EdgeList/degree/AdjacencyList), matrix views
                    (Adjacency/Incidence/Kirchhoff/Distance), generators
                    (Complete/Cycle/Path/Petersen/Random/...), algorithms
                    (shortest path, components, spanning tree, connectivity,
                    Euler/Hamilton, cliques, colouring), centrality (PageRank,
                    Katz, betweenness, clustering) and GraphPlot/Graph3D
                    (one builtin per file, mirrors linalg/)
  simp/             Simplify, trig simplification, trig rationalisation
  calculus/         D / Dt / Derivative, Series, Limit, Integrate, Risch-Norman
  special_functions/ Higher transcendental & special functions: Gamma, LogGamma,
                    PolyGamma, Pochhammer, Zeta, StieltjesGamma, EulerGamma,
                    BernoulliB, EulerE, PolyLog, HypergeometricPFQ, Fibonacci,
                    LucasL
  graphics/         2D graphics engine: Graphics[]/Show[]/Plot[] primitives,
                    adaptive sampler, Raylib renderer (USE_GRAPHICS), vector font
  internal/         .m bootstrap scripts (init.m, deriv.m, integral tables)
  external/ecm/     Vendored GMP-ECM (DO NOT MODIFY)
tests/              CMake-built unit suite (test_*.c)
makefile            Primary build system
```

The full per-category built-in reference is in [`docs/spec/builtins/`](docs/spec/builtins/).

---

## 3. Core Architecture

The pipeline:

```
String ─► Parser ─► Evaluator ─► Printer ─► String
              │
              ├─ Symbol Table   (OwnValues, DownValues, builtins)
              ├─ Attribute Sys  (Hold*, Flat, Orderless, Listable, ...)
              ├─ Pattern Matcher (structural + backtracking)
              ├─ Rule Engine    (Replace / ReplaceAll)
              └─ Built-in C library
```

### 3.1 Expression representation (`expr.{c,h}`)

Every value is an `Expr`, a tagged union of `EXPR_INTEGER` (`int64_t`),
`EXPR_REAL` (`double`), `EXPR_SYMBOL`, `EXPR_STRING`, `EXPR_FUNCTION`
(head + args), and `EXPR_BIGINT` (GMP `mpz_t`). Compound expressions like
`f[x, y]` use `EXPR_FUNCTION` with `head = Symbol("f")`; lists `{a, b}` are
`List[a, b]` internally. Bigints auto-promote from int64 on overflow and
demote back via `expr_bigint_normalize()` when they fit.

Expressions are immutable by convention during evaluation; transformations
return new trees.

| Function | Purpose |
|----------|---------|
| `expr_new_integer / _real / _symbol / _string / _function / _bigint_from_mpz` | Construct nodes |
| `expr_copy(e)` | Deep copy |
| `expr_free(e)` | Deep free |
| `expr_eq(a, b)` | Structural equality |
| `expr_compare(a, b)` | Canonical ordering |
| `expr_hash(e)` | Hash for tables |
| `expr_to_mpz(e, out)` | Coerce int-like to GMP |
| `expr_is_integer_like(e)` | True for `EXPR_INTEGER` or `EXPR_BIGINT` |

### 3.2 Parser (`parse.{c,h}`)

Pratt parser (top-down operator precedence). Lexing is inline — no separate
tokenizer. Precedences mirror Mathematica's standard values. Pattern syntax
(`_`, `__`, `___`, `x_`, `x_h`, `/;`) is recognised at parse time.

Selected precedences: `1000` `f[x]`, `730` blanks, `590` `^`, `400` `*`,
`310` `+ -`, `290` comparisons, `100` `[[…]]`, `90` `&`, `40` `= := -> :>`,
`10` `;`. The full table is in [`docs/spec/operators.md`](docs/spec/operators.md).

Public API: `parse_expression(input)` and `parse_next_expression(input_ptr)`.

### 3.3 Symbol table (`symtab.{c,h}`)

A `SymbolDef` carries the symbol's `own_values`, `down_values`, optional
builtin C function pointer, attribute bitflags, and docstring. Rules are
linked lists (newest first for DownValues).

- **OwnValues** — immediate symbol assignments (`x = 5`).
- **DownValues** — pattern rules on a function head (`f[x_] := x^2`).

| Function | Purpose |
|----------|---------|
| `symtab_init / symtab_clear` | Lifecycle |
| `symtab_get_def(name)` | Get or create a `SymbolDef` |
| `symtab_add_builtin(name, func)` | Register a C builtin |
| `symtab_set_docstring(name, doc)` | `?name` text |
| `symtab_add_own_value / _down_value` | Add a rule |
| `apply_own_values / apply_down_values` | Try matching rules during evaluation |

### 3.4 Attribute system (`attr.{c,h}`)

Bitflags on symbols that the evaluator consults before processing a call.
Selected flags:

| Flag | Effect |
|------|--------|
| `ATTR_HOLDFIRST / HOLDREST / HOLDALL / HOLDALLCOMPLETE` | Suppress arg evaluation |
| `ATTR_FLAT` | Associative: flatten nested same-head calls |
| `ATTR_ORDERLESS` | Commutative: canonical sort of args |
| `ATTR_LISTABLE` | Thread element-wise over `List` args |
| `ATTR_NUMERICFUNCTION` | Numeric operation |
| `ATTR_ONEIDENTITY` | `f[x]` ≡ `x` for pattern matching |
| `ATTR_PROTECTED` | Cannot be redefined |
| `ATTR_LOCKED / READPROTECTED / TEMPORARY / SEQUENCEHOLD / NHOLDREST` | Specialized |

Typical bundles: `Plus`/`Times` are `FLAT | LISTABLE | NUMERICFUNCTION |
ONEIDENTITY | ORDERLESS`; `Sin` is `LISTABLE | NUMERICFUNCTION | PROTECTED`;
`Table` is `HOLDALL | PROTECTED`.

### 3.5 Evaluator (`eval.{c,h}`)

Mathematica's infinite-evaluation semantics: `evaluate(e)` calls
`evaluate_step(e)` in a loop and stops when the expression stabilises (or a
recursion limit is hit).

For a function `f[a1, …, aN]` each step:

1. Evaluate the head `f`.
2. Read its attributes.
3. Evaluate arguments, except those suppressed by Hold flags.
4. **Listable** — thread over `List` arguments if `ATTR_LISTABLE`.
5. **Flat** — flatten nested same-head calls if `ATTR_FLAT`.
6. **Orderless** — sort args canonically if `ATTR_ORDERLESS`.
7. Call the registered C builtin if any. A builtin returns either a new
   `Expr*` (success) or `NULL` ("can't evaluate, leave alone").
8. Otherwise try `DownValues` — the first matching pattern's replacement wins.

For a bare symbol, step is: try `OwnValues`; return the replacement if any.

Helpers: `evaluate_step`, `eval_flatten_args(e, head_name)`,
`eval_sort_args(e)`.

### 3.6 Pattern matcher (`match.{c,h}`)

Structural tree unification with sequence matching and backtracking.

| Surface | Internal | Matches |
|---------|----------|---------|
| `_` | `Blank[]` | any single expression |
| `_h` | `Blank[h]` | any expression with head `h` |
| `__` / `___` | `BlankSequence[]` / `BlankNullSequence[]` | 1+ / 0+ expressions |
| `x_` / `x__` / `x___` | `Pattern[x, Blank…[]]` | bound version of the above |
| `p ? test` | `PatternTest[p, test]` | `p` if `test[match]` is `True` |
| `p /; cond` | `Condition[p, cond]` | `p` if `cond` evaluates to `True` |
| `p..` / `p...` | `Repeated[p]` / `RepeatedNull[p]` | sequence of `p`-matchers |
| `p:def` | `Optional[p, def]` | `p`, else `def` if absent |

Bindings live in a `MatchEnv`. API: `match(expr, pattern, env)`,
`replace_bindings(expr, env)`, `env_new / _free / _set / _get`. Sequence
backtracking handles different argument partitions when the pattern list
contains `__` or `___`.

### 3.7 Rule engine (`replace.{c,h}`)

- `Rule` (`->`) — immediate; RHS evaluated at rule creation.
- `RuleDelayed` (`:>`) — delayed; RHS evaluated after pattern binding.
- `ReplaceAll` (`/.`) — top-down traversal; matched subexpressions are not
  re-traversed into.
- `ReplaceRepeated` (`//.`) — `ReplaceAll` to a fixed point.
- `Replace` — applies rules only at specified levels.

---

## 4. Memory Management

Mathilda uses explicit manual memory management. The single most important
contract is the builtin ownership rule.

### Builtin ownership contract

```c
Expr* builtin_myfunc(Expr* res);
```

- The builtin **takes ownership** of `res`.
- Return a new `Expr*` on success — the **evaluator** frees `res`. Do not call
  `expr_free(res)` yourself.
- Return `NULL` if you cannot evaluate. The evaluator retains ownership and the
  expression remains unevaluated. Do **not** free `res` in this case.

When a builtin reuses parts of its input, NULL them out before the evaluator
frees the wrapper:

```c
Expr* arg0 = res->data.function.args[0];
res->data.function.args[0] = NULL;  /* prevent double-free */
return arg0;
```

### Safety checklist

- Never touch an `Expr*` after `expr_free()`.
- Never free the same `Expr*` twice — use the NULL-out-before-free pattern.
- Every `expr_new_*` or `expr_copy` you own must eventually be paired with
  `expr_free`.
- Run `valgrind --leak-check=full ./Mathilda` regularly.

---

## 5. Module Initialization

`main()` in `repl.c` calls `symtab_init()` then `core_init()`. `core_init()` in
`core.c` registers its own builtins and then calls each subsystem's
`*_init()` function: `parfrac`, `modular`, `facint`, `comparisons`, `boolean`,
`list`, `replace`, `patterns`, `cond`, `iter`, `complex`, `trig`, `simp`,
`hyperbolic`, `logexp`, `piecewise`, `attr`, `purefunc`, `stats`, `poly`,
`facpoly`, `rat`, `expand`, `info`, `datetime`, `linalg`, `load`, `graphics`,
`graph`.

After C-side init, `main()` loads `src/internal/init.m`, which `Get[]`s the
remaining `.m` bootstrap files.

---

## 6. Internal Definition Files

`src/internal/` contains Mathematica-syntax definition files loaded at startup:

- `init.m` — bootstrap, loads the other `.m` files via `Get[]`.
- `deriv.m` — derivative rules for `D[expr, x]`, defined as DownValues. Pure
  pattern-matching implementation of chain/product/quotient rules and all
  elementary derivatives.
- `CRCMathTablesIntegrals.m` — reference integral tables from the CRC
  Mathematical Tables.

This is the system's two-tier pattern: performance-sensitive primitives in C,
higher-level mathematical identities as rewrite rules in Mathilda's own
language.

---

## 7. Built-in Function Reference

The per-category built-in reference is the dedicated index, not this file.
Start from [`Mathilda_spec.md`](Mathilda_spec.md), which links into the
category-specific pages under [`docs/spec/builtins/`](docs/spec/builtins/).
Weekly changes (Mon – Sun, keyed by Monday's date) are recorded in
[`docs/spec/changelog/`](docs/spec/changelog/).

---

## 8. Extending Mathilda

Recipes for adding new builtins, modules, pattern constructs, internal `.m`
rules, and parser operators live in [`docs/extending.md`](docs/extending.md).

---

## 9. Build & Test

```bash
make -j$(nproc)        # builds ./Mathilda; -std=c99 -O3, links -lreadline -lgmp -lm
./Mathilda             # REPL

cd tests
mkdir -p build && cd build
cmake .. && make -j$(nproc)
for t in *_tests; do ./$t; done

valgrind --leak-check=full ./Mathilda
```

The makefile auto-discovers `src/*.c`. ECM is built and linked by default
(`USE_ECM=1`).

---

## 10. Coding Standards

- **C99 strictly.** No C11+ features. No POSIX-only types (`ssize_t`) or
  functions (`strdup`, `getline`, `asprintf`, `popen`, `fileno`, …) without
  C99-safe fallbacks. `<math.h>` constants like `M_PI` are POSIX, not C99 —
  guard with `#ifndef M_PI` fallbacks (see `src/trig.c`, `src/numeric.c`).
- **Memory safety.** Trace ownership; valgrind regularly. See §4.
- **No edits under `src/external/`.** ECM is vendored.
- **Docstrings.** Every builtin must have one via `symtab_set_docstring()`. Keep
  them terse — examples live in `docs/spec/...`.
- **Attributes.** Every builtin gets appropriate attributes in its module
  `_init()`.
- **Docs in sync.** Every new or modified builtin updates the matching file in
  `docs/spec/builtins/` and adds a changelog note under the current week's
  `docs/spec/changelog/<YYYY-MM-DD>.md` (where `<YYYY-MM-DD>` is the Monday
  of the ISO week the change lands in).

`CLAUDE.md` carries the full contributor workflow.

---

## 11. REPL

`repl.c` reads with GNU Readline (history, multiline via trailing `\`), parses
to an `Expr*`, evaluates to a fixed point, stores `In[n]`/`Out[n]` as
DownValues for back-reference, prints the result, and frees both trees. Exits
on `Quit[]` or EOF.

---

## 12. Printing (`print.{c,h}`)

The printer converts trees back into Mathematica-like syntax with
precedence-aware parenthesisation. `Plus[a, b]` prints as `a + b`,
`Times[a, b]` as `a b` or `a*b`, `Rational[1, 2]` as `1/2`, `Complex[a, b]` as
`a + b I`, `List[a, b]` as `{a, b}`. `FullForm[expr]` prints the raw tree;
`InputForm[expr]` prints in a form suitable for re-parsing.

---

## 13. Key Design Patterns

- **Everything is an expression.** Lists, matrices, rules, patterns are all
  `Expr*`. A matrix is a `List` of `List`s; a rule `a -> b` is `Rule[a, b]`;
  `x_` is `Pattern[x, Blank[]]`. Uniformity means the same generic tools
  (`Part`, `Map`, `ReplaceAll`, …) work on everything.
- **Attributes drive evaluation.** Generic evaluator + per-symbol attribute
  bits keeps the evaluator small and lets new builtins opt into listing,
  flattening, ordering, holding by setting a bit.
- **C for performance, rules for mathematics.** Hot paths (parser, evaluator,
  matcher, arithmetic) are C; higher-level identities (derivatives, integral
  tables) are DownValues in `.m` files.
- **Fixed-point evaluation.** Rule chains compose automatically: if
  `f[x_] := g[x]` and `g[x_] := x^2`, then `f[3]` evaluates to `9` over two
  passes, with no explicit chaining logic.
- **`NULL` means "I can't evaluate this".** Builtins return `NULL` for inputs
  they don't handle, letting symbolic arguments flow through the evaluator
  unchanged.
