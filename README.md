# Mathilda

Mathilda is a small, open source computer algebra system (CAS) heavily inspired by the core architecture and evaluation semantics of Mathematica (the Wolfram Language). Written entirely in C99 and its own language, it implements a recursive expression model, structural pattern matching with backtracking, rewriting rules, and an extensive library of built-in mathematical functions. 

Today Mathilda spans roughly **232,000 lines of C99** across **340 source modules**, exposing **~540 built-in functions** organized into **29 functional categories** — from arbitrary-precision arithmetic and symbolic calculus to polynomial factorization, dense linear algebra, and integer factorization.

## 🌟 Key Features

**Evaluation engine**
* **Infinite evaluation semantics:** expressions are repeatedly evaluated top-down until a fixed point is reached.
* **Attribute-driven evaluator:** a small generic core consults per-symbol bitflags (`HoldFirst`/`HoldAll`, `Flat`, `Orderless`, `Listable`, `OneIdentity`, `Protected`, …) to decide how to process each call.

**Pattern matching & rules**
* **First-class pattern matching:** `Blank` (`_`), `BlankSequence` (`__`), `BlankNullSequence` (`___`), named bindings (`x_`, `x_h`), `Condition` (`/;`), `PatternTest` (`?`), `Optional`, and `Repeated` — with full sequence backtracking.
* **Rule engine:** transformation rules (`->`, `:>`) and replacement operators (`/.`, `//.`, `Replace`).

**Numbers**
* **Arbitrary-precision integers** via the GNU Multiple Precision Arithmetic Library (GMP), with automatic promotion/demotion from machine integers.
* **Exact rationals and complex numbers**, plus **MPFR-backed arbitrary-precision reals** with precision/accuracy tracking (`N[expr, prec]`, precision literals such as `` 3.98`50 ``).

**Symbolic mathematics**
* **Calculus:** symbolic differentiation (`D`, `Dt`, `Derivative`); multi-method integration (`Integrate`) cascading rational-function, Risch–Norman, radical-substitution, derivative-divides, and CRC integral-table methods; `Series`, `Limit`, and symbolic summation (`Sum` via polynomial, geometric, and Gosper algorithms).
* **Polynomials:** univariate and multivariate arithmetic, factorization (square-free decomposition, Hensel lifting, irreducibility testing), algebraic-number factoring over ℚ(α), Gröbner bases, GCD/LCM, and partial fractions.
* **Linear algebra:** `Det`, `Inverse`, `Dot`, `Cross`; LU / QR / Cholesky / SVD / Schur decompositions; eigenvalues and eigenvectors via multiple algorithm kernels; norms, rank, and condition numbers — with optional LAPACK acceleration for machine-precision work.
* **Simplification:** `Simplify` with a complexity-driven search, trigonometric identities, radical denesting, and an assumptions framework (`$Assumptions`, `Element`).

**Number theory & factorization**
* **Number theory:** `GCD`, `LCM`, `ExtendedGCD`, `PowerMod`, `Divisors`, `EulerPhi`, `MoebiusMu`, `PrimitiveRoot`, continued fractions, and more.
* **Integer factorization:** a unified, automatic pipeline alongside explicit algorithms — Pollard's Rho, Pollard's $P-1$, Williams' $P+1$, Fermat, CFRAC, Dixon's Method, and the Elliptic Curve Method (ECM).

**Programming & utilities**
* **Functional programming:** `Map`, `Apply`, `Fold`, `Nest`, `Through`, `Composition`, and pure functions (`&` / `#`).
* **Scoping & control flow:** `Module`, `Block`, `With`; `If`, `Which`, `Switch`, `Do`, `For`, `While`, `Piecewise`.
* **Standard library:** lists and iteration, strings, statistics, random numbers, date/time, and file I/O.

---

## 📚 Function Categories

The complete reference (~540 functions) lives in [`Mathilda_spec.md`](Mathilda_spec.md), which indexes the per-category pages under [`docs/spec/builtins/`](docs/spec/builtins/):

* Arithmetic and Algebra
* Calculus
* Elementary Functions
* Linear Algebra
* Structural Manipulation
* Expression Information
* Functional Programming
* Lists and Iteration
* Pattern Matching
* Control Flow
* Assignment and Rules
* Scoping Constructs
* Simplification
* Power Series
* String Operations
* File I/O
* Statistics
* Random Number Generation
* Time and Date

Weekly change summaries are recorded under [`docs/spec/changelog/`](docs/spec/changelog/).

---

## 🚀 Getting Started

### Prerequisites

To build and run Mathilda you need:

* A C99-compliant compiler (`gcc` or `clang`)
* **GMP** (`libgmp` / `gmp-dev`) — arbitrary-precision integers *(required)*
* **GNU Readline** (`libreadline` / `readline-dev`) — interactive line editing *(required)*
* **MPFR** (`libmpfr` / `mpfr-dev`) — arbitrary-precision reals *(enabled by default)*
* **FLINT** ≥ 3.0 (`libflint` / `flint-dev`) — fast, rigorous polynomial arithmetic over algebraic extensions and rigorous `acb` numerics *(optional, auto-detected)*
* **GMP-ECM** (`gmp-ecm` / `libecm-dev`) — Elliptic Curve Method integer factorization *(optional, auto-detected)*
* **LAPACK / BLAS** — fast machine-precision linear algebra *(optional, auto-detected)*
* **CMake** — only required to build the test suite

The optional backends are controlled by build-time flags and **degrade gracefully** when disabled or absent:

| Flag | Default | Effect when on |
|------|---------|----------------|
| `USE_MPFR`  | `1` | Arbitrary-precision reals: `N[expr, prec]`, `Precision`/`Accuracy`, precision literals. Build without it via `make USE_MPFR=0`. |
| `USE_FLINT` | `1` | Fast, rigorous FLINT (≥ 3.0) kernels: multivariate polynomial GCD/factoring over ℚ, univariate GCD/factoring over number fields ℚ(α) (via the `gr` layer + ANTIC), the finite-field workhorse behind parametric ℚ(t)(α) work, and rigorous `acb` numerics (`Zeta`, `HurwitzZeta`, `PolyGamma`, `StieltjesGamma`). Auto-detected via `pkg-config` with a ≥ 3.0 version floor. Falls back to the classical (slower but still rigorous) path (`USE_FLINT=0`) when absent. |
| `USE_LAPACK`| `1` | Fast machine-precision linear algebra. Auto-detected: Apple **Accelerate** on macOS, `lapacke`/`lapack`/`blas` on Linux. Falls back to the pure-C path (`USE_LAPACK=0`) if none is found. |
| `USE_ECM`   | `1` | Elliptic Curve Method factorization via the system GMP-ECM library. Auto-detected via a compile-link probe; install `gmp-ecm` / `libecm-dev`. Falls back to disabled (`USE_ECM=0`) when absent. |

#### Installing dependencies

**Linux (Debian / Ubuntu):**

```bash
# Required libraries
sudo apt install libgmp-dev        # GMP — arbitrary-precision integers
sudo apt install libmpfr-dev       # MPFR — arbitrary-precision reals
sudo apt install libreadline-dev   # GNU Readline — interactive REPL

# Optional: FLINT (>= 3.0) for fast, rigorous algebraic-extension arithmetic
sudo apt install libflint-dev      # Debian Bookworm+/Ubuntu 24.04+ ship >= 3.0

# Optional: GMP-ECM for advanced integer factorization
sudo apt install libecm-dev

# Optional: LAPACK / BLAS for fast machine-precision linear algebra
sudo apt install liblapacke-dev libopenblas-dev

# Optional: CMake, only needed to build the test suite
sudo apt install cmake
```

On Fedora/RHEL the equivalents are `gmp-devel`, `mpfr-devel`, `readline-devel`,
`flint-devel` (≥ 3.0), `gmp-ecm-devel`, `lapack-devel`/`openblas-devel`, plus
`cmake`.

> **Note on FLINT versions.** Mathilda requires **FLINT ≥ 3.0** (the release that
> merged ANTIC for number-field arithmetic). Distributions that only package
> FLINT 2.x — e.g. Ubuntu 22.04 or Debian Bullseye — are detected as too old and
> the build automatically falls back to `USE_FLINT=0`. Install a newer FLINT from
> source or a backport if you want the accelerated paths on those systems.

**macOS (Homebrew):**

```bash
brew install gmp mpfr readline cmake
# Optional: FLINT (>= 3.0) for fast, rigorous algebraic-extension arithmetic:
brew install flint
# Optional: GMP-ECM for advanced integer factorization:
brew install gmp-ecm
```

LAPACK/BLAS need not be installed on macOS — the build auto-detects Apple's
**Accelerate** framework.

### Building Mathilda

The `makefile` auto-discovers `src/*.c`, configures and compiles internal dependencies, then links the main executable (`-std=c99 -O3`).

1. Clone the repository:
   ```bash
   git clone https://github.com/stblake/Mathilda.git
   cd Mathilda
   ```
   Install GMP-ECM (used for advanced integer factorization) from your package
   manager — `brew install gmp-ecm` on macOS or `sudo apt install libecm-dev`
   on Debian/Ubuntu. The build autodetects it and links `-lecm`; if it is
   absent, the build still succeeds with advanced factorization disabled
   (equivalent to `make USE_ECM=0`).
2. Build the project:
   ```bash
   make -j$(nproc)
   ```
   To build a leaner binary, disable optional backends, e.g.
   `make -j$(nproc) USE_LAPACK=0 USE_MPFR=0`.
3. Start the interactive REPL:
   ```bash
   ./Mathilda
   ```

### Running the Test Suite

Mathilda ships a comprehensive C-based unit suite — **216 `test_*.c` files**
covering evaluation, parsing, pattern matching, arithmetic, polynomials,
calculus, linear algebra, and more. CMake auto-detects the same optional
backends (MPFR, LAPACK, ECM) during configuration.

```bash
cd tests
mkdir -p build && cd build
cmake ..
make -j$(nproc)

# Run all test binaries
for t in *_tests; do ./$t; done
```

---

## 🛠️ Developer Guide & Architecture

Everything in Mathilda is represented by an immutable-by-convention `Expr` AST
node — a tagged union over `Integer`, `Real`, `BigInt`, `Symbol`, `String`,
`Function`, and (when built with MPFR) arbitrary-precision `MPFR` reals. Compound
values are `Function` nodes: `Rational[n, d]`, `Complex[a, b]`, and lists
(`List[…]`) are all expressions, so the same generic tools (`Part`, `Map`,
`ReplaceAll`, …) operate uniformly on everything.

The system is modularized into several independent subsystems:

1. **Parser (`parse.c`)** — a Pratt parser mirroring Mathematica's operator precedences (inline lexing, no separate tokenizer).
2. **Evaluator (`eval.c`)** — the fixed-point evaluation loop; applies `Hold*`, `Flat`, `Orderless`, and `Listable` before recursively evaluating arguments.
3. **Symbol table (`symtab.c`)** — `OwnValues`, `DownValues` (user rules), attributes, docstrings, and native C built-in function pointers.
4. **Pattern matcher & rule engine (`match.c`, `replace.c`)** — structural tree unification with sequence segmenting and backtracking.

Larger mathematical domains live in dedicated subdirectories of `src/`:

| Subsystem | Responsibility |
|-----------|----------------|
| `poly/`     | Univariate/multivariate polynomial arithmetic, factorization, algebraic-number fields, Gröbner bases |
| `linalg/`   | Dense linear algebra; eigen kernels split by algorithm |
| `calculus/` | `D`/`Dt`/`Derivative`, `Series`, `Limit`, `Integrate` (incl. Risch–Norman) |
| `simp/`     | `Simplify`, trigonometric simplification, radical denesting, assumptions |
| `sum/`      | Symbolic summation (polynomial, geometric, Gosper) |
| `internal/` | Mathematica-syntax bootstrap `.m` files (init, integral tables) loaded at startup |

A recurring design pattern is **C for performance, rules for mathematics**: hot
paths (parser, evaluator, matcher, arithmetic) are C, while higher-level
identities (integral tables, etc.) are expressed as `DownValues` in Mathilda's
own language. The full architecture is documented in [`SPEC.md`](SPEC.md), and
extension recipes in [`docs/extending.md`](docs/extending.md).

### Extending Mathilda: Adding a New Built-in Function

Adding new functionality to Mathilda is straightforward:

1. **Write the C implementation.**
   Create your evaluation logic in the appropriate `.c` module (e.g., `core.c`).
   Your function signature must be `Expr* builtin_myfunc(Expr* res)`.

   * **Memory rule (ownership contract):** the builtin **takes ownership** of
     `res`. On success, return a **new** `Expr*` — the **evaluator** frees `res`
     for you, so you must **not** call `expr_free(res)` yourself (doing so causes
     a double-free). If you cannot evaluate the input (e.g. symbolic arguments to
     a purely numeric function), return `NULL` **without freeing `res`**, and the
     evaluator retains ownership, leaving the expression unevaluated. When you
     reuse parts of `res` in your result, NULL them out first so the evaluator's
     cleanup doesn't free them twice.

   ```c
   Expr* builtin_myfunc(Expr* res) {
       if (res->data.function.arg_count != 1) return NULL;  /* leave unevaluated */
       // ... mathematical logic ...
       return expr_new_integer(42);  /* evaluator frees res — do NOT free it here */
   }
   ```

2. **Register the function.**
   In the module's initialization routine (e.g., `core_init()`), register the
   function and assign a documentation string so it is available via `?MyFunc`:
   ```c
   symtab_add_builtin("MyFunc", builtin_myfunc);
   symtab_set_docstring("MyFunc", "MyFunc[x]\n\tComputes the ultimate answer.");
   ```

3. **Assign attributes.**
   If your function threads over lists, operates symmetrically, or holds its
   arguments unevaluated, set the corresponding attributes during initialization:
   ```c
   symtab_get_def("MyFunc")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
   ```

4. **Test and document.**
   * Add test cases to the appropriate suite in `tests/` using the `TEST(...)` macro.
   * Document the function in the matching file under [`docs/spec/builtins/`](docs/spec/builtins/) and add a note to the current week's [`docs/spec/changelog/`](docs/spec/changelog/) entry.

See [`docs/extending.md`](docs/extending.md) for the full recipes (modules,
patterns, internal `.m` rules, operators) and [`CLAUDE.md`](CLAUDE.md) for the
contributor workflow.

---

## 📜 Open Source & License

Mathilda is open-source software licensed under the **GNU General Public License v3.0 (GPLv3)**.

You are heavily encouraged to explore the codebase, submit pull requests, report issues, and expand the CAS with new mathematical algorithms! Please see the `LICENSE` file for more details.
