# Mathilda

**A small, highly efficient computer algebra system in written in C99.**

[![Build](https://github.com/stblake/mathilda/actions/workflows/build.yml/badge.svg)](https://github.com/stblake/mathilda/actions/workflows/build.yml)
[![Tests](https://github.com/stblake/mathilda/actions/workflows/tests.yml/badge.svg)](https://github.com/stblake/mathilda/actions/workflows/tests.yml)
[![License: GPL v3](https://img.shields.io/badge/license-GPLv3-blue.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-stblake.github.io%2Fmathilda-informational)](https://stblake.github.io/mathilda/)
[![Downloads](https://img.shields.io/endpoint?url=https%3A%2F%2Fraw.githubusercontent.com%2Fstblake%2Fmathilda%2Ftraffic-stats%2Fclones-badge.json&cacheSeconds=3600)](.github/workflows/traffic.yml)

Mathilda is built around a small, general core: a recursive expression model,
infinite (fixed-point) evaluation driven by attributes, structural pattern
matching with backtracking, and rewriting rules. On top of that core it builds a
broad mathematical library: exact and arbitrary-precision arithmetic, symbolic
calculus, polynomial algebra, equation solving, special functions, numerical
analysis, and a `Compile[]` bytecode compiler over transparently packed arrays.

The kernel is about 580,000 lines of C99 across 850 source files, plus 14,000
lines of Mathilda's own language for the rule-based parts of the library. It
exposes more than 1,250 built-in symbols and is checked by over 600 C unit-test
suites. Big-number arithmetic uses GMP, MPFR and (optionally) FLINT. 

## A taste of Mathilda

```
In[1]:= Integrate[1/(x^3 + 1), x]
Out[1]= ArcTan[(-1 + 2 x)/Sqrt[3]]/Sqrt[3] + 1/3 Log[1 + x] - 1/6 Log[1 - x + x^2]

In[2]:= Series[Tan[x], {x, 0, 7}]
Out[2]= x + 1/3 x^3 + 2/15 x^5 + 17/315 x^7 + O[x]^8

In[3]:= DSolve[y''[x] + y[x] == Sin[x], y[x], x]
Out[3]= {{y[x] -> -1/2 x Cos[x] + C[1] Cos[x] - C[2] Sin[x]}}

In[4]:= Solve[x^2 - 2 y^2 == 1 && 0 < x < 100 && y > 0, {x, y}, Integers]
Out[4]= {{x -> 3, y -> 2}, {x -> 17, y -> 12}, {x -> 99, y -> 70}}

In[5]:= Sum[1/k^2, {k, 1, Infinity}]
Out[5]= 1/6 Pi^2

In[6]:= FactorInteger[2^128 + 1]
Out[6]= {{59649589127497217, 1}, {5704689200685129054721, 1}}

In[7]:= fib[0] = 0; fib[1] = 1; fib[n_] := fib[n] = fib[n - 1] + fib[n - 2]; fib[200]
Out[7]= 280571172992510140037611932413038677189525

In[8]:= N[Pi, 50]
Out[8]= 3.1415926535897932384626433832795028841971693993751
```

## Features

- **Evaluation engine.** Fixed-point evaluation; attribute-driven argument
  handling (`Hold*`, `Flat`, `Orderless`, `Listable`, `OneIdentity`,
  `Protected`, …); `OwnValues`, `DownValues` and `UpValues`; `Module`, `Block`,
  `With`; and non-local control flow (`Catch`/`Throw`, `Sow`/`Reap`).
- **Patterns and rules.** `_`, `__`, `___`, `x_h`, `/;`, `?`, `:`, `..`,
  `Alternatives`, `Except` and `HoldPattern`, with full sequence backtracking.
  Rules (`->`, `:>`) and replacement (`/.`, `//.`, `Replace`).
- **Numbers.** Arbitrary-precision integers and exact rationals,
  complex numbers, and MPFR reals with precision and accuracy tracking
  (`` 3.98`50 ``). Also `Interval` arithmetic.
- **Calculus.** `D`/`Dt`; `Integrate`, which cascades through Risch–Norman,
  Cherry/Liouvillian, algebraic (Trager), Weierstrass, CRC-table and other
  methods, plus definite and contour integration. Also `Series`, `Limit` (a
  Gruntz engine), `Sum`/`Product` (Gosper, hypergeometric and zeta families),
  and `DSolve`.
- **Algebra.** Univariate and multivariate factorization (over ℚ and number
  fields ℚ(α)), `Together`/`Apart`/`Cancel`, Gröbner bases, resultants,
  `Simplify`/`FullSimplify` with assumptions, and `Root`/`RootReduce`.
- **Equation solving.** `Solve`, `Reduce` (cylindrical algebraic decomposition
  and quantifier elimination) and `FindInstance`. `Solve[…, Integers]` handles
  Diophantine equations: Hermite normal form, Pell and Thue equations, ternary
  forms and more.
- **Numerics.** `NIntegrate`, `NSum`, `NDSolve` (ODEs and PDEs), `NMinimize`,
  `FindMinimum` (L-BFGS-B, SLSQP, trust-region, …), `FindRoot`, `NSolve`.
- **High-performance arrays.** Lists of machine numbers are transparently packed
  into dense buffers, and element-wise, reduction and structural kernels
  (optionally multithreaded) run on the buffer directly. `Compile[]` lowers
  numeric code to a register VM and auto-compiles `Table`/`Sum`/`NIntegrate`
  bodies.
- **Linear algebra.** LU, QR, Cholesky, SVD and Schur decompositions, and
  eigensystems, with LAPACK/BLAS acceleration. Also `SparseArray`.
- **Special functions.** Gamma, Zeta, PolyLog, the Bessel and Airy families,
  hypergeometric `pFq`, `MeijerG`, elliptic integrals and orthogonal
  polynomials. FLINT/arb gives rigorous numerics.
- **Number theory.** Factorization (Pollard ρ and *p*−1, Williams *p*+1, CFRAC,
  ECM, …), modular arithmetic, discrete logarithms and Dirichlet characters.
- **And more.** Associations with O(1) lookup, strings and regular expressions,
  graphs, statistics and machine learning, Fourier transforms, image
  processing, and interactive 2D/3D graphics.

See the [function reference](Mathilda_spec.md) for the complete list, organised
by category.

## Quick start

Mathilda needs GCC (Apple clang is rejected), GMP and MPFR. GNU Readline is
recommended for the REPL. FLINT, GMP-ECM, LAPACK, SuiteSparse (UMFPACK,
CHOLMOD), ARPACK, FFTW, PCRE2 and Raylib are optional and auto-detected.

```bash
# macOS
brew install gcc gmp mpfr readline flint gmp-ecm suite-sparse arpack
# Debian / Ubuntu
sudo apt install build-essential libgmp-dev libmpfr-dev libreadline-dev libflint-dev libecm-dev \
                 libsuitesparse-dev libarpack2-dev

git clone https://github.com/stblake/mathilda.git
cd mathilda
make -j$(nproc)
./Mathilda                  # interactive REPL
./Mathilda -file script.m   # run a script
```

[`docs/building.md`](docs/building.md) covers the rest:

- Windows (MSYS2)
- every build flag
- `make install` and relocatable deployment
- embedding the kernel through `libmathilda` and its C API
- the desktop app

## Documentation

| | |
|---|---|
| [Documentation site](https://stblake.github.io/mathilda/) | Tutorials and the per-function reference |
| [`Mathilda_spec.md`](Mathilda_spec.md) | Built-in function reference, by category |
| [`SPEC.md`](SPEC.md) | System architecture: parser, evaluator, pattern matcher, memory model |
| [`docs/extending.md`](docs/extending.md) | Adding builtins, modules, rules and operators |
| [`docs/design/`](docs/design/) | Design notes (packed arrays, `Compile[]`, message routing, …) |
| [`docs/design/performance.md`](docs/design/performance.md) | Performance benchmarks, method and known gaps |
| [`docs/spec/changelog/`](docs/spec/changelog/) | Weekly changelog |
| [`book/`](book/) | *The Mathilda Book* (LaTeX source) |

## Testing

```bash
cd tests && mkdir -p build && cd build
cmake .. && make -j$(nproc)
ctest -LE nightly -j4        # unit suites (what CI runs on every push)
ctest -j4 --timeout 3600     # everything, including the nightly corpora
```

The source tree is also checked by audits, such as `make check-c99`,
`check-messages`, `check-packed-aware` and `check-compile-coverage`.
[`SPEC.md` §9](SPEC.md#9-build--test) describes each one.

## Contributing

Contributions are welcome. See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the
build and test workflow, the builtin ownership contract, and the conventions
around documentation and versioning.

## License

Mathilda is free software, released under the
[GNU General Public License v3.0](LICENSE).
