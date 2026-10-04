# Compile

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Compile[{x, ...}, expr] or Compile[{{x, _Real}, ...}, expr] builds a CompiledFunction that evaluates expr over machine numbers (types _Real, _Integer, _Complex; default _Real), falling back to the interpreter for symbolic arguments or non-compilable bodies. With RuntimeAttributes -> Listable the object threads over List arguments; the default is RuntimeAttributes -> {}. RuntimeOptions -> {"CatchMachineIntegerOverflow" -> False} (or the shorthand RuntimeOptions -> "Speed") lets machine-integer arithmetic wrap instead of falling back to the interpreter, which is faster and gives a different answer from the interpreter once a result leaves the machine-integer range; the default True never does. WorkingPrecision -> n compiles real/complex arithmetic in MPFR at n decimal digits (one fixed precision for the whole function), for the straight-line arithmetic + elementary-function subset; MachinePrecision (the default) keeps the machine path unchanged. "BigIntegers" -> True makes integer arithmetic exact (GMP) instead of int64.`**

## Examples (17)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= f = Compile[{{x, _Real}}, x^2 + 1]
Out[1]= CompiledFunction[{x}, x^2 + 1]

In[2]:= f[3.0]
Out[2]= 10.0
```

Symbolic argument -> interpreter fallback

```mathematica
In[3]:= f[a]
Out[3]= 1 + a^2
```

```mathematica
In[4]:= g = Compile[{{n, _Integer}}, Module[{s = 0.}, Do[s = s + 1/i^2, {i, 1, n}]; s]]; g[100]
Out[4]= 1.63498

In[5]:= Compile[{{z, _Complex}}, z^2][1.0 + 2.0 I]
Out[5]= -3.0 + 4.0*I

In[6]:= Compile[{{m, _Real, 2}}, m[[All, 1]]][{{1., 2.}, {3., 4.}}]
Out[6]= {1.0, 3.0}
```

A 5-point stencil: read an argument grid, write a local copy

```mathematica
In[7]:= Compile[{{a, _Real, 2}}, Module[{n = Length[a], b = a}, Do[b[[i, j]] = (a[[i - 1, j]] + a[[i + 1, j]] + a[[i, j - 1]] + a[[i, j + 1]])/4, {i, 2, n - 1}, {j, 2, n - 1}]; b]][Table[1.0 (10 i + j), {i, 1, 3}, {j, 1, 3}]]
Out[7]= {{11.0, 12.0, 13.0}, {21.0, 22.0, 23.0}, {31.0, 32.0, 33.0}}
```

### Scope (1)

```mathematica
In[8]:= Compile[{{n, _Integer}}, n^20, "BigIntegers" -> True][99]
Out[8]= 8179069375972308708891986605443361898001
```

### Options (4)

RuntimeAttributes -> Listable: the object threads over lists

```mathematica
In[9]:= h = Compile[{{x, _Real}}, If[x > 0, 1., -1.], RuntimeAttributes -> Listable]; h[{1., -2., 3.}]
Out[9]= {1.0, -1.0, 1.0}
```

A rank-1 parameter consumes one level, so this maps over the rows

```mathematica
In[10]:= Compile[{{v, _Real, 1}}, Total[v], RuntimeAttributes -> Listable][ {{1., 2.}, {3., 4.}}]
Out[10]= {3.0, 7.0}
```

```mathematica
In[11]:= f = Compile[{{x, _Real}}, Sin[x] Cos[x] + x^3, WorkingPrecision -> 40]; f[N[7/5, 40]]
Out[11]= 2.9114940750779524597719268763562110530148

In[12]:= Compile[{{z, _Complex}}, Exp[z] + z^2, WorkingPrecision -> 45][N[1/2 + I/3, 45]]
Out[12]= 1.696859506173835946311071541529964410492750645 + 0.8727861896014286091867555415942467494395774456*I
```

### Applications (5)

A machine-number function applied to a Real

```mathematica
In[13]:= f = Compile[{{x, _Real}}, x^2 + 1]; f[3.0]
Out[13]= 10.0
```

A symbolic argument transparently falls back to the interpreter

```mathematica
In[14]:= f[a]
Out[14]= 1 + a^2
```

Integer arithmetic stays exact past the int64 range

```mathematica
In[15]:= Compile[{{n, _Integer}}, n^3][3000000]
Out[15]= 27000000000000000000
```

Complex machine arithmetic

```mathematica
In[16]:= Compile[{{z, _Complex}}, z^2][1.0 + 2.0 I]
Out[16]= -3.0 + 4.0*I
```

The object threads over a list

```mathematica
In[17]:= Compile[{{x, _Real}}, If[x > 0, 1., -1.], RuntimeAttributes -> Listable][{1., -2., 3.}]
Out[17]= {1.0, -1.0, 1.0}
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| NDSolve Van der Pol mu=10 | 0.441 s | 0.566 s | 3.32 s |
| NDSolve harmonic oscillator | 0.216 s | 0.214 s | 7.49 s |
| NDSolve long horizon, t to 200 | 0.111 s | 0.152 s | 2.89 s |
| NDSolve y'=-y on [0,10] | 0.065 s | 0.156 s | 2.41 s |
| NDSolve y'=y^2 t nonlinear | 0.034 s | 0.135 s | 0.3 s |

## Implementation notes

**Algorithm.** `Compile` is `HoldAll, Protected`, so the body is compiled in its
raw unevaluated form and the argument symbols stay local (a global value for one
does not leak in). `builtin_compile` seeds `RuntimeAttributes`/`RuntimeOptions`
from the registered defaults (so `SetOptions[Compile, …]` takes effect), then
strips trailing options right-to-left — `RuntimeAttributes`, `RuntimeOptions`,
`WorkingPrecision`, `"BigIntegers"` — the rightmost setting winning; an option it
has no meaning for leaves `Compile[…]` unevaluated rather than being silently
ignored. The remaining two arguments (argspec, body) go to
`compiled_function_new`, which parses the spec and type-infers and lowers the body
to typed bytecode; a malformed argspec makes it return `NULL`, again leaving
`Compile[…]` unevaluated. Success wraps the program in an `EXPR_COMPILED` atom.

**Data structures.** The result is an `EXPR_COMPILED` holding a reference-counted,
immutable `CompiledFunction` (a bare pointer in the `Expr` union, so no node-size
growth). The VM is a register machine with three banks — scalar `R`, array handle
`V`, strip-mining tile `T` — and a reusable per-program frame (no per-call
`malloc`). It shares the machine-precision kernel registry (`ndkernels.c`) with
NDSolve/Plot/NIntegrate, so every special function with a machine kernel lowers
through the same `KERN_*` opcodes.

**Complexity / limits.** Applying the object to numeric arguments runs the
bytecode with no `Expr` allocation; a symbolic argument, or a body outside the
compilable subset, transparently falls back to the interpreter, so the value is
always what the uncompiled body would give. Integer arithmetic is exact — an
`int64` overflow abandons the compiled call and the interpreter promotes to a
bignum — unless `RuntimeOptions -> "Speed"` opts into wrapping. `WorkingPrecision`
and `"BigIntegers"` opt into MPFR / GMP arithmetic; both are off by default, so
the machine path is unchanged when neither is given.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Mod](../../arithmetic/Mod/), [Quotient](../../arithmetic/Quotient/), [Power](../../arithmetic/Power/), [Gamma](../../special-functions/Gamma/), [Erf](../../special-functions/Erf/), [BesselJ](../../special-functions/BesselJ/), [Zeta](../../special-functions/Zeta/)

- Source: [`src/compile/compiled_function.c`](https://github.com/stblake/mathilda/blob/main/src/compile/compiled_function.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_bitwise.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bitwise.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_arbprec.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_arbprec.c)

## Notes & additional examples

### Notes

`Compile[{arg, …}, expr]` builds a `CompiledFunction` that evaluates `expr` over
raw machine numbers, bypassing the symbolic evaluator. Argument types are `_Real`
(the default for a bare symbol), `_Integer` or `_Complex`, with an optional rank
for array parameters (`{v, _Real, 1}`). It is `HoldAll`, so the body compiles in
its raw form and the argument symbols stay local.

A symbolic argument, or a body outside the compilable subset, transparently falls
back to the interpreter (the second example) — so a `CompiledFunction` always
returns what the original expression would. Integer arithmetic is exact: an
`int64` overflow re-runs the body through the interpreter and promotes to a
bignum. Use `CompileDiagnostics` to find out whether, and why, a body does not
lower; `WorkingPrecision -> n` and `"BigIntegers" -> True` opt into MPFR/GMP
arithmetic.
