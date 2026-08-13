# Compile

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Compile[{x, ...}, expr] or Compile[{{x, _Real}, ...}, expr] builds a CompiledFunction that evaluates expr over machine numbers (types _Real, _Integer, _Complex; default _Real), falling back to the interpreter for symbolic arguments or non-compilable bodies. With RuntimeAttributes -> Listable the object threads over List arguments; the default is RuntimeAttributes -> {}. RuntimeOptions -> {"CatchMachineIntegerOverflow" -> False} (or the shorthand RuntimeOptions -> "Speed") lets machine-integer arithmetic wrap instead of falling back to the interpreter, which is faster and gives a different answer from the interpreter once a result leaves the machine-integer range; the default True never does. WorkingPrecision -> n compiles real/complex arithmetic in MPFR at n decimal digits (one fixed precision for the whole function), for the straight-line arithmetic + elementary-function subset; MachinePrecision (the default) keeps the machine path unchanged. "BigIntegers" -> True makes integer arithmetic exact (GMP) instead of int64.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

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

### Options (2)

RuntimeAttributes -> Listable: the object threads over lists

```mathematica
In[7]:= h = Compile[{{x, _Real}}, If[x > 0, 1., -1.], RuntimeAttributes -> Listable]; h[{1., -2., 3.}]
Out[7]= {1.0, -1.0, 1.0}
```

A rank-1 parameter consumes one level, so this maps over the rows

```mathematica
In[8]:= Compile[{{v, _Real, 1}}, Total[v], RuntimeAttributes -> Listable][ {{1., 2.}, {3., 4.}}]
Out[8]= {3.0, 7.0}
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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Power](../../arithmetic/Power/), [Gamma](../../special-functions/Gamma/), [Erf](../../special-functions/Erf/), [BesselJ](../../special-functions/BesselJ/), [Zeta](../../special-functions/Zeta/), [If](../../control-flow/If/), [Sum](../../calculus/Sum/)

- Source: [`src/compile/compiled_function.c`](https://github.com/stblake/mathilda/blob/main/src/compile/compiled_function.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_arbprec.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_arbprec.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
