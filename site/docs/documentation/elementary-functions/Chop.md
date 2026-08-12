# Chop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Chop[expr]`**

replaces approximate real numbers in expr that are close to zero by the exact integer 0.

**`Chop[expr, delta]`**

replaces numbers smaller in absolute magnitude than delta by 0.

<details>
<summary>Notes</summary>

The default tolerance is 10^-10. Chop walks the entire expression tree, so small real-valued subterms inside arbitrary heads, lists, and held forms are all chopped. Exact numbers -- integers, bigints, rationals, and symbolic constants -- pass through untouched. For machine complex numbers Complex\[re, im\] whose real and imaginary parts are both machine reals: if only the imaginary part is below tolerance the whole Complex wrapper is dropped and the real part is returned; if only the real part is below tolerance the result is Complex\[0., im\], preserving the machine-complex shape with a machine zero. If both parts are below tolerance the result is the exact integer 0.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Chop[Exp[N[Range[4] Pi I]]]
Out[1]= {-1.0, 1.0, -1.0, 1.0}

In[2]:= Chop[N[Pi] - Rationalize[N[Pi], 10^-12]] === 0
Out[2]= True

In[3]:= Chop[N[Pi] - Rationalize[N[Pi], 10^-12], 10^-14] === 0
Out[3]= False

In[4]:= Chop[10.^-12 + 2. I]
Out[4]= 0.0 + 2.0*I

In[5]:= Chop[2. + 10.^-12 I]
Out[5]= 2.0
```

### Applications (5)

```mathematica
In[1]:= Chop[N[Sin[Pi]]]
Out[1]= 0
```

```mathematica
In[1]:= Chop[{1.5, 1.0*^-12, 3.2}]
Out[1]= {1.5, 0, 3.2}
```

```mathematica
In[1]:= Chop[N[Exp[I Pi] + 1]]
Out[1]= 0
```

```mathematica
In[1]:= Chop[3.0 + 1.0*^-15 I]
Out[1]= 3.0
```

```mathematica
In[1]:= Chop[1.0*^-15 + 2.5 I]
Out[1]= 0.0 + 2.5*I
```

## Options & behaviour

**Complex handling**.  `Complex[re, im]` whose components are both
machine reals is the "machine complex" case and gets Mathematica's
special treatment:

| `re` below tolerance | `im` below tolerance | result |
|----------------------|----------------------|--------|
| yes                  | yes                  | exact integer `0` |
| no                   | yes                  | `re` (Complex wrapper dropped, machine real survives) |
| yes                  | no                   | `Complex[0., im]` &mdash; the real part is the machine zero `0.0`, not an exact `0`, preserving the machine-complex shape |
| no                   | no                   | unchanged |

When `Complex[re, im]` has at least one exact component (e.g.
`Complex[1, 1.*^-12]`), Chop recurses into each part normally: a tiny
`Real` becomes the exact integer `0`, and `builtin_complex` then
collapses `Complex[r, 0]` to `r` on the next evaluator pass.

## Implementation notes

**Algorithm.** `builtin_chop` reads an optional `delta` (default `1.0e-10`) via
`chop_extract_delta` (accepting Integer, BigInt, Real, MPFR, and `Rational[n,d]`,
returning the absolute value; a non-coercible argument makes the whole call stay
unevaluated). It then recurses with `chop_recursive`. For a "machine real" leaf
(detected by `chop_is_machine_real`), it replaces the value by the exact integer
`0` when `fabs(value) < delta`. `Complex` leaves are handled specially: when both
parts are machine reals it independently chops each, returning `0` if both are
small, the real part alone if only the imaginary part chops, or
`Complex[0., im]` (built with `make_complex` to avoid re-evaluation) if only the
real part chops. Non-machine `Complex` nodes recurse into each part and let
`builtin_complex` collapse `Complex[r, 0]` on the next evaluator pass. Function
nodes are rebuilt with each argument chopped; all other atoms (integers, bigints,
symbols, strings, rational components) are copied unchanged.

**Data structures.** Pure `Expr*` tree rewrite — no auxiliary structure. The
comparison is `double`-based (`chop_to_double`), so only inexact numeric leaves
are ever affected; exact integers and rationals are never chopped.

- `Protected`.
- Walks the entire expression tree, so small real-valued subterms inside
  arbitrary heads, lists, and held forms are all chopped.
- Exact numbers (`Integer`, BigInt, `Rational`, symbolic constants) and
  symbolic input pass through untouched.
- `delta` may be supplied as `Integer`, `Real`, `Rational[n, d]`,
  BigInt, or (when `USE_MPFR=1`) `MPFR`; its absolute value is the
  effective tolerance.
- **Packed/NDArray fast path** (packed-aware): a machine buffer is scanned and
  kept untouched when nothing chops (the common cleanup call); a chopping rank-1
  buffer builds the exact mixed `{…, 0, …}` result in one pass, else it degrades
  to the list path. See [packed arrays](../packed-arrays/index.md).
- **`Compile[]`/auto-compilation**: `Chop[x]` / `Chop[x, δ]` (literal `δ`) lower to
  the branchless `x * (Abs[x] >= δ)`, returning a machine `0.` where the
  interpreter gives the exact `Integer 0`; a machine integer is the identity.

**Attributes:** `Protected`.

## See also

[Rational](../../arithmetic/Rational/)

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_chop.c`](https://github.com/stblake/mathilda/blob/main/tests/test_chop.c)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)

## Notes & additional examples

### Notes

`Chop[expr]` replaces approximate reals within `delta` of zero (default
`10^-10`) by the exact integer `0`, walking the whole expression tree. It is the
standard tool for cleaning numerical noise: `Chop[N[Sin[Pi]]]` and Euler's
identity `Chop[N[Exp[I Pi] + 1]]` both collapse to a clean `0`. For machine
complex numbers it chops parts independently — dropping a negligible imaginary
part recovers a pure real (`3.0`), while a negligible real part leaves a
`Complex[0., im]` that preserves the machine-complex shape. Exact numbers and
symbolic constants pass through untouched.
