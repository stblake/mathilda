# Minus

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Minus[x] is the arithmetic negation of x, equivalent to -x (Times[-1, x]).`**

Listable; Minus\[x, y\] (any count other than one) is left unevaluated with Minus::argx.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Minus[3]
Out[1]= -3

In[2]:= Minus[{1, 2.5, 1/2}]
Out[2]= {-1, -2.5, -1/2}

In[3]:= SortBy[{3, 1, 2}, Minus]
Out[3]= {3, 2, 1}
```

## Algorithm

minus.c -- Minus[x], the functional form of unary negation.

The parser already reads `-x` as Times[-1, x], so Minus only has to exist as a callable head: Minus[x] rewrites to Times[-1, x] and lets Times do the arithmetic (numbers, Rationals, Complex, Plus distribution, packed arrays). This is what makes SortBy[list, Minus] and KeySortBy[a, Minus] work.

```text
Attributes match Mathematica: Listable, NumericFunction, Protected.  Any
```

other argument count emits Minus::argx and stays unevaluated (WL prints Minus[x, y] as x − y but does not evaluate it).

## Implementation notes

- `Listable`, `NumericFunction`, `Protected`, as in Mathematica.
- Rewrites to `Times[-1, x]`, so numbers, rationals, complex numbers, `Plus`
  distribution and symbolic arguments all behave exactly like `-x`.
- Usable as a sort key: `SortBy[list, Minus]`, `KeySortBy[a, Minus]`.
- Packed arrays and visible `NDArray`s stay on the buffer (`Minus` is on the
  packed-aware list and hands the buffer to `Times`); `Compile[]` lowers it at
  scalar and array shapes.
- Any other argument count emits `Minus::argx` and stays unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Plus](../../arithmetic/Plus/), [NDArray](../../linear-algebra/NDArray/), [Times](../../arithmetic/Times/)

- Source: [`src/minus.c`](https://github.com/stblake/mathilda/blob/main/src/minus.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
