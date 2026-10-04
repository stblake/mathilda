# UnitVector

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`UnitVector[k]`**

gives the 2-D unit vector in the k-th direction.

**`UnitVector[n, k]`**

gives the n-D unit vector: a length-n list with a 1 in position k and 0s elsewhere.

<details>
<summary>Notes</summary>

Components are exact integers unless WorkingPrecision requests MachinePrecision or a digit count.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= UnitVector[1]
Out[1]= {1, 0}

In[2]:= UnitVector[3, 2]
Out[2]= {0, 1, 0}
```

### Options (1)

```mathematica
In[3]:= UnitVector[2, WorkingPrecision -> MachinePrecision]
Out[3]= {0.0, 1.0}
```

### Applications (4)

A 1 in position 2 of a length-4 vector

```mathematica
In[4]:= UnitVector[4, 2]
Out[4]= {0, 1, 0, 0}
```

One argument means dimension 2

```mathematica
In[5]:= UnitVector[2]
Out[5]= {0, 1}
```

The last basis vector

```mathematica
In[6]:= UnitVector[5, 5]
Out[6]= {0, 0, 0, 0, 1}
```

Exactly one nonzero component

```mathematica
In[7]:= Total[UnitVector[100, 7]]
Out[7]= 1
```

## Algorithm

UnitVector — the n-dimensional unit vector in the k-th direction.

```text
  UnitVector[k]       the 2-D unit vector in the k-th direction
                      (equivalent to UnitVector[2, k]).
  UnitVector[n, k]    the n-D unit vector: a length-n list with a 1 in
                      position k and 0s in every other position.
```

Components are exact integers by default (WorkingPrecision -> Infinity). The WorkingPrecision option selects the component representation, mirroring HilbertMatrix (src/linalg/hilbertmat.c):

```text
  WorkingPrecision -> Infinity          exact integers (default)
  WorkingPrecision -> MachinePrecision  machine-precision Reals
  WorkingPrecision -> d                  d-digit MPFR Reals (d above machine
                                         precision; otherwise machine Reals)
```

Diagnostics mirror Wolfram's surface text:

```text
  - zero arguments               -> UnitVector::argt  (1 or 2 expected)
  - non-option trailing argument -> UnitVector::nonopt
```

Non-integer or out-of-range (k < 1 or k > n) arguments leave the call unevaluated (return NULL), matching Mathilda's "can't evaluate" convention.

## Implementation notes

**Algorithm.** `builtin_unit_vector` builds the length-`n` vector with a `1` in
position `k` and `0` elsewhere. `UnitVector[k]` is the 2-D unit vector
(`UnitVector[2, k]`); `UnitVector[n, k]` the general form. It counts the leading
required (non-rule) arguments, reports a trailing non-option via
`UnitVector::nonopt` and a zero-argument call via `UnitVector::argt`, and
requires `1 <= k <= n` with both positive machine integers — out-of-range or
non-integer arguments leave the call unevaluated (symbolic arguments flow
through).

**Component precision.** The `WorkingPrecision` option selects the component
representation, mirroring `HilbertMatrix`: `Infinity` (default) gives exact
integer components, `MachinePrecision` gives machine `Real`s, and a digit count
`d` above machine precision gives `d`-digit MPFR reals (degrading to machine
precision when `USE_MPFR=0`, with one `UnitVector::wprec` warning). Last valid
setting wins; an unparseable value is ignored.

**Buffer fast path.** For the exact-integer and machine-real modes every
component shares one dtype, so `ndbuild_open` opens a packed `NDArray`
(`NDT_INT64` / `NDT_FLOAT64`), writes zeros, and sets the single `1` — making
`UnitVector[10^6, 1]` a buffer fill rather than 10⁶ boxed nodes for one nonzero
(130 ms → the NumPy-grade path). MPFR components have no buffer form and keep the
`List`. `ndbuild_open` declining (packing off, small `n`) falls through to the
boxed `List`. `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## References

- Source: [`src/vectors.c`](https://github.com/stblake/mathilda/blob/main/src/vectors.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_unit_vector.c`](https://github.com/stblake/mathilda/blob/main/tests/test_unit_vector.c)

## Notes & additional examples

### Notes

`UnitVector[n, k]` is the length-`n` vector with a `1` in position `k` and `0`
everywhere else; `UnitVector[k]` is the two-dimensional case (`UnitVector[2, k]`).
Both `n` and `k` must be positive integers with `1 <= k <= n`.

Components are exact integers by default. `WorkingPrecision -> MachinePrecision`
gives machine reals, and a digit count gives arbitrary-precision reals. At scale
the integer and machine-real forms are returned as a packed array, so
constructing a very long unit vector costs a buffer fill rather than one boxed
element per zero.
