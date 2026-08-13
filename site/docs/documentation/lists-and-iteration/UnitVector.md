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

## Examples (3)

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

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_unit_vector.c`](https://github.com/stblake/mathilda/blob/main/tests/test_unit_vector.c)
