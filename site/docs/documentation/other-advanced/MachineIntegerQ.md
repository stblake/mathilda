# MachineIntegerQ

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`MachineIntegerQ[expr]`**

gives True if expr is a machine-word (64-bit) integer, False otherwise.

<details>
<summary>Notes</summary>

Unlike IntegerQ, returns False for a BigInt: MachineIntegerQ\[2^100\] is False because that value has been promoted out of a 64-bit word.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

A 64-bit machine-word integer

```mathematica
In[1]:= MachineIntegerQ[5]
Out[1]= True
```

Too big for a word -- a BigInt, so False

```mathematica
In[2]:= MachineIntegerQ[2^100]
Out[2]= False
```

Overflows a signed 64-bit word, promoted to BigInt

```mathematica
In[3]:= MachineIntegerQ[2^63]
Out[3]= False
```

A real, not an integer

```mathematica
In[4]:= MachineIntegerQ[3.0]
Out[4]= False
```

## Implementation notes

**Algorithm.** `builtin_machineintegerq` (`src/core.c`) is a one-line type test: it takes a
single argument and returns `True` exactly when that argument's tagged-union type is
`EXPR_INTEGER` — Mathilda's 64-bit machine-word integer — and `False` otherwise. It does no
arithmetic and never looks at the value's magnitude, because the representation already
encodes the distinction: a value too large for a machine word has been promoted to
`EXPR_BIGINT` (a GMP `mpz_t`) by the arithmetic heads, so it is no longer `EXPR_INTEGER`.

**Data structures.** None beyond the input `Expr`. The function inspects `res->data.function
.args[0]->type` and constructs a fresh `True`/`False` symbol. It is registered `Protected`
in `core_init`; a non-unary call returns `NULL` and is left unevaluated.

**Complexity / limits.** `O(1)`. This is the sharp difference from `IntegerQ`:
`MachineIntegerQ[2^100]` is `False` (the value is a `BigInt`), whereas `IntegerQ[2^100]` is
`True`. It is likewise `False` for a real such as `3.0`, and `False` for `2^63`, which
overflows a signed 64-bit word and so is stored as a `BigInt`. Use it to decide whether a
value will travel on the packed/fixnum fast paths rather than the arbitrary-precision path.

**Attributes:** `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`MachineIntegerQ[expr]` gives `True` only when `expr` is a 64-bit machine-word integer, and
`False` otherwise. It is an `O(1)` type test — it reads the value's internal representation
and does no arithmetic.

It differs from `IntegerQ` precisely on big integers: `MachineIntegerQ[2^100]` is `False`
because that value has been promoted out of a machine word into a `BigInt`, while
`IntegerQ[2^100]` is `True`. The `2^63` case shows the boundary — it overflows a signed
64-bit word and so is already a `BigInt`. Use it to tell whether a value will ride the
machine-integer fast paths rather than the arbitrary-precision path.
