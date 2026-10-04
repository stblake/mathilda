### Worked examples

```mathematica
In[1]:= MachineIntegerQ[5]  (* a 64-bit machine-word integer *)
```

```mathematica
In[1]:= MachineIntegerQ[2^100]  (* too big for a word -- a BigInt, so False *)
```

```mathematica
In[1]:= MachineIntegerQ[2^63]  (* overflows a signed 64-bit word, promoted to BigInt *)
```

```mathematica
In[1]:= MachineIntegerQ[3.0]  (* a real, not an integer *)
```

### Notes

`MachineIntegerQ[expr]` gives `True` only when `expr` is a 64-bit machine-word integer, and
`False` otherwise. It is an `O(1)` type test — it reads the value's internal representation
and does no arithmetic.

It differs from `IntegerQ` precisely on big integers: `MachineIntegerQ[2^100]` is `False`
because that value has been promoted out of a machine word into a `BigInt`, while
`IntegerQ[2^100]` is `True`. The `2^63` case shows the boundary — it overflows a signed
64-bit word and so is already a `BigInt`. Use it to tell whether a value will ride the
machine-integer fast paths rather than the arbitrary-precision path.
