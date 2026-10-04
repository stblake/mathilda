# ApplyTo

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ApplyTo[x, f]`**

Sets x to f\[x\] and returns the new value. x may be a symbol with a value, a part s\[\[i\]\], or an association entry s\[key\].

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= x = 5; ApplyTo[x, #^2 &]; x
Out[1]= 25

In[2]:= r = <|"n" -> 1|>; ApplyTo[r["n"], # + 10 &]; r
Out[2]= <|"n" -> 11|>
```

### Applications (3)

```mathematica
In[3]:= x = {3, 1, 2};
```

Rewrites x to f[x] and returns the new value

```mathematica
In[4]:= ApplyTo[x, Sort]
Out[4]= {1, 2, 3}
```

```mathematica
In[5]:= x
Out[5]= {1, 2, 3}
```

## Implementation notes

**Algorithm.** `builtin_applyto` is the in-place update `x = f[x]`, returning the
new value. It is `HoldFirst`, so the first argument arrives as an unevaluated
l-value. `lvalue_root` walks it to the underlying symbol, accepting a bare symbol
`s`, a part `s[[i]]`, or an association entry `s[key]`; the builtin then
evaluates the l-value, applies `f`, and writes the result back by evaluating
`Set[lhs, newvalue]` — so the write-back works for all three l-value shapes
through the one `Set` path.

**Data structures.** Plain `Expr` trees; the current value is adopted into the
`f[...]` call, and the write-back is an ordinary `Set` evaluation. No hashing of
its own.

**Complexity / limits.** The cost is that of evaluating `f[x]` plus the
assignment. The l-value must already have a value — otherwise
`ApplyTo::rvalue` is emitted (via the message funnel) and nothing is changed;
a wrong argument count raises `ApplyTo::argrx`.

**Attributes:** `HoldFirst`, `Protected`.

## References

**See also:** [Set](../../assignment-and-rules/Set/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)

## Notes & additional examples

### Notes

`ApplyTo[x, f]` sets `x = f[x]` in place and returns the new value — the function
analogue of `AddTo`/`TimesBy`. It is `HoldFirst`, so `x` is not evaluated before
the update, and the left-hand side may be a plain symbol, a part `x[[i]]`, or an
association entry `x[key]`; the write-back goes through `Set`, which knows all
three shapes. The variable must already have a value, otherwise the update is
declined with an `ApplyTo::rvalue` message.
