# ApplyTo

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ApplyTo[x, f]`**

Sets x to f\[x\] and returns the new value. x may be a symbol with a value, a part s\[\[i\]\], or an association entry s\[key\].

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= x = 5; ApplyTo[x, #^2 &]; x
Out[1]= 25

In[2]:= r = <|"n" -> 1|>; ApplyTo[r["n"], # + 10 &]; r
Out[2]= <|"n" -> 11|>
```

## Implementation notes

**Attributes:** `HoldFirst`, `Protected`.

## References

**See also:** [Set](../../assignment-and-rules/Set/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
