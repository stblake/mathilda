### Worked examples

```mathematica
In[1]:= IntegerQ[5]  (* an exact integer *)
```

```mathematica
In[1]:= IntegerQ[5.0]  (* a Real is not an integer, even at an integral value *)
```

```mathematica
In[1]:= IntegerQ[1/2]  (* a Rational is not an integer *)
```

```mathematica
In[1]:= IntegerQ[x]  (* a symbol is not known to be one, so False *)
```

### Notes

`IntegerQ[e]` is `True` exactly when `e` is an exact integer — a machine `Integer` or an
arbitrary-precision bigint — and `False` for everything else, including reals at integral
values, rationals, and symbols. Like the other `*Q` predicates it always returns a
boolean; it never stays unevaluated.

Use it to test exactness: `IntegerQ[5.0]` is `False` because `5.0` is a machine real, not
an integer.
