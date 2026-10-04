### Worked examples

```mathematica
In[1]:= x = {3, 1, 2};
```

```mathematica
In[1]:= ApplyTo[x, Sort]  (* rewrites x to f[x] and returns the new value *)
```

```mathematica
In[1]:= x
```

### Notes

`ApplyTo[x, f]` sets `x = f[x]` in place and returns the new value — the function
analogue of `AddTo`/`TimesBy`. It is `HoldFirst`, so `x` is not evaluated before
the update, and the left-hand side may be a plain symbol, a part `x[[i]]`, or an
association entry `x[key]`; the write-back goes through `Set`, which knows all
three shapes. The variable must already have a value, otherwise the update is
declined with an `ApplyTo::rvalue` message.
