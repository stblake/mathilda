### Worked examples

```mathematica
In[1]:= x = 5
In[2]:= x *= 4
In[3]:= x
```

### Notes

`x *= dx` (`TimesBy`) is shorthand for `x = x dx`: it multiplies the current
value of `x` in place and returns the **new** value. It shares its whole
machinery with `AddTo`/`SubtractFrom`/`DivideBy`, differing only in that the
combining head is `Times`.

The target must already hold a value; applied to an unassigned symbol the call
is left unevaluated with a `TimesBy::rvalue` message. The left-hand side may be a
part target (`v[[i]] *= c`), and the product is formed and re-evaluated, so
threading over lists and symbolic simplification behave as for `x = x dx`.
