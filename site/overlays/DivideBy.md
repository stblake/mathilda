### Worked examples

```mathematica
In[1]:= x = 12
In[2]:= x /= 3
In[3]:= x
```

### Notes

`x /= dx` (`DivideBy`) is shorthand for `x = x/dx`: it divides the current value
of `x` in place and returns the **new** value. The target must already be a
variable with a value — applying it to an unassigned symbol leaves the
expression unevaluated and issues `DivideBy::rvalue`.

Because the update runs through `Set`, the left-hand side may be a compound
target such as `v[[i]] /= c`, and because the new value is formed by evaluating
`Times[old, dx^-1]`, list threading and symbolic simplification apply exactly as
they would to the written-out quotient.
