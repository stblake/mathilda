### Worked examples

`Switch` returns the value for the first matching pattern; `_` is the default:

```mathematica
In[1]:= Switch[5, _Integer, "int", _Real, "real", _, "other"]
Out[1]= "int"
```

```mathematica
In[1]:= Switch[2.5, _Integer, "int", _Real, "real", _, "other"]
Out[1]= "real"
```

Forms can be arbitrary patterns, including `PatternTest` predicates:

```mathematica
In[1]:= Switch[7, _?PrimeQ, "prime", _, "composite"]
Out[1]= "prime"
```

Combined with `Table` and alternative (`|`) patterns it expresses FizzBuzz in a
single dispatch:

```mathematica
In[1]:= Table[Switch[Mod[n, 15], 0, "FizzBuzz", 3 | 6 | 9 | 12, "Fizz", 5 | 10, "Buzz", _, n], {n, 1, 15}]
Out[1]= {1, 2, "Fizz", 4, "Buzz", "Fizz", 7, 8, "Fizz", "Buzz", 11, "Fizz", 13, 14, "FizzBuzz"}
```

If no form matches and there is no default, the call is returned unevaluated:

```mathematica
In[1]:= Switch[x, 1, "a", 2, "b"]
Out[1]= Switch[x, 1, "a", 2, "b"]
```

### Notes

`Switch[expr, form_1, value_1, ...]` evaluates `expr`, compares it against each
`form_i` in turn, and returns the corresponding `value_i` for the first match.
Because of the `HoldRest` attribute, the form/value pairs are held: each form is
evaluated only when its match is tried, and only the chosen value is evaluated.
A trailing `_` form acts as a catch-all default; with no matching form and no
default, the `Switch` stays unevaluated. `Break`, `Return`, and `Throw` inside
the selected value behave as in any other held context.
