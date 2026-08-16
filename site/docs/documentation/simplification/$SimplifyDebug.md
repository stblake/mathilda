# $SimplifyDebug

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$SimplifyDebug`**

When set to True, Simplify prints one stderr line per transform invocation: /Name/: \<input\> -\> \<output\> \[\<ms\> ms\]. Defaults to False. Useful for diagnosing slow Simplify calls.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= $SimplifyDebug = True; Simplify[a x + b x]; $SimplifyDebug = False;
```

### Applications (1)

```mathematica
In[2]:= $SimplifyDebug
Out[2]= False
```

## Implementation notes

A diagnostic flag (not a builtin), given an OwnValue defaulting to `False` in
`simp_init`. When set to `True`, `simp_debug_enabled` (read directly off the
OwnValue list to avoid re-evaluation) causes `traced_call_unary` /
`simp_debug_log` to emit one stderr line per transform invocation inside the
Simplify search, in the form `/<TransformName>/: <input> -> <output> [<ms> ms]`,
used to diagnose slow Simplify calls and runaway candidate explosion.

- Default `False`. When set to `True`, `Simplify` prints one line per transform
  invocation to **stderr**, in the form
  `/<TransformName>/: <input> -> <output> [<elapsed> ms]`. Useful for diagnosing
  slow or hanging `Simplify` calls and runaway candidate-set growth. The value is
  read directly off the `OwnValue`, so there is no cost when it is `False`.

**Attributes:** none registered.

## References

**See also:** [Simplify](../../simplification/Simplify/)

- Source: [`src/simp/simp_util.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp_util.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)

## Notes & additional examples

### Notes

`$SimplifyDebug` is a global flag, default `False`. When set to `True`,
`Simplify` prints one stderr line per transform invocation
(`/Name/: <input> -> <output> [<ms> ms]`), which is useful for diagnosing slow
`Simplify` calls.
