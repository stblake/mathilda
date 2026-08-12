# Unique

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Unique[] generates a new symbol; Unique["x"] or Unique[x] uses a name prefix; Unique[{x, ...}] gives a list of fresh symbols. Each is Temporary and never previously used.`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= {Unique[], Unique["x"], Unique[{a, b}]}
Out[1]= {$1, x2, {a3, b3}}
```

## Implementation notes

- `Protected`.
- The numeric suffix is drawn from the shared `$ModuleNumber` counter (the same
  source `Module` uses), advanced until the generated name is unused, so every
  result is genuinely fresh and distinct.
- Created symbols have the `Temporary` attribute.

**Attributes:** `Protected`.

## See also

[Module](../../scoping-constructs/Module/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/scoping-constructs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/scoping-constructs.md)
