# Flat

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

```text
Flat is an attribute that can be assigned to a symbol f to indicate that all expressions involving nested functions f should be flattened out. This property is accounted for in pattern matching.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`Flat` is an attribute marker symbol, not a function. `attr.c` maps the name to the `ATTR_FLAT` bitflag (via `SetAttributes`/the attribute name table). When a head carries `ATTR_FLAT`, the evaluator's flattening step (`eval_flatten_args` in `eval.c`) splices nested same-head calls into the parent argument list, giving associative behaviour (`f[a, f[b, c]]` → `f[a, b, c]`). The pattern matcher (`match.c`) also consults `ATTR_FLAT` when matching sequence patterns against such heads. Plus and Times set this bit. The symbol `Flat` itself only ever appears as an argument to `Attributes`/`SetAttributes` or inside a `Function[..., Flat]` attribute spec (`pure_function_attributes` in `purefunc.c` maps `SYM_Flat` → `ATTR_FLAT`).

**Attributes:** none registered.

## Implementation status

**Partial** — implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## References

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= SetAttributes[f, Flat]
Out[1]= Null

In[2]:= f[f[a, b], c]
Out[2]= f[a, b, c]

In[3]:= f[a, f[b, f[c, d]]]
Out[3]= f[a, b, c, d]
```

### Notes

`Flat` marks a head as associative, so nested calls with the same head are flattened into a single argument sequence. `Plus` and `Times` carry this attribute by default; the flattening is also accounted for during pattern matching.
