# Orderless

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

```text
Orderless is an attribute that can be assigned to a symbol f to indicate that the elements e_i in expressions of the form f[e_1, e_2, ...] should automatically be sorted into canonical order. This property is accounted for in pattern matching.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`Orderless` is an attribute marker symbol, not a function. `attr.c` maps the name to the `ATTR_ORDERLESS` bitflag. A head carrying this bit makes the evaluator's ordering step (`eval_sort_args` in `eval.c`) sort the arguments into canonical order (`expr_compare`), giving commutative behaviour. The pattern matcher also accounts for `ATTR_ORDERLESS` so a pattern can match commuted arguments. Plus and Times set this bit. The symbol appears only as an argument to `Attributes`/`SetAttributes`, or inside a `Function[..., Orderless]` attribute spec (`purefunc.c` maps `SYM_Orderless` → `ATTR_ORDERLESS`).

**Attributes:** none registered.

## Implementation status

**Partial** — implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## References

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= SetAttributes[g, Orderless]
Out[1]= Null

In[2]:= g[3, 1, 2]
Out[2]= g[1, 2, 3]

In[3]:= g[c, a, b]
Out[3]= g[a, b, c]
```

### Notes

`Orderless` marks a head as commutative, so its arguments are automatically sorted into canonical order. `Plus` and `Times` carry this attribute by default; canonical ordering is also used when matching patterns against `Orderless` heads.
