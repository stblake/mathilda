# Attributes

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Attributes[s] gives the list of attributes for s.`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Attributes[Plus]
Out[1]= {Flat, Listable, NumericFunction, OneIdentity, Orderless, Protected}
```

## Implementation notes

`builtin_attributes` (`src/attr.c`) reads the symbol's attribute bitflags via `get_attributes(name)` and builds a `List` of attribute symbols from them, e.g. `ATTR_FLAT` -> `Flat`, the `HoldFirst|HoldRest` pair collapsing to `HoldAll`, `ATTR_LISTABLE` -> `Listable`, `ATTR_PROTECTED` -> `Protected`, etc. The symbol itself is held unevaluated (`Attributes` carries `ATTR_HOLDALL`).

- Common attributes include `Flat` (associativity), `Orderless` (commutativity), `Listable` (automatic threading over lists), `HoldFirst`, `HoldRest`, `HoldAll` (evaluation control), and `Protected`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [Orderless](../../expression-information/Orderless/), [HoldFirst](../../other-advanced/HoldFirst/), [HoldRest](../../other-advanced/HoldRest/), [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_array_flatten.c`](https://github.com/stblake/mathilda/blob/main/tests/test_array_flatten.c)
- Tests: [`tests/test_chop.c`](https://github.com/stblake/mathilda/blob/main/tests/test_chop.c)
- Tests: [`tests/test_clearall_remove_protect.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clearall_remove_protect.c)
- Tests: [`tests/test_clip.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clip.c)
