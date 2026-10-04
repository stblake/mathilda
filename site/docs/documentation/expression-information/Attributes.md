# Attributes

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Attributes[s] gives the list of attributes for s.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Attributes[Plus]
Out[1]= {Flat, Listable, NumericFunction, OneIdentity, Orderless, Protected}
```

### Applications (3)

The full bundle an arithmetic head carries

```mathematica
In[2]:= Attributes[Plus]
Out[2]= {Flat, Listable, NumericFunction, OneIdentity, Orderless, Protected}
```

A listable numeric function

```mathematica
In[3]:= Attributes[Sin]
Out[3]= {Listable, NumericFunction, Protected}
```

HoldAll collapses the HoldFirst and HoldRest bits into one token

```mathematica
In[4]:= Attributes[Hold]
Out[4]= {HoldAll, Protected}
```

## Implementation notes

`builtin_attributes` (`src/attr.c`) reads the symbol's attribute bitflags via `get_attributes(name)` and builds a `List` of attribute symbols from them, e.g. `ATTR_FLAT` -> `Flat`, the `HoldFirst|HoldRest` pair collapsing to `HoldAll`, `ATTR_LISTABLE` -> `Listable`, `ATTR_PROTECTED` -> `Protected`, etc. The symbol itself is held unevaluated (`Attributes` carries `ATTR_HOLDALL`).

- Common attributes include `Flat` (associativity), `Orderless` (commutativity), `Listable` (automatic threading over lists), `HoldFirst`, `HoldRest`, `HoldAll` (evaluation control), and `Protected`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [Orderless](../../expression-information/Orderless/), [HoldFirst](../../other-advanced/HoldFirst/), [HoldRest](../../other-advanced/HoldRest/), [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_algebraicnumberdenominator.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumberdenominator.c)
- Tests: [`tests/test_algebraicnumbernorm.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumbernorm.c)
- Tests: [`tests/test_algebraicnumberpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumberpolynomial.c)
- Tests: [`tests/test_algebraicnumbertrace.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumbertrace.c)

## Notes & additional examples

### Notes

`Attributes[sym]` returns the sorted list of attribute symbols set on `sym`, read from its
attribute bitflags. The `HoldFirst | HoldRest` pair is reported as the single token
`HoldAll`; `Protected` marks every builtin.

`Attributes` holds its argument (`HoldAll`), so the symbol is not evaluated first — you
get the attributes of `sym` itself, not of its value. Set and clear them with
`SetAttributes` and `ClearAttributes`.
