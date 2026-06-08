# Function

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
body & or Function[body]
    represents a pure function with formal parameters #, #1, #2, ... and ##, ##1, ##2, ... for sequences of arguments.
Function[x, body] or Function[{x1, x2, ...}, body]
    represents a pure function with named formal parameters x or x1, x2, ....
Function[params, body, attrs]
    is a pure function that is treated as having attributes attrs for purposes of evaluation.
    attrs can be a single attribute or a list of attributes; recognised attributes include HoldFirst, HoldRest, HoldAll, HoldAllComplete, Listable, Flat, Orderless, OneIdentity, NumericFunction, SequenceHold, and NHoldRest.
Function[Null, body, attrs]
    represents a function in which the parameters in body are given using # etc.

Parameter binding is lexical: named parameters are substituted into the body before evaluation. Nested Function expressions shadow their own parameters.
By default Function has no Hold attributes; the arguments are evaluated before substitution. Adding HoldAll (or HoldFirst / HoldRest / HoldAllComplete) in the 3-arg form holds arguments in the chosen positions.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `Function` is held as an inert head; the real work is `apply_pure_function` (called from `eval.c` when the evaluator sees `Function[...][args]`). Application is by **lexical substitution**, never by mutating the symbol table — this keeps `Unevaluated`-wrapped references correct and avoids variable capture (nested `Function`s are treated as opaque). Supported shapes: `Function[body]` and `Function[Null, body, attrs]` use the Slot form — `substitute_slots` replaces `#`/`#n` (`Slot`) and `##`/`##n` (`SlotSequence`) with the supplied arguments; `Function[x, body]` and `Function[{x1,…}, body]` use named parameters — `substitute_names` replaces each parameter symbol with its bound argument (missing arguments stay symbolic). The substituted body is then `evaluate`d, and `trap_return` strips a `Return[v]` targeting this `Function` boundary.

**Attributes.** By default `Function` has no Hold attributes, so its call arguments are evaluated before substitution. The 3-arg form `Function[params, body, attrs]` opts into attributes: `pure_function_attributes` maps the attribute spec (a symbol or `List` of symbols) onto the standard `ATTR_*` bits (`HoldAll`/`HoldFirst`/`HoldRest`/`HoldAllComplete`, plus `Listable`/`Flat`/`Orderless`/`OneIdentity`/`NumericFunction`/…), which the evaluator consults when preparing the arguments.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Harold Abelson and Gerald Jay Sussman, *Structure and Interpretation of Computer Programs*, 2nd ed., §1.3.2 (lambda; constructing procedures).
- Source: [`src/purefunc.c`](https://github.com/stblake/mathilda/blob/main/src/purefunc.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= (# + 1 &)[10]
Out[1]= 11
```

```mathematica
In[1]:= Function[x, x^2][5]
Out[1]= 25
```

```mathematica
In[1]:= Function[{x, y}, x + y][3, 4]
Out[1]= 7
```

```mathematica
In[1]:= f = #1 - #2 &; f[10, 3]
Out[1]= 7
```

### Notes

`body &` is the anonymous (pure) function shorthand, where `#`/`#1` denotes the
first argument, `#2` the second, and so on. The named forms `Function[x, body]`
and `Function[{x1, x2, ...}, body]` bind explicit parameter symbols. Parameter
binding is lexical: arguments are substituted into the body before evaluation,
and nested `Function`s shadow their own parameters. A pure function can be stored
in a symbol (`f = #1 - #2 &`) and called like any other function head.
