# Rationalize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Rationalize[x]
    converts an approximate number x to a nearby rational with small denominator.
Rationalize[x, dx]
    yields the rational number with smallest denominator that lies within dx of x.

Rationalize[x] yields x unchanged if there is no rational number close enough to x to satisfy |p/q - x| < c/q^2, with c = 10^-4.
Rationalize[x, dx] works with exact numbers x: the value is first numericalised, then rationalised.
Rationalize[x, 0] forces conversion of any inexact number x to rational form, using a tolerance derived from the precision of x.
Rationalize threads over compound expressions and Complex[re, im], so e.g. Rationalize[1.2 + 6.7 x] gives 6/5 + (67 x)/10.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
