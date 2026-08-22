# I

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`I`**

is the imaginary unit Sqrt\[-1\].

<details>
<summary>Notes</summary>

I represents the imaginary unit; I^2 evaluates to -1 and complex numbers are written a + b I. It has attribute Protected, and N\[I\] is 0. + 1. I.

</details>

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attribute `Protected`. `Attributes[I] = {Protected}`; the symbol cannot be
  reassigned.
- Carries the OwnValue `Complex[0, 1]`, so `I` evaluates to the imaginary unit;
  `I^2 = -1` and complex numbers are written `a + b I`.
- `N[I] = 0. + 1. I`.

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)
