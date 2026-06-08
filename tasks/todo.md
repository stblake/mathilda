# Task: Implement EulerGamma as a first-class constant

## Background (already present in tree)
- [x] `N[EulerGamma]` / `N[EulerGamma, n]` — numeric.c constant table + MPFR fill
- [x] `D[EulerGamma, x] = 0` — deriv.c constant list
- [x] `NumericQ[EulerGamma] = True` — core.c is_numeric_quantity
- [x] `SYM_EulerGamma` interned — sym_names.c / sym_names.h
- [x] LaTeX printing `\gamma` — print.c

## To do
- [ ] Add `ATTR_CONSTANT` flag to attr.h (next free bit, 1<<14)
- [ ] Teach attr.c about "Constant": parse in string_to_attribute,
      count+emit in builtin_attributes (alphabetically first)
- [ ] New module src/eulergamma.{c,h}: `eulergamma_init()` registers the
      EulerGamma symbol with attributes Constant | Protected
- [ ] Wire eulergamma_init() into core_init() (core.c) + include header
- [ ] Docstring for EulerGamma in info.c (terse, no examples)
- [ ] Docs: add EulerGamma section to docs/spec/builtins/special-functions.md
- [ ] Changelog note in docs/spec/changelog/2026-06-08.md
- [ ] Build + verify: Attributes, ?EulerGamma, N, D, NumericQ

## Review
All items done. EulerGamma is now a first-class constant symbol.

Files changed:
- src/attr.h         — new ATTR_CONSTANT flag (1<<14)
- src/attr.c         — parse "Constant"; count + emit it (alphabetically first)
- src/eulergamma.{c,h} — NEW module; eulergamma_init() sets {Constant, Protected}
- src/core.c         — include + call eulergamma_init()
- src/info.c         — terse EulerGamma docstring
- tests/CMakeLists.txt — add ../src/eulergamma.c to COMMON_SRC
- docs/spec/builtins/special-functions.md — EulerGamma section
- docs/spec/changelog/2026-06-08.md — changelog entry

Verified (matches Mathematica exactly):
- Attributes[EulerGamma] = {Constant, Protected}
- ?EulerGamma shows docstring
- N[EulerGamma]=0.577216; N[EulerGamma,50]/N[EulerGamma,250] full precision
- D[EulerGamma,x]=0; NumericQ[EulerGamma]=True
- Round[1/EulerGamma^100] = 734833795660954410469466
- RealDigits[EulerGamma,10,50,-10^4] = digits 10000-10049
- EulerGamma=5 blocked by Set::wrsym (Protected)
- core_tests pass (incl. test_clear_attributes, test_information)

Note: numeric eval (numeric.c, mpfr_const_euler) and derivative (deriv.c)
were already wired before this task; left untouched. Only EulerGamma was given
the Constant attribute (minimal impact); Pi/E/Catalan/etc. were not changed.
