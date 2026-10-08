# Fix: Integrate emits two-argument ArcTan instead of one-arg form

`Integrate[x^(5/2) ArcTan[Sqrt[x]], x]` rendered its arctan part as
`6 ArcTan[1,√x] - 6 ArcTan[1,-√x]` (= `12 ArcTan[√x]`). Root cause: the Risch
trig/exp front-end's complex-split `cx_reim` emits `Arg[a+bi]` as the two-arg
`ArcTan[a,b]` (`src/calculus/risch_trig_frontend.c:153`).

## A. Code
- [ ] `src/calculus/risch_trig_frontend.c` — `cx_reim` `Log[a+bi]` case: collapse to
      one-arg `ArcTan[b/a]` (= `ArcTan[b]` when a==1) when `a` is a positive real
      constant; keep two-arg form otherwise (branch-sensitive). Use `evaluate()` +
      `expr_numeric_sign()`.

## B. Release / docs
- [ ] `src/version.h` — bump 0.302 → 0.303 (number + string)
- [ ] `docs/spec/changelog/2026-10-05.md` — integrate note
- [ ] `docs/spec/builtins/calculus.md` — one line on one-arg ArcTan rendering

## C. Test
- [ ] `tests/test_integrate_risch_transcendental.c` — regression assertions

## D. Verify
- [ ] `make -j` clean build; `make check-c99`
- [ ] repro family → one-arg; diff-back `Simplify[D-f]===0`; numeric unchanged @x=2
- [ ] shared-path regressions byte-identical (Sec, 1/(2+Cos), Sec^3, Tan, ArcTanh, x ArcTan[x], ArcSin)
- [ ] integration test suites green
- [ ] valgrind leak-free on repro

## Review

Done. One-line root-cause fix in `cx_reim` (`src/calculus/risch_trig_frontend.c`):
collapse `Arg[a+bi]` to one-arg `ArcTan[b/a]` (= `ArcTan[b]` when a==1) only when `a`
is a positive real constant (tested via `evaluate()` + `expr_numeric_sign()`); keep the
two-arg form otherwise (branch-sensitive). Fixes the whole `x^(n/2) ArcTan[√x]` /
`∫ArcTan[…]` family.

Verified:
- `∫x^(5/2) ArcTan[√x]` → `1/42(-6x+3x²-2x³+6Log[1+x]+12 x^(7/2) ArcTan[√x])`; `FreeQ[.,ArcTan[_,_]]==True`; diff-back `0`; numeric @2 = 2.86404 (unchanged).
- Family (`x^(3/2)`,`√x`,`∫ArcTan[x]`,`∫ArcTan[√x]`) all one-arg + diff-back 0.
- Regression set (`Sec`,`1/(2+Cos)`,`Sec³`,`Tan`,`ArcTanh`,`x ArcTan[x]`,`ArcSin`,`1/(1+x⁴)`,`1/(x³+1)`) byte-identical.
- Clean `make` build; `make check-c99` exit 0.
- 12 integration suites green incl. new `test_arctan_one_arg_output`; CRC corpus (1522) passes (3 pre-existing non-arctan `1/Sqrt[.Tan²]` diffs).
- valgrind: no new `definitely/indirectly lost` (byte-identical to trivial-eval baseline).
- v0.302 → v0.303; changelog + calculus.md updated.

Pending: git commit + tag `v0.303` (awaiting user go-ahead).
