# Task: Split number-theory builtins out of arithmetic.c into numbertheory.c

## Rationale
`arithmetic.c` (2318 lines) is a grab-bag of two concerns: (1) core
rational/complex constructors + widely-shared numeric helpers, and (2)
number-theoretic builtins. Split #2 into a sibling top-level file
`numbertheory.c`, matching the flat convention of `modular.c`/`facint.c`.
Introduce `numbertheory_init()` so registration moves out of `core.c`,
matching every other subsystem.

## Scope (moves to numbertheory.c)
Builtins: GCD, LCM, ExtendedGCD, PowerMod, Factorial, Factorial2,
FactorialPower, Binomial, PrimitiveRoot, PrimitiveRootList,
MultiplicativeOrder + their private static helpers.

## Stays in arithmetic.c
int64 `gcd()`/`lcm()` (used repo-wide), `g_arith_warnings_muted`,
make_rational/builtin_rational, make_complex/builtin_subtract/
builtin_complex/builtin_divide, and all shared helpers in arithmetic.h.

## Steps
- [ ] 1. Create `src/numbertheory.h` (guard, `numbertheory_init()`, 11 builtin protos)
- [ ] 2. Create `src/numbertheory.c` (preamble incl. + region A 34-410 + region B 931-EOF + `numbertheory_init()`)
- [ ] 3. Trim `src/arithmetic.c` to lines 1-33 + 411-930
- [ ] 4. Remove 11 NT builtin protos from `src/arithmetic.h` (keep builtin_add)
- [ ] 5. `src/core.c`: include numbertheory.h, replace registration block (401-411) with `numbertheory_init();`, remove attr/docstring block (439-450)
- [ ] 6. Build clean (`make -j`), smoke test REPL (GCD/Factorial/Binomial/PowerMod)
- [ ] 7. Run relevant tests; valgrind smoke
- [ ] 8. Docs: changelog note + SPEC.md layout line; update graph

## Notes
- internal.c uses local `extern` decls (lines 149-153) — unaffected.
- Regions A and B are independent; no static crosses the staying block.

## Review
Done. Outcome:
- New `src/numbertheory.{c,h}` (~2070 LoC) owns the 11 number-theory builtins
  + private helpers (modular-root, egcd, binomial-poly, primitive-root).
- `arithmetic.c`: 2318 → 305 LoC. Retains rational/complex constructors,
  int64 gcd/lcm, and shared numeric predicates (is_infinity_sym,
  expr_numeric_sign, is_neg_infinity_form, ...).
- `numbertheory_init()` added, called from `core_init()` after `facint_init()`;
  registration + attrs + FactorialPower docstring removed from `core.c`.
- 11 builtin protos moved arithmetic.h → numbertheory.h; `internal.c`'s local
  `extern` decls unaffected.
- `numbertheory.c` added to tests `COMMON_SRC`.
- Verified: clean `-Wall -Wextra` build; REPL smoke (all 12 forms incl.
  rational-exponent PowerMod, single-arg GCD, attrs Flat/Orderless); suites
  extended_gcd / primitive_root / modular / factorial_simplify / core PASS.

## Gotcha hit (logged to lessons.md)
Two non-`Expr*`-returning groups lived *inside* the number-theory span and my
first `^Expr\*|^static Expr` boundary grep missed both: PowerMod's `static int`
modular-root helpers (orig 639-929) and the `bool`/`int` numeric predicates
(orig 1406-1478, which are core helpers and had to stay). Always grep ALL
return types when carving regions.
