# Implement ChineseRemainder

## Tasks
- [x] Create `src/numbertheory/chineseremainder.c` (streaming GMP CRT merge)
- [x] Add prototype to `src/numbertheory/numbertheory.h`
- [x] Register + ATTR_PROTECTED in `src/numbertheory/numbertheory.c`
- [x] Add `SYM_ChineseRemainder` to `src/sym_names.h` + `src/sym_names.c` (3 sites)
- [x] Add docstring in `src/info.c`
- [x] Docs: `docs/spec/builtins/number-theory.md` + changelog `docs/spec/changelog/2026-08-24.md`
- [x] Build REPL, smoke-test doc examples (all 12 pass; bignum + error msgs verified)
- [x] `make check-c99` — clean
- [ ] Build & run unit test (in progress); valgrind clean
- [ ] `make check-fastpath-sweep` (OFF_BUFFER only if flagged)

## Review

Implemented `ChineseRemainder` as a streaming pairwise CRT fold (all GMP) in
`src/numbertheory/chineseremainder.c`. Two-arg and offset (3-arg) forms; handles
non-coprime-but-consistent moduli, returns unevaluated on inconsistent systems /
zero modulus / length mismatch / non-integer / symbolic args. `Protected`, not
`Listable`. Registered symbol + docstring + spec/changelog.

Verification:
- Main `make` build: clean (gcc-16, no warnings). `make check-c99`: clean.
- REPL smoke test: all doc examples correct incl. 27-digit bignum + Mod
  round-trip; `ChineseRemainder::argt` on wrong arity.
- Unit test (`chinese_remainder_tests`): 27/27 pass.
- valgrind: leak-free. Definitely-lost delta vs untouched `extended_gcd_tests`
  baseline is +16 B / +2 blocks, attributable to the PRE-EXISTING `LCM` /
  `rational_like_to_mpz_pair` double-init (my tests call `LCM[...]`); zero leak
  frames touch `chineseremainder`.
- `make check-packed-aware`: OK. `make check-fastpath-sweep`: (running).

FIXED: the pre-existing `GCD`/`LCM` leak in `src/numbertheory/nt_util.c`
`rational_like_to_mpz_pair` (called `expr_to_mpz` on already-init'd targets;
`expr_to_mpz` re-inits, leaking prior limbs each fold step). Now extracts via a
temporary + `mpz_set`. valgrind: ChineseRemainder now matches the ExtendedGCD
baseline exactly (13,504 B / 422 blocks, pure macOS ObjC noise); GCD/LCM results
unchanged (incl. rational and 28-digit bignum paths).
