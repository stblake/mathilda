# Advanced tutorial: Compilation & auto-compilation (+ two enabling fixes)

## Fixes made first (both verified, unit-tested)

- [x] **N over a packed integer array** now stays packed. `N[Range[10^6]]` used
      to unpack to a list of boxed reals (the gate handed the int64 buffer to N
      as a plain list; the element-by-element rebuild dropped packing). Fix:
      `N` claims `packed_int64_ok` (src/pack.c INT64_OK) + `numericalize_rec`'s
      EXPR_NDARRAY case widens an int64 buffer to float64 in one pass, inheriting
      presentation (src/numeric.c). Test: `test_n_over_integer_packs`.
- [x] **Integer auto-compiled loops** now run at machine speed. numloop's double
      VM refused an exact-integer accumulator (Real result, no overflow guard).
      Added an int64 twin execution of the imperative block (overflow-checked
      ci_*_i64; bail -> interpreter on overflow/rational/transcendental), wired
      into do_count / do_range / For. `Do[s=s+i,{i,1,10^7}]`: 6.2 s -> 0.070 s
      (89x), exact Integer, overflow -> bignum. Test: `test_int64_exact_loops`.

## Tutorial

- [x] Research: compile facilities, $AutoCompilation, $AutoArrayPacking, pipe
      protocol, verifier, machine = i9-9880H, Python 3.11.15 + NumPy 2.4.4.
- [x] All measurements taken (integer/double-recurrence/array/association),
      parity checked vs interpreter and vs NumPy/Python.
- [x] Write site/docs/tutorials/16-compilation.md (verified In/Out + measured tables).
- [x] Fix verify_tutorial.py (was vacuous — pipe mode emits NDJSON, 0 Out[]; now
      drives the pipe protocol and matches by id).
- [x] Verify the tutorial: `OK: 45 inputs, 0 mismatch(es)`.
- [x] Wire into .pages + index.md.
- [x] Build the site: `mkdocs build --strict` clean, no link warnings.
- [x] Changelog + docs/spec notes (control-flow.md, packed-arrays.md).

## Review

- **Both fixes verified end to end.** check-c99, check-packed-aware,
  check-array-exactness (0 MIXED) all clean; suites numeric, iter, association,
  packed_list (+new N test), ndarray*, eval_timestamps, compile_assoc, numloop
  (+new int64 test), autocompile, eval, compile, compiledfunction all pass.
- **Integer auto-compile**: `Do[s=s+i,{i,1,10^7}]` 6.2 s → 0.070 s (89×), exact
  Integer, overflow → bignum, rational → interpreter. int64 twin runner reuses
  the existing bytecode; strictly additive.
- **N packing**: `N[Range[10^6]]` now packed float64, values exact, MPFR path
  and real/complex buffers unchanged, `$AutoArrayPacking=False` respected.
- **Tutorial**: 45 verified transcripts; every perf claim measured vs NumPy 2.4.4
  / Python 3.11 on the i9-9880H. The verifier itself was broken (silently passing
  everything) and is now real.

## Lessons captured
- A "verification" tool that never sees output passes everything — verify the
  verifier. (verify_tutorial.py drove In[]/Out[]; the binary serves NDJSON.)
- Tutorial `In[k]:=` expressions must be single-line: the pipe protocol (and the
  verifier) send one line per expression; a wrapped `Compile[...]` truncates.
