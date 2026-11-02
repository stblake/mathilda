# Contributing to Mathilda

Thank you for your interest in Mathilda. This page is the short version of the
workflow. [`SPEC.md`](SPEC.md) describes the architecture and the coding
standards, and [`docs/extending.md`](docs/extending.md) has step-by-step recipes
for new builtins, modules, pattern constructs and operators.

## Build and test

```bash
make -j$(nproc)                                   # ./Mathilda
cd tests && mkdir -p build && cd build
cmake .. && make -j$(nproc) && ctest -LE nightly -j4
```

Any change to C code should also pass the portability and surface audits that
CI runs:

```bash
make check-c99            # strict C99 / glibc / MinGW portability
make check-messages       # every diagnostic routes through mth_message()
make check-packed-aware   # numeric heads with a buffer fast path opt in to packing
make check-compile-coverage
```

Check memory changes under `valgrind --leak-check=full ./Mathilda`, or use
AddressSanitizer.

## Adding a builtin

1. **Implement it in C** with the signature `Expr* builtin_myfunc(Expr* res)`.
   The builtin **takes ownership** of `res`:
   - On success, return a **new** `Expr*`. The evaluator frees `res`, so don't
     free it yourself.
   - If the input can't be evaluated (for example, symbolic arguments to a
     numeric function), return `NULL` and don't free `res`. The expression is
     then left unevaluated.
   - If you reuse a part of `res` in the result, set that slot to `NULL`
     first, so the part isn't freed twice.

   ```c
   Expr* builtin_myfunc(Expr* res) {
       if (res->data.function.arg_count != 1) return NULL;  /* leave unevaluated */
       /* ... */
       return expr_new_integer(42);
   }
   ```

2. **Register it** in the module's `*_init()`, with a docstring and attributes:

   ```c
   symtab_add_builtin("MyFunc", builtin_myfunc);
   symtab_set_docstring("MyFunc", "MyFunc[x]\n\tComputes the ultimate answer.");
   symtab_get_def("MyFunc")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
   ```

   Interned internal symbols belong in `src/sym_names.{c,h}`.

3. **Emit diagnostics** through `mth_message(head, tag, fmt, ...)`
   (`src/message.h`). Never write `fprintf(stderr, ...)` directly: a raw write
   escapes both `Quiet[]` and `Check[]`.

4. **Cover the numeric surfaces.** A builtin that operates on numbers also
   needs:
   - a packed/`NDArray` kernel
   - a `Compile[]` lowering
   - a machine-precision scalar path that reuses the kernel's helper

   See [`docs/design/packed_arrays.md`](docs/design/packed_arrays.md). If a
   surface genuinely doesn't apply, record that in the audit tool's exemption
   list.

5. **Test and document.** Add tests under `tests/`, document the function in
   the matching page under [`docs/spec/builtins/`](docs/spec/builtins/), and
   add an entry to the current week's
   [`docs/spec/changelog/`](docs/spec/changelog/) file. The file is named
   after that week's Monday (`YYYY-MM-DD.md`).

## Conventions

- **Strict C99.** The code must compile cleanly with
  `gcc -std=c99 -Wall -Wextra` on Linux, macOS and Windows (MinGW).
  `SPEC.md` §10 lists the common traps: POSIX-only symbols, `int64_t` versus
  `long long`, and `long` being 32 bits on Windows.
- **Versioning.** Every commit that changes behaviour bumps
  `MATHILDA_VERSION_NUMBER`/`MATHILDA_VERSION_STRING` in `src/version.h` (by
  `+0.0001`) and is tagged `v<version>`. Commits that only touch docs, CI or
  tooling are not bumped or tagged.
- **Commit messages.** Start with a type (`fix:`, `feat:`, `perf:`, `docs:`,
  `test:`, `chore:`) and a one-line summary, then explain what changed and
  why. End the subject with the version for bumped commits, as in `; v0.6141`.
- `src/external/` contains third-party code. Don't modify it.

[`CLAUDE.md`](CLAUDE.md) gives the full contributor checklist used when
developing with AI assistance.

## Reporting issues

Please open an issue at <https://github.com/stblake/mathilda/issues>. Include
a minimal input, the output you expected, and the output of `./Mathilda -v`.
