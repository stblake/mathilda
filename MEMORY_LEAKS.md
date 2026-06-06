# Known Memory Leaks

A running log of **pre-existing** memory leaks discovered in the codebase but
not yet fixed (because the fix is out of scope for the change that surfaced
them). Each entry should give enough detail to fix it later without re-doing
the investigation.

How to reproduce a leak check on macOS (valgrind reports a constant ~12,800 B /
400 blocks of dyld/libobjc/Accelerate init noise — always diff against a
`Sin[1.0]` baseline and only count stacks that reference `src/` files):

```bash
printf 'Sin[1.0]\nQuit[]\n'      | valgrind --leak-check=full ./Mathilda 2>baseline.txt
printf '<expr>\nQuit[]\n'        | valgrind --leak-check=full ./Mathilda 2>run.txt
grep "definitely lost:" baseline.txt run.txt        # compare totals
grep -E "\.c:[0-9]+\)" run.txt                       # inspect src frames
```

On macOS without valgrind, `leaks --atExit -- ./Mathilda` is the faster check;
diff the `N leaks for M total leaked bytes` line against a `1+1` baseline. Note
that an `-O3` build inlines aggressively and mis-attributes the leak root to the
outermost inlined frame — rebuild the suspect TUs with `-O0 -fno-inline` before
trusting the stack (this is how the 2026-06-06 `draw_int_range` leak, originally
mis-blamed on `numericalize`, was correctly located).

---

_No known unfixed leaks. The three previously tracked here
(`builtin_inequality`, `intsimp_log_to_arctanh`, and the `PossibleZeroQ`
Schwartz–Zippel sampler) were all fixed on 2026-06-06 — see
`docs/spec/changelog/2026-06-01.md`._
