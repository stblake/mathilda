---
source: src/numberform.c
---
**Algorithm.** `NumberForm` is a **print wrapper**: `builtin_numberform`
(`src/numberform.c`) is inert — it returns `NULL` so the `NumberForm[...]` head
survives in the tree (`FullForm[NumberForm[1.23, 2]]` is itself), and all the work
happens at print time. `print.c` installs an active `NumberFormCtx` built from the
precision spec and option rules (`build_ctx` / `nf_parse_spec` /
`nf_apply_option_rule`) and routes every numeric leaf through
`numberform_render_number`. Per number, `nf_format_parts`: rejects non-finite,
special-cases zero and exact integers; extracts `count` significant base-10 digits
plus a decimal exponent (`mpfr_get_str` for MPFR, `"%.*e"` for machine reals);
chooses scientific vs decimal (`ScientificNotationThreshold` / `ExponentFunction`
/ `ExponentStep`); lays the digits into integer and fractional strings (a `{n, f}`
spec re-rounds to `f` fractional digits, a plain `n` drops trailing zeros); then
applies `DigitBlock` grouping and assembles sign / `NumberPoint` /
`NumberMultiplier` / `NumberFormat`.

**Data structures.** A small growable string buffer (`SB`) builds each rendered
number; the `NumberFormCtx` holds the parsed spec, the option strings, and the
measured alignment field. The *same* `nf_format_parts` drives both the measure
pass (`nf_measure`, for padding/alignment widths) and the output pass, so the
printed widths can never disagree.

**Complexity / limits.** Linear in the digit count per number and in the leaf
count of the wrapped expression. A requested precision below the integer-digit
count issues `NumberForm::reqsigz` (once per print, via `mth_message`) and pads
with zeros. Because the head survives, an intervening `NumberForm` blocks
arithmetic on the surrounding expression. Attributes `NHoldRest`, `Protected` —
the first argument evaluates, the spec and options are held under numeric
evaluation.
