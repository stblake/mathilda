---
source: src/eval.c
---
**Definition.** `$AutoCompilation` is a boolean system variable that controls
whether Mathilda compiles numeric bodies to bytecode behind the scenes. It is
`True` by default; set it to `False` to force every such body through the
interpreter. It covers **both** automatic mechanisms at once: the auto-compile
adapter that compiles a held body once for many sample points (`Plot`, `Plot3D`,
`Table`, `NIntegrate`, `NSum`, `FindRoot`, the plot samplers) and the numeric-loop
compiler for `Do`, `For`, `While`, `Map`, `Nest`, `Fold` and `FixedPoint` bodies. A
`CompiledFunction` the user built with `Compile[]` is *not* affected — that was
asked for.

**Representation.** It is one row of the `EVAL_SYSFLAGS` table in `src/eval.c`,
which pairs the name with a setter and getter. The setter
(`eval_set_autocompilation`) flips the two underlying C flags —
`autocompile_set_enabled` and `numloop_set_enabled` — so one user-facing switch
drives both paths. The value is surfaced as an OwnValue (registered by
`eval_init_sysflags` with the live default) so it can be read back as well as
assigned; the assignment hook in `apply_assignment` routes `$`-prefixed writes
through `eval_sync_sysflag`. Only `True` and `False` are accepted; any other value
raises `$AutoCompilation::flagset` and leaves the flag unchanged.

**Usage & limits.** A compiled body is contracted to give the interpreter's answer,
so this changes **speed and nothing else** — it exists so the two paths can be
compared. It reads back `False` in a session started with the environment variable
`MATHILDA_NO_AUTOCOMPILE` set (the OwnValues are registered after the environment
overrides are read, so the reported value is never a lie).
