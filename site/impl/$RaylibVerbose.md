---
source: src/eval.c
---
**Definition.** `$RaylibVerbose` is a boolean system variable that controls whether
the Raylib graphics backend prints its diagnostic log — window and OpenGL
initialisation lines — to the terminal. It is `False` by default, so opening a
`Show`/`Plot` window or exporting an image stays quiet; set it to `True` to see
Raylib's full trace log.

**Representation.** It is one row of the `EVAL_SYSFLAGS` table in `src/eval.c`,
alongside `$AutoCompilation` and `$AutoArrayPacking`. Its setter/getter
(`raylib_verbose_set` / `raylib_verbose_enabled`) live in
`src/graphics/plot_common.c` and toggle a single raylib-free backing flag, so the
switch exists and is readable even in a build compiled without graphics. The value
is surfaced as an OwnValue by `eval_init_sysflags` so it can be read back; writes go
through the `$`-assignment hook (`eval_sync_sysflag`), which accepts only `True` or
`False` and otherwise raises `$RaylibVerbose::flagset` and keeps the current value.

**Usage & limits.** Affects only the terminal chatter from the renderer — never the
graphics produced. It has no visible effect in a build without graphics
(`USE_GRAPHICS=0`), where there is no Raylib log to gate, though the variable still
reads and assigns.
