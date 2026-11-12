# Interactive REPL

What the interactive session (`./Mathilda` with a terminal on stdin) does
beyond reading and evaluating lines: the line editor, completion, history,
colour themes, slash commands, the live status line while evaluating, and
user mods. None of it applies to the NDJSON pipe protocol, to `-file`
scripts, or to builds without readline (`USE_READLINE=0`, Windows): those
modes are byte-for-byte unchanged, install no signal handler, print no
colour, complete nothing, read no config file and write no history file.

```
In[1]:= Integrate[Sin[x]^2, {x, 0, Pi}]
Out[1]= Pi/2

In[2]:= Table[Pause[0.01]; Prime[k], {k, 400}];
⠼ Table 140/400 ▕███▌      ▏ 35% · eta 2.6s · 1.4s · cpu 0.0s · mem 96 MB (0.6%)
```

Code map: `src/term/` holds the terminal layer -- `term_editor` (the editor
state machine) with `term_editor_tty` (raw mode, signals), `term_render` (the
screen diff), `term_highlight` and `term_lex` (colouring), `term_width`,
`term_complete`, `term_history`, `term_interrupt` (Ctrl-C), `term_status`,
`term_spinner` and `term_progress` (the status line), `term_color` and
`term_theme` (colour), `term_config`, `term_slash`, `term_input` (the reader
seam) and `term_readline` (the only file that includes readline).
`src/repl_cmds.c` runs the slash commands and loads mods; `src/frontend.c`
holds the `SlashCommand` and `Monitor` builtins.

## Reading input

On a terminal (stdin and stdout both ttys, `TERM` set and not `dumb`) the
session uses Mathilda's own line editor. `MATHILDA_EDITOR=readline` (also
`libedit`, `plain`, `none`, `off`, `0`) or `TERM=dumb` selects GNU readline
(or libedit) instead; history, completion, Ctrl-C, colour, slash commands and
the status line work there too, without the editor's popups and colouring
while typing. Builds without readline keep their plain line reader.

Return evaluates when the input is complete; an open bracket, string or
comment continues it on the next line, and Esc then Return evaluates it
anyway.

### Keys (editor)

| Key | Action |
|---|---|
| Return | evaluate when the input is complete, else a new line |
| Esc then Return, Alt-Return | evaluate now, complete or not |
| Up / Down, Ctrl-P / Ctrl-N | previous / next screen row of a wrapped line, then line; from the first / last row, previous / next history entry |
| Alt-< / Alt-> | oldest history entry / back to what you were typing |
| Left / Right, Ctrl-B / Ctrl-F | character left / right |
| Alt-B / Alt-F, Ctrl- or Alt-arrows | word left / right |
| Home / End, Ctrl-A / Ctrl-E | start / end of the line |
| PageUp / PageDown | start / end of the input |
| Backspace, Delete, Ctrl-D | delete before / at the cursor |
| Ctrl-K, Ctrl-U | kill to the end / start of the line (Ctrl-K at a line's end joins the next) |
| Ctrl-W, Alt-Backspace, Alt-D | kill the word before (to whitespace / word characters), the word after |
| Ctrl-Y | yank the last kill (consecutive kills accumulate) |
| Ctrl-T | transpose characters |
| Ctrl-_ (Ctrl-/) | undo (a typing run or a deletion run is one step) |
| Alt-_ | redo what undo took back (until the next edit) |
| Tab | completion popup |
| Ctrl-R | history search |
| `/` in an empty input | the command menu |
| Ctrl-L | clear the screen |
| Ctrl-C | discard the input, also a multi-line one; the same `In[n]:=` |
| Ctrl-D on an empty input | end the session |
| Ctrl-Z | suspend (`fg` resumes and redraws the input) |

- A pasted block (bracketed paste) is inserted as typed text: its line breaks
  never submit it, tabs become four spaces.
- Continuation lines are indented under the prompt on screen only; the
  input sent to the evaluator and to the history is exactly what was typed.
- Long lines wrap; an input taller than the screen scrolls with the cursor
  and is written out whole when it is evaluated.
- Keys typed while an evaluation runs become the start of the next input.

### Colouring as you type

Builtins, symbols with a value or rules, numbers, strings, comments and
patterns (`x_`, `_Integer`, `x__`, `#2`) get their theme's colours. The
bracket at or before the cursor and its partner are bold and underlined; a
closing bracket that matches nothing, or the wrong opener, is red. An
unclosed opener is not marked: it is still being typed. Results (`Out[n]=`)
are coloured the same way, with exactly the characters of uncoloured output.

### Completion

- **Popup.** Tab with a single candidate inserts it. With several it inserts
  their common prefix and opens a list under the word (above it when the
  input sits at the bottom of the screen and has room above). The row under
  the list shows the first line of the selected name's documentation
  (`?Name`) and a position counter when the list scrolls. Tab / Down / Ctrl-N
  and Shift-Tab / Up / Ctrl-P move; Enter inserts (it does not evaluate); Esc
  closes. Typing letters narrows the list; any other key closes it and then
  does what it normally does. The list opens only on Tab. Below 40 columns
  the description row is left out; below 16 columns or 5 rows each Tab
  inserts the next candidate instead.
- **What completes.** Completion reads the symbol table; it never evaluates
  and never creates a symbol. A prefix completes to defined symbols
  (Module temporaries hidden); a prefix with a backtick (``System`Integr``)
  completes fully qualified names. Just after the `[` of a call or a `,` at
  that depth, the head's options come first as `Name -> `
  (`NIntegrate[f, {x, 0, 1}, ` lists `AccuracyGoal -> `, ...);
  string-named options complete inside an opening quote. Nothing completes
  inside a comment or an ordinary string, and an empty word completes nothing
  outside an argument position.
- With readline, an ambiguous prefix completes to the common part and a
  second Tab lists the candidates.

### History

- **Ctrl-R** shows `search:` and the history, newest first, each distinct
  input once; typing filters it (substring match, highlighted) and a
  multi-line input is shown whole. Ctrl-R / Up go to older matches, Ctrl-S /
  Down to newer ones. Enter puts the selected input into the editor (Enter
  again evaluates it); Esc or Ctrl-G puts back what was there before.
- **The file.** Every accepted input, slash commands and Dialog lines
  included, is appended at once to `$HOME/.mathilda_history`, so a crash or a
  forced quit loses nothing. `MATHILDA_HISTFILE` names another file (set but
  empty: no file, Up still recalls the session's inputs); `MATHILDA_HISTSIZE`
  is the number of entries kept (default 1000, `0` disables history, at most
  1,000,000). Blank and comment-only inputs and an input equal to the
  previous one are not recorded.
- **Format:** the first line is `#mathilda-history v1`; each following line
  is one entry with `\` written `\\`, a line break `\n` and a carriage return
  `\r`. New files are created with mode 0600. A non-empty file without the
  header is never written to; that, or a write failure, gives one
  `Mathilda: ...` note and history is kept for the session only. Concurrent
  sessions append to the same file; a session's start-up truncation can drop
  entries another session appends at that moment.

## Slash commands

A line that starts with `/` is a command. It is never evaluated: it gets no
`In[n]`/`Out[n]` and leaves `$Line` alone. (No Mathilda input can start with
`/`: `/.`, `//.`, `//`, `/@` and `/;` all need a left operand, so they are
unaffected, as is `/` anywhere else.)

| Command | |
|---|---|
| `/help` | the commands, the keys and the per-user files |
| `/theme [name]` | the theme picker, or choose a theme by name (below) |
| `/history [n]` | the last `n` (20) inputs, numbered |
| `/timing [on\|off\|always]` | the took note: for results of 1 s or more (on, the default), never, or after every result; no argument toggles on/off |
| `/statusline [on\|off]` | the summary line after each result (below); no argument toggles |
| `/clear` | clear the screen and the scrollback |
| `/mods` | the mods loaded, their status, commands and hooks |
| `/quit` | end the session, as `Quit[]` (`$Epilog` runs) |

Mods add their own (see Mods). In the editor, `/` typed into an empty input
opens the command menu: the completion popup over the commands, each with its
summary. Typing filters it (Backspace widens it again), Up/Down select, Tab
completes the name and adds a space for an argument, Return runs the selected
command, Esc closes the menu. Return always runs a line that starts with `/`,
complete or not. An unknown command says so; `/` alone shows `/help`.

## Colour and themes

**When there is colour.** `NO_COLOR` set and non-empty turns colour off;
`CLICOLOR_FORCE` set and not `0` turns it on for both streams; `TERM` unset,
empty or `dumb` turns it off; otherwise each of stdout and stderr is coloured
only when it is a terminal (`./Mathilda 2>err.txt` leaves `err.txt` plain).
Without colour the session writes no SGR sequence at all.

**Depth.** `MATHILDA_COLOR_DEPTH` (`24`, `256`, `16`, `0`) if set; else 24-bit
when `COLORTERM` is `truecolor` or `24bit` or `TERM` ends in `-direct`; 256
colours when `TERM` contains `256`; else the 16 basic colours. A 24-bit theme
colour becomes the nearest xterm-256 colour at 256, and the theme's own
16-colour choice below that.

**Themes.** Every colour the session draws -- prompts, results, message heads,
the status line, syntax colouring, brackets, popups -- comes from the theme.

| Theme | |
|---|---|
| `dark` | the default, for dark backgrounds |
| `light` | for light backgrounds (the default when `COLORFGBG` says the background is light) |
| `high-contrast` | bright, bold, never dim |
| `mono` | no colour: bold, dim, underline and reverse only |
| `ansi` | the terminal's own 16-colour palette |

`dark` uses its 24-bit colours only on a 24-bit terminal; below that it is
exactly the REPL's original 16-colour palette, which follows the terminal's
own (light or dark) colours. `ansi` is that palette at every depth. On a
light terminal with 24-bit colour, choose `light` once with `/theme`.

**Choosing.** `/theme` opens a picker: Up/Down (or `j`/`k`, or a digit) move
through the themes, and the picker and its preview -- a prompt, a result, a
message, the status line and a took note -- switch as you move. Return
applies the theme and saves it; Esc (or Ctrl-C, `q`) returns to the theme
you had. `/theme light` sets one directly. The choice is saved as
`theme = light` in `~/.mathilda/config`; `MATHILDA_THEME=name` overrides the
file at start-up. Readline sessions, without a picker, list the themes and
take `/theme name`.

## While evaluating: the status line

An evaluation that runs longer than 200 ms gets one live line at the bottom
of the terminal, redrawn about ten times a second and erased before the
result:

```
⠹ Integrate · Together · 3.2s · cpu 3.1s · mem 412 MB (2.5%) · peak 1.17 GB
⠼ Table 3400/10000 ▕███▍      ▏ 34% · eta 6.4s · 2.1s · cpu 2.1s · mem 96 MB (0.6%)
⠧ Do 41/100 ▕████      ▏ 41% · {41, 0.73} · 1.1s · cpu 1.1s · mem 20 MB (<0.1%)
```

- **What it shows:** the head of the input (the right-hand side of `x = ...`,
  the last call of a `;` sequence, the watched expression of a `Monitor`),
  the builtin running now, loop progress, elapsed and CPU time, resident
  memory as a share of the machine's RAM, the peak, and the number of running
  worker threads (FactorInteger, NDArray and multiprecision linear-algebra
  kernels). Memory turns yellow at 60 % and red at 85 %.
- **Progress:** `Do`, `Table`, `Sum`, `Product`, `Nest`, `NestList`, `Fold`,
  `FoldList`, `NestWhile` and `FixedPoint` report their step and total. For
  nested loops the outer loop and the one directly inside it are shown
  (`Table 2/3 ▸ Table 4/10`), and the bar is their combined fraction. An ETA
  appears after one second of progress; an unknown total (`Do[..., Infinity]`)
  shows a count and no bar. `Map` is not tracked (it builds `f[e1], ...` and
  the evaluator evaluates them afterwards), nor are numeric fast paths.
- **Monitor.** `Monitor[expr, mon]` evaluates `expr` and shows the current
  value of `mon` in the line (`Monitor[Do[..., {i, 100}], {i, x}]` above).
  `mon` is evaluated on the main thread when the line asks for a fresh value
  (about ten times a second), at the next loop step or the evaluator's
  periodic check, with its messages suppressed; a quick evaluation never
  evaluates it. Outside the interactive session `Monitor` just evaluates
  `expr`.
- **Output while it runs:** `Print`, messages and anything else written to
  stdout or stderr (including C libraries and `Run[]` children) go through the
  status thread, which erases the line, writes the output and redraws, so they
  never overwrite each other. The order of stdout and stderr is kept. A
  partial line (`WriteString` without a newline) is never drawn over.
- **Quick results** draw nothing: the output is byte-for-byte what it was
  without the status line.
- **Narrow terminals:** the line never exceeds the width minus one column.
  Segments are dropped in this order: running builtin, peak, cpu, threads,
  eta, memory, progress, monitor (first cut to half the line), elapsed; then
  the label is cut. The width is re-read on every frame.
- **When it is off:** `MATHILDA_SPINNER=0`; `TERM` unset, empty or `dumb`; any
  of stdin, stdout and stderr not a terminal, or stdout and stderr on
  different terminals; builds without threads or readline; pipe mode and
  `-file`. With `NO_COLOR` the line is plain and uses no escape sequences.
  Without a UTF-8 locale it uses ASCII (`| Integrate | 3.2s | ...`,
  `[###.......]`).

## After a result

- **The took note.** After a result that took at least one second, a dim line
  gives the time, the CPU time and the peak resident memory during the
  evaluation: `took 3.21 s · cpu 3.18 s · peak 412 MB` (`took 1 min 15 s` for
  longer ones). `/timing off` drops it, `/timing always` shows it after
  every result.
- **The status line** (`/statusline`, or `MATHILDA_STATUSLINE=1` at start-up)
  prints a dim summary after every result:
  `last 0.02 s · mem 96 MB · Global``. It is off by default.

## Interrupting and suspending

- **Ctrl-C while evaluating** ends the evaluation: the result is
  `Out[n]= $Aborted` and the session continues. The interrupt is the abort
  `Abort[]` raises, so every frame's cleanup runs (`Block`/`Module` scoping
  is restored), `CheckAbort` catches it and `AbortProtect` defers it to the
  end of the protected body. `TimeConstrained` does not take it for a
  timeout. `Pause` ends within 50 ms. C kernels that poll the
  TimeConstrained deadline (CAD, Groebner F4/FLINT, FactorInteger, LU over
  MPFR, discrete logs, NSum, sparse solvers, compiled numeric loops) stop
  within their poll period; messages raised while the abort unwinds are not
  printed. An open Raylib window (`Plot`, `Manipulate`, `Animate`) closes.
  With the status line showing, the line is erased and `^C` printed before
  the result.
- **Not reached** by Ctrl-C: code that never returns to the evaluator and
  never polls -- the compiled VM's inner loop, a packed `Table`, LAPACK/BLAS
  calls, GMP-ECM stage 2, the classical Buchberger fallback. For those, **a
  second Ctrl-C within one second** ends the process with status 130, the
  terminal restored.
- **Ctrl-C at the prompt** discards the input and shows the same `In[n]:=`.
  Inside `Dialog[]` it aborts the dialog line being evaluated.
- **Ctrl-Z** at the prompt restores the terminal and stops; `fg` puts the
  editor back in raw mode and redraws the input. During an evaluation the
  status line is erased and the terminal's echo restored before the stop;
  on `fg` the echo is turned off again and the line comes back with the next
  frame.

## Mods

At start-up, after the banner, the interactive session evaluates
`~/.mathilda/init.m` and then every `~/.mathilda/mods/NAME.m`, in name order.
A mod is ordinary Mathilda code. Besides defining functions it can:

- set the session hooks `$PreRead` (the input text), `$Pre` (each input
  before evaluation), `$Post` (each result) and `$PrePrint` (each result
  before printing), with their usual meaning;
- add commands: `SlashCommand["name", f, "summary"]` makes `/name rest`
  evaluate `f["rest"]` (the trimmed rest of the line) and print a non-`Null`
  result (a string as plain text). `SlashCommand["name", None]` removes it,
  `SlashCommand[]` lists all commands, and the built-in ones cannot be
  redefined. The command shows in the menu and in `/help` with its mod's
  name.

```
(* ~/.mathilda/init.m *)
SlashCommand["vars", Function[s, Names["Global`*"]], "list my symbols"];
SlashCommand["n", Function[s, N[ToExpression[s], 30]], "30 digits: /n Pi"];
```

A mod that reports messages, is aborted (Ctrl-C works during the load) or
evaluates `Quit[]` is reported in one dim line ("mods/x.m reported messages;
the session continues"), and the session goes on; `Quit[]` in a mod is
ignored. A failing command prints its messages, or `/name: aborted`, and
returns to the prompt. `/mods` lists each file with its status, the commands
it added and the hooks it set. `MATHILDA_MODS=0` skips the mods;
`MATHILDA_CONFIG_DIR` moves the whole directory (config, `init.m`, `mods/`),
and set but empty turns it off.

## The terminal

- While an input is being read the terminal is in raw mode, so Ctrl-C and
  Ctrl-Z are keys. Before an input is evaluated the settings are restored
  exactly; the status line then turns echo off for the evaluation (Ctrl-C
  still raises the interrupt) and back on after it.
- The terminal is restored when the session exits, and if it is killed or
  crashes, at the prompt or during an evaluation (SIGTERM, SIGHUP, SIGQUIT,
  SIGSEGV, SIGBUS, SIGFPE, SIGILL, SIGABRT), and after a double Ctrl-C.
- A resize redraws the input at the new width. Wide (CJK, emoji) and
  combining characters take their proper columns.
- The editor asks the terminal where the cursor is (`ESC [ 6 n`) to place
  popups. It never waits for a terminal that has not answered yet: the
  first query's reply is read among the keys, and a terminal that leaves two
  queries unanswered is not asked again.

## libedit (macOS without Homebrew readline, `MATHILDA_EDITOR=readline`)

Completion and the history file work. A Ctrl-C at the prompt does not clear
the line (it is ignored, and never exits or submits). Prompts are not
coloured; `Out` and messages are. Multi-line input does not work on libedit.

## Environment

| Variable | Effect |
|---|---|
| `MATHILDA_EDITOR` | `readline` (and synonyms) uses readline instead of the editor |
| `MATHILDA_HISTFILE`, `MATHILDA_HISTSIZE` | the history file and its size |
| `NO_COLOR`, `CLICOLOR_FORCE` | colour off / on |
| `MATHILDA_COLOR_DEPTH`, `COLORTERM` | colour depth |
| `MATHILDA_THEME`, `COLORFGBG` | the theme at start-up |
| `MATHILDA_CONFIG_DIR` | the per-user directory (default `~/.mathilda`; empty: none) |
| `MATHILDA_MODS` | `0` skips `init.m` and the mods |
| `MATHILDA_SPINNER`, `MATHILDA_SPINNER_DELAY_MS` | `0` turns the status line off; its delay (default 200) |
| `MATHILDA_STATUSLINE` | `1` turns the after-result status line on |

## For contributors

- **The terminal hand-over.** The editor owns the terminal from the start of
  `term_read()` until it returns (raw mode, bracketed paste, its resize and
  crash handlers). `term_status_begin()` / `term_status_end()`, called
  immediately around `evaluate` in `process_input`, give it to the status
  thread for the evaluation, which points fd 1 and fd 2 at a pipe it drains.
  At both points the tty is cooked, fds 1/2 are the terminal and the cursor
  is at column 0. A read during an evaluation (`Dialog[]`) is bracketed by
  `term_status_suspend()` / `term_status_resume()` in `dialog_read_term`.
  The status module's fatal-signal and SIGTSTP handlers are installed once
  at start-up; the editor's per-read crash handlers chain to them.
- **Progress hooks.** A loop reports with `mth_progress_push(head, total)`,
  `mth_progress_step(tok, done)` once per iteration (a store and one load of
  the Monitor flag) and `mth_progress_pop(tok)` (`src/term/term_progress.h`).
  A parallel driver brackets its threads with `mth_progress_workers_add(n)`
  and `(-n)`. `mth_progress_set_memory_limit()` is there for a future
  `MemoryConstrained`.
- **Subprocesses.** POSIX `system()` ignores SIGINT in the parent while it
  waits, so a Ctrl-C only kills the child. Code that waits on a child must
  turn the child's death into an interrupt:

  ```c
  #include "term/term_interrupt.h"
  if ((WIFSIGNALED(st) && WTERMSIG(st) == SIGINT) ||
      (WIFEXITED(st) && WEXITSTATUS(st) == 130))
      mth_interrupt_request();      /* the evaluation then ends as $Aborted */
  ```

  `mth_interrupt_pending()` and `mth_interrupt_request()` are a stable API in
  `libmathilda.a`.
- **The editor** (`term_editor.c`) is terminal-free: bytes in through
  `term_editor_feed`, drawing out through a write callback, one write per
  batch. `term_editor_tty.c` puts it on a terminal; `tests/term_vt.h` is the
  small VT emulator the unit tests draw into. `term_editor_pick_begin` runs a
  list picker (the theme menu) in the same frame. Performance gate:
  `tests/test_term_editor_latency.c` (200-line buffer, slowest keystroke
  under 2 ms, at most one write per keystroke, no allocation while typing).
- **Colours** are theme roles (`TermRole` in `term_theme.h`); add a role to
  every theme and map the style to it rather than writing an SGR literal.
- **A new built-in command** is one `term_slash_add_builtin` line and a case
  in `repl_cmds_run` (`src/repl_cmds.c`).
- **Tests:** `term_*_tests` unit suites, and `tests/pty/repl_pty_test.py`,
  which drives the real binary on a pseudo-terminal (it pins
  `MATHILDA_COLOR_DEPTH=16`; its Ctrl-Z cases run Mathilda under a small
  job-control wrapper, because a pty child's orphaned process group ignores
  stop signals).
