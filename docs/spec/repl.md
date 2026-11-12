# Interactive REPL

What the interactive session (`./Mathilda` with a terminal on stdin, built with
GNU Readline) does beyond reading and evaluating lines. None of it applies to
the NDJSON pipe protocol, to `-file` scripts, or to builds without readline
(`USE_READLINE=0`, Windows): those modes are byte-for-byte unchanged, install
no signal handler, print no colour, complete nothing and write no history
file.

The code is in `src/term/`: `term_lex` (input lexer), `term_interrupt` (the
interrupt flag and the SIGINT handler), `term_input` (the line-reader seam),
`term_history`, `term_complete`, `term_color`, and `term_readline`, the only
file that includes readline.

## Interrupting an evaluation (Ctrl-C)

- **While evaluating**, Ctrl-C ends the evaluation: the result is
  `Out[n]= $Aborted` and the session continues. The interrupt is the same
  abort `Abort[]` raises, so every frame's cleanup runs (`Block`/`Module`
  scoping is restored), `CheckAbort` catches it, and `AbortProtect` defers it
  to the end of the protected body.
- `TimeConstrained` does not mistake it for a timeout: its fail expression is
  not used.
- `Pause` ends within 50 ms (not inside `AbortProtect`).
- C kernels that poll the TimeConstrained deadline (CAD, Groebner F4/FLINT,
  FactorInteger, LU over MPFR, discrete logs, NSum, sparse solvers, compiled
  numeric loops) decline within their poll period. Messages raised while the
  abort unwinds are not printed (outside `AbortProtect`).
- An open Raylib window (`Plot`, `Manipulate`, `Animate`) closes.
- **Not reached:** code that never returns to the evaluator and never polls --
  the compiled VM's inner loop, a packed `Table`, LAPACK/BLAS calls, GMP-ECM
  stage 2, and the classical Buchberger fallback. For those:
- **A second Ctrl-C within one second**, still evaluating, ends the process
  with exit status 130.
- **At the prompt**, Ctrl-C discards the line being typed (including a
  multi-line one) and shows the same `In[n]:=` again. Two Ctrl-C at a prompt
  only clear the line twice.
- Inside `Dialog[]`, Ctrl-C aborts the dialog line being evaluated, and the
  dialog continues.

### For contributors: subprocesses

POSIX `system()` ignores SIGINT in the parent while it waits, so a Ctrl-C only
kills the child. Code that waits on a child (`system()`, `popen()`/`pclose()`,
any blocking wait) must turn the child's death into an interrupt:

```c
#include "term/term_interrupt.h"
if ((WIFSIGNALED(st) && WTERMSIG(st) == SIGINT) ||
    (WIFEXITED(st) && WEXITSTATUS(st) == 130))
    mth_interrupt_request();      /* the evaluation then ends as $Aborted */
```

`mth_interrupt_pending()` and `mth_interrupt_request()` are a stable API; both
are in `libmathilda.a`.

## History

- Every accepted input is appended at once to `$HOME/.mathilda_history`, so a
  crash or a forced quit loses nothing. Dialog lines are included.
- `MATHILDA_HISTFILE` names another file; set but empty, it disables the file
  (Up-arrow still recalls the session's own inputs).
- `MATHILDA_HISTSIZE` is the number of entries kept (default 1000; `0`
  disables history; invalid values give the default; at most 1,000,000). A
  longer file is cut to its newest entries when the session starts.
- Blank and comment-only inputs, and an input equal to the previous one, are
  not recorded.
- Format: the first line is `#mathilda-history v1`; each following line is one
  entry with `\` written `\\`, a line break `\n` and a carriage return `\r`,
  so multi-line inputs come back whole. New files are created with mode 0600.
  A non-empty file without the header is never written to; that, or a write
  failure, gives one `Mathilda: ...` note and history is kept for the session
  only.
- Concurrent sessions append to the same file; a session's start-up
  truncation can drop entries another session appends at that moment.

## Tab completion

Completion reads the symbol table directly. It never evaluates and never
creates a symbol.

- A prefix completes to defined symbols (builtins, symbols with attributes,
  values, rules or options); Module temporaries (`x$123`) are hidden. A prefix
  with a backtick (``System`Integr``) completes fully qualified names.
- Just after the `[` of a call, or after a `,` at that depth, the call head's
  options come first, as `Name -> ` (`ReplaceRepeated[x, {}, MaxI` gives
  `MaxIterations -> `), then symbols. String-named options complete inside an
  opening quote, as `Name" -> `.
- Nothing completes inside a comment or an ordinary string, and an empty word
  completes nothing outside an argument position.
- An ambiguous prefix completes to the common part; a second Tab lists the
  candidates.

## Colour and the elapsed-time note

- `In[n]:=` and `Out[n]=` and the `Head::tag` of messages are coloured.
- Rules, in order: `NO_COLOR` set and non-empty turns colour off;
  `CLICOLOR_FORCE` set and not `0` turns it on for both streams; `TERM` unset,
  empty or `dumb` turns it off; otherwise each of stdout and stderr is coloured
  only when it is a terminal (`./Mathilda 2>err.txt` leaves `err.txt` plain).
- After a result that took at least one second, a dim line such as
  `took 2.41 s` (or `took 1 min 15 s`) follows, when stdout is a terminal.

## libedit (macOS without Homebrew readline)

Completion and the history file work. A Ctrl-C at the prompt does not clear
the line (it is ignored, and never exits or submits). Prompts are not
coloured; `Out` and messages are.

## The line editor

On a terminal (stdin and stdout both ttys, `TERM` set and not `dumb`) the
session reads its input with Mathilda's own editor rather than readline. The
rest of this page up to here applies to both; the sections below are the
editor's.

- `MATHILDA_EDITOR=readline` (also `libedit`, `plain`, `none`, `off`, `0`)
  selects readline instead, as does `TERM=dumb`. Builds without readline
  (`USE_READLINE=0`) keep their plain line reader.
- The banner says "Tab completes, Ctrl-R searches the history, Ctrl-_
  undoes" when the editor is in use.

### Keys

| Key | Action |
|---|---|
| Return | evaluate when the input is complete, else a new line |
| Esc then Return, Alt-Return | evaluate now, complete or not |
| Up / Down, Ctrl-P / Ctrl-N | previous / next line; on the first / last line, previous / next history entry |
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
| Ctrl-_ (Ctrl-/) | undo (typing runs and deletion runs undo as one step) |
| Tab | completion (below) |
| Ctrl-R | history search (below) |
| Ctrl-L | clear the screen |
| Ctrl-C | discard the input, also a multi-line one; the same `In[n]:=` |
| Ctrl-D on an empty input | end the session |
| Ctrl-Z | suspend (`fg` resumes and redraws) |

- A pasted block (bracketed paste) is inserted as typed text: its line breaks
  never submit it, tabs become four spaces.
- Continuation lines are indented under the prompt on screen only; the
  input sent to the evaluator and to the history is exactly what was typed.
- Long lines wrap; an input taller than the screen scrolls with the cursor
  and is written out whole when it is evaluated.

### Colouring

As you type: builtins (cyan), symbols with a value or rules (blue), numbers
(magenta), strings (green), comments (grey) and patterns such as `x_`,
`_Integer`, `x__`, `#2` (yellow). The bracket at or before the cursor and its
partner are bold and underlined; a closing bracket that matches nothing, or the
wrong opener, is red. An unclosed opener is not marked: it is still being
typed. Results (`Out[n]=`) are coloured the same way, with the same characters
as uncoloured output. Colour follows the rules of "Colour" above; without it
the editor emits no SGR at all.

### Completion popup

- Tab with a single candidate inserts it. With several it inserts their
  common prefix and opens a list under the word (above it when the input
  sits at the bottom of the screen and has room above; otherwise the screen
  scrolls up).
- The row under the list shows the first line of the selected name's
  documentation (`?Name`), and a position counter when the list scrolls.
- Tab / Down / Ctrl-N and Shift-Tab / Up / Ctrl-P move the selection; Enter
  inserts it (it does not evaluate); Esc closes the list.
- Typing letters narrows the list; any other key closes it without
  inserting and then does what it normally does, so typing straight on never
  picks a name.
- Inside a call, the call's option names come first (`NIntegrate[f, {x, 0,
  1}, ` then Tab lists `AccuracyGoal -> `, ...).
- The list opens only on Tab, never by itself while typing.
- Below 40 columns the description row is left out; below 16 columns or 5
  rows there is no list, and each Tab inserts the next candidate instead.

### History search (Ctrl-R)

Ctrl-R shows `search:` and the history, newest first, each distinct input
once. Typing filters it (substring match, highlighted); a multi-line input is
shown whole. Ctrl-R / Up go to older matches, Ctrl-S / Down to newer ones.
Enter puts the selected input into the editor for editing (Enter again
evaluates it); Esc or Ctrl-G puts back what was there before; any other key
takes the selected input and then acts on it.

### The terminal

- While an input is being read the terminal is in raw mode, so Ctrl-C and
  Ctrl-Z are keys. Before an input is evaluated the terminal settings are
  restored exactly, so the evaluation (and Ctrl-C to interrupt it) sees the
  usual terminal.
- The terminal is also restored when the session exits, and if it is killed
  or crashes while reading (SIGTERM, SIGHUP, SIGQUIT, SIGSEGV, SIGBUS,
  SIGFPE, SIGILL, SIGABRT).
- A resize redraws the input at the new width. Wide (CJK, emoji) and
  combining characters take their proper columns.
- Keys typed while an evaluation runs become the start of the next input.

### For contributors

- The state machine (`term_editor.c`) is terminal-free: bytes in through
  `term_editor_feed`, drawing out through a write callback, one write per
  batch. `term_editor_tty.c` puts it on a terminal. `tests/term_vt.h` is a
  small VT emulator the unit tests draw into.
- Hand-over: the editor owns the terminal from the start of `term_read()` to
  its return. Anything that draws during an evaluation (a status line) owns
  it from then until the next `term_read()`, which can also come from a
  `Dialog[]` inside the evaluation.
- Performance gate: `tests/test_term_editor_latency.c` (200-line buffer,
  slowest keystroke under 2 ms, at most one write per keystroke, no
  allocation in steady-state editing).
