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

## Status line while evaluating (spinner)

An evaluation that runs longer than 200 ms gets one live line at the bottom
of the terminal, redrawn about ten times a second and erased before the
result:

```
⠹ Integrate · Together · 3.2s · cpu 3.1s · mem 412 MB (2.5%) · peak 1.17 GB
⠼ Table 3400/10000 ▕███▍      ▏ 34% · eta 6.4s · 2.1s · cpu 2.1s · mem 96 MB (0.6%)
```

- **What it shows:** the head of the input (the right-hand side of `x = ...`,
  the last call of a `;` sequence), the builtin running now, loop progress,
  elapsed and CPU time, resident memory as a share of the machine's RAM, the
  peak, and the number of running worker threads (FactorInteger, NDArray and
  multiprecision linear-algebra kernels). Memory turns yellow at 60 % and red
  at 85 %.
- **Progress:** `Do`, `Table`, `Sum`, `Product`, `Nest`, `NestList`, `Fold`,
  `FoldList`, `NestWhile` and `FixedPoint` report their step and total. For
  nested loops the outer loop and the one directly inside it are shown
  (`Table 2/3 ▸ Table 4/10`), and the bar is their combined fraction. An ETA
  appears after one second of progress. An unknown total (`Do[..., Infinity]`)
  shows a count and no bar. `Map` is not tracked: it builds `f[e1], f[e2], ...`
  and the evaluator evaluates them afterwards. Numeric fast paths (compiled
  `Do`/`Table` bodies) finish without reporting.
- **Output while it runs:** `Print`, messages and anything else written to
  stdout or stderr (including C libraries and `Run[]` children) go through the
  status thread, which erases the line, writes the output and redraws, so they
  never overwrite each other. The order of stdout and stderr is kept. A
  partial line (`WriteString` without a newline) is never drawn over.
- **Quick results** draw nothing: the output is byte-for-byte what it was
  without the spinner.
- **Narrow terminals:** the line never exceeds the width minus one column.
  Segments are dropped in this order: running builtin, peak, cpu, threads,
  eta, memory, progress, elapsed; then the label is cut. The width is re-read
  on every frame, so resizing works.
- **Ctrl-C** behaves exactly as described above. The line is erased and `^C`
  is printed before `Out[n]= $Aborted`. Typed-ahead text is not echoed while
  the line is shown; it appears at the next prompt.
- **When it is off:** `MATHILDA_SPINNER=0`; `TERM` unset, empty or `dumb`; any
  of stdin, stdout and stderr not a terminal, or stdout and stderr on
  different terminals; builds without threads or readline; pipe mode and
  `-file`. With `NO_COLOR` the line is plain and uses no escape sequences.
  Without a UTF-8 locale it uses ASCII (`| Integrate | 3.2s | ...`,
  `[###.......]`).

## The took note and the status line

- The dim note under a result that took at least one second also gives the
  CPU time and the peak resident memory during the evaluation:
  `took 3.21 s · cpu 3.18 s · peak 412 MB`.
- `MATHILDA_STATUSLINE=1` prints a dim line after every result:
  `last 0.02 s · mem 96 MB · Global``. It is off by default.

## For contributors: the terminal handover and progress hooks

- The line editor owns the terminal from the prompt until it returns a line.
  `term_status_begin()` and `term_status_end()` (`src/term/term_status.c`),
  called immediately around `evaluate` in `process_input`, own it during the
  evaluation. At both handover points the tty is in its cooked mode, fd 1 and
  fd 2 are the terminal, and the cursor is at column 0. A reader called
  during an evaluation (Dialog) calls `term_status_suspend()` and
  `term_status_resume()` around the read.
- A loop reports progress with `mth_progress_push(head, total)`, once per
  iteration `mth_progress_step(tok, done)` (one store), and
  `mth_progress_pop(tok)` (`src/term/term_progress.h`). A parallel driver
  brackets its threads with `mth_progress_workers_add(n)` and `(-n)`.
  `mth_progress_set_memory_limit()` is there for a future
  `MemoryConstrained`.
