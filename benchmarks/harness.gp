\\ ==========================================================================
\\ benchmarks/harness.gp -- the PARI/GP twin of harness.m / harness.py.
\\ ==========================================================================
\\
\\ Loaded by every .gp experiment as its first line:
\\
\\     read("../harness.gp");
\\
\\ `read` resolves against the current directory; run_all.py sets cwd to the
\\ experiment folder (as it does for the .m half).  By hand:
\\
\\     cd benchmarks/100-rsa-keygen && gp -q -f rsa_keygen.gp
\\
\\ THE EMIT CONTRACT is identical to the other two harnesses -- one tagged,
\\ tab-separated line per case, joined across systems by <label>:
\\
\\     BENCH   <label>  <milliseconds>
\\     CHECK   <label>  <value>
\\     REQUIRE <head>   present|absent
\\     SKIP    <label>  <head>
\\
\\ Timing: one untimed warm-up, then the MINIMUM of reps timed runs, wall clock
\\ (getabstime), matching AbsoluteTiming / perf_counter on the other sides.
\\ getabstime() has MILLISECOND resolution, so every group-F case is sized to
\\ run >= 20 ms (loop inside the case) -- a sub-ms PARI case would read 0.000
\\ and drop out of the ratio as "no baseline".
\\
\\ A case body is a closure, ()->(...).  Every experiment file must end with
\\ quit(), or gp drops into the interactive prompt and the runner times it
\\ out.
\\ ==========================================================================

BENCH_REPS = 3;

bench(label, f, reps = BENCH_REPS) =
{
  my(t, best = -1);
  f();                                            \\ warm-up, untimed
  for (i = 1, reps,
    t = getabstime(); f(); t = getabstime() - t;
    if (best < 0 || t < best, best = t));
  printf("BENCH\t%s\t%.3f\n", label, best);
}

\\ Sub-millisecond cases outside group F: repeat each timed run until it spans
\\ >= 100 ms of the ms clock and divide, so the case does not read 0.000.
benchLoop(label, f, reps = BENCH_REPS) =
{
  my(t, n, best = -1.);
  f();                                            \\ warm-up, untimed
  for (i = 1, reps,
    n = 0; t = getabstime();
    until (getabstime() - t >= 100, f(); n++);
    t = (getabstime() - t) / n;
    if (best < 0 || t < best, best = t));
  printf("BENCH\t%s\t%.3f\n", label, best);
}

\\ Single timed run, no warm-up: for cases already expensive on their own.
benchOnce(label, f) =
{
  my(t = getabstime());
  f();
  printf("BENCH\t%s\t%.3f\n", label, getabstime() - t);
}

\\ The value gate.  Values must print on one line; group F reduces every check
\\ to an integer (or a short list of integers) so text compares exactly.
check(label, v) = printf("CHECK\t%s\t%s\n", label, v);

\\ PARI ships everything group F asks of it; `require` exists so the three
\\ halves stay symmetrical and the register reads the same.
require(heads) =
{
  if (type(heads) != "t_VEC", heads = [heads]);
  for (i = 1, #heads, printf("REQUIRE\t%s\tpresent\n", heads[i]));
}

\\ A case PARI genuinely cannot run (no such primitive): emit SKIP so the
\\ runner records the PARI column as empty rather than "never measured".
skip(label, head) = printf("SKIP\t%s\t%s\n", label, head);

\\ LAST, on purpose: raising parisizemax reallocates the stack, and gp abandons
\\ the rest of the file being read when it does -- placed first, every
\\ function above silently failed to define.
default(debugmem, 0);
default(parisizemax, 2^32);
