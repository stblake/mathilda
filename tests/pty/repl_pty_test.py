#!/usr/bin/env python3
"""Drive the Mathilda binary through a pseudo-terminal.

This is the only test of the real interactive terminal path (readline, the
SIGINT handler, the history file, Tab completion, colour).  It also pins the
non-interactive modes -- the NDJSON pipe protocol and `-file` scripts -- to
golden outputs recorded from the binary before the interactive work began, so
those modes are proven byte-identical on every run.

Usage:
    python3 tests/pty/repl_pty_test.py [--bin ./Mathilda] [-k NAME]
                                       [--golden DIR] [--record DIR]
                                       [--wrapper "valgrind ..."]

  --bin      the binary under test (default ./Mathilda).  Missing -> exit 77.
  --record   write the golden files for the invariance cases into DIR, then exit.
  --golden   where the golden files live (default tests/pty/golden).
  --wrapper  a command prefix (valgrind, leaks); every timeout is multiplied by 10.
  -k NAME    run only the cases whose name matches the regex NAME.

Exit status: 0 when every selected case passes, 1 on any failure, 77 when the
binary is missing or the host has no pty (ctest SKIP_RETURN_CODE).

MATHILDA_PTY_SELFTEST=1 adds the `harness_can_fail` case, which expects a wrong
answer and must report FAIL: it proves the harness is able to fail.

Only the Python 3 standard library is used.
"""

import argparse
import fcntl
import os
import re
import select
import shlex
import signal
import struct
import subprocess
import sys
import tempfile
import termios
import threading
import time

if os.name != "posix":
    print("repl_pty_test: no pty on this host; skipping")
    sys.exit(77)

import pty  # noqa: E402  (POSIX only)

ROWS, COLS = 24, 100
TIME_SCALE = 1.0          # multiplied by 10 under --wrapper
WRAPPER = []              # command prefix list
BIN = "./Mathilda"


class CaseFailure(Exception):
    pass


def check(cond, what):
    if not cond:
        raise CaseFailure(what)


# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------

def base_env(home, term="xterm-256color", extra=None, drop=()):
    """A minimal, reproducible environment: HOME in a tempdir, no windows."""
    env = {}
    for k in ("PATH", "LANG", "LC_ALL", "LC_CTYPE", "SDKROOT", "TMPDIR",
              "DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"):
        if k in os.environ:
            env[k] = os.environ[k]
    env["HOME"] = home
    env["MATHILDA_NO_WINDOW"] = "1"
    # The colour cases pin the 16-colour SGR codes (the default theme at 16
    # colours); the theme cases set their own depth.
    env["MATHILDA_COLOR_DEPTH"] = "16"
    if term is not None:
        env["TERM"] = term
    if extra:
        env.update(extra)
    for k in drop:
        env.pop(k, None)
    return env


def argv(args):
    return WRAPPER + [BIN] + list(args)


# ---------------------------------------------------------------------------
# Interactive sessions on a pty
# ---------------------------------------------------------------------------

class Session:
    """A child process whose stdin/stdout/stderr are one pty (pty.fork)."""

    def __init__(self, env, args=(), rows=ROWS, cols=COLS):
        cmd = argv(args)
        pid, fd = pty.fork()
        if pid == 0:  # child
            try:
                fcntl.ioctl(0, termios.TIOCSWINSZ,
                            struct.pack("HHHH", rows, cols, 0, 0))
                os.execvpe(cmd[0], cmd, env)
            finally:
                os._exit(127)
        self.pid = pid
        self.fd = fd
        self.buf = b""        # everything read so far (the transcript)
        self.pos = 0          # expect() searches from here
        self.status = None

    def _read_some(self, timeout):
        r, _, _ = select.select([self.fd], [], [], max(0.0, timeout))
        if not r:
            return False
        try:
            data = os.read(self.fd, 65536)
        except OSError:      # EIO: the child closed the pty
            data = b""
        if not data:
            raise EOFError
        self.buf += data
        return True

    def expect(self, pattern, timeout=10.0):
        """Wait for regex `pattern` (bytes) after the previous match."""
        rx = re.compile(pattern if isinstance(pattern, bytes)
                        else pattern.encode())
        deadline = time.monotonic() + timeout * TIME_SCALE
        while True:
            m = rx.search(self.buf, self.pos)
            if m:
                self.pos = m.end()
                return m
            left = deadline - time.monotonic()
            if left <= 0:
                break
            try:
                self._read_some(min(left, 0.2))
            except EOFError:
                m = rx.search(self.buf, self.pos)
                if m:
                    self.pos = m.end()
                    return m
                raise CaseFailure("EOF while waiting for %r\n--- transcript tail ---\n%s"
                                  % (pattern, tail(self.buf)))
        raise CaseFailure("timeout (%.1fs) waiting for %r\n--- transcript tail ---\n%s"
                          % (timeout * TIME_SCALE, pattern, tail(self.buf)))

    def drain(self, seconds):
        """Read whatever arrives within `seconds`."""
        deadline = time.monotonic() + seconds * TIME_SCALE
        while True:
            left = deadline - time.monotonic()
            if left <= 0:
                return
            try:
                self._read_some(min(left, 0.1))
            except EOFError:
                return

    def send(self, data):
        if isinstance(data, str):
            data = data.encode()
        os.write(self.fd, data)

    def finish(self, timeout=10.0):
        """Wait for the child to exit; returns the waitpid status."""
        deadline = time.monotonic() + timeout * TIME_SCALE
        while self.status is None:
            pid, st = os.waitpid(self.pid, os.WNOHANG)
            if pid:
                self.status = st
                break
            if time.monotonic() > deadline:
                self.kill()
                raise CaseFailure("child did not exit within %.1fs\n%s"
                                  % (timeout * TIME_SCALE, tail(self.buf)))
            try:
                self._read_some(0.05)
            except EOFError:
                time.sleep(0.02)
        self.drain(0.1)
        return self.status

    def exit_code(self, timeout=10.0):
        st = self.finish(timeout)
        if os.WIFEXITED(st):
            return os.WEXITSTATUS(st)
        return 128 + os.WTERMSIG(st)

    def kill(self):
        try:
            os.kill(self.pid, signal.SIGKILL)
            os.waitpid(self.pid, 0)
        except OSError:
            pass

    def close(self):
        if self.status is None:
            self.kill()
            self.status = -1
        try:
            os.close(self.fd)
        except OSError:
            pass

    def text(self):
        return self.buf.decode("utf-8", "replace")


def tail(buf, n=1500):
    return buf[-n:].decode("utf-8", "replace")


def wait_prompt(s, n, timeout=30.0):
    return s.expect(r"In\[%d\]:= " % n, timeout)


def new_session(home, **kw):
    # Plain output unless a case asks for colour (drop=("NO_COLOR",)), so the
    # regexes below match the text, not the SGR codes around it.
    extra = dict({"NO_COLOR": "1"}, **(kw.pop("extra", None) or {}))
    term = kw.pop("term", "xterm-256color")
    drop = kw.pop("drop", ())
    return Session(base_env(home, term=term, extra=extra, drop=drop), **kw)


# ---------------------------------------------------------------------------
# Non-interactive runs (pipe protocol and -file)
# ---------------------------------------------------------------------------

PIPE_INPUTS = [
    "1+1",
    'Print["hi"]; 3',
    "1/0",
    "Abort[]",
    "x = 5; x^2",
    "f[[",                     # parse error
    "Pause[0.01]",
    "?Plus",
]

FILE_SCRIPT = """\
Print[1 + 1]
Print["hi"]; 3
1/0
x = 5; Print[x^2]
Pause[0.01]
Print[CheckAbort[Abort[], caught]]
?Plus
Print["before abort"]
Abort[]
Print["never printed"]
"""

FILE_SYNTAX_SCRIPT = """\
Print["first"]
f[[
Print["never printed"]
"""


def pipe_request_text():
    import json
    return "".join(json.dumps({"id": i + 1, "expr": e}) + "\n"
                   for i, e in enumerate(PIPE_INPUTS))


def normalise(data):
    data = data.replace(b"\r\n", b"\n")
    data = re.sub(rb'"memory":\s*\d+', b'"memory":N', data)
    return data


def run_capture(args, env, stdin_data=b"", stderr_pty=False, timeout=60.0):
    """Run with piped stdin/stdout; stderr a pipe or (stderr_pty) a pty."""
    if stderr_pty:
        master, slave = pty.openpty()
        p = subprocess.Popen(argv(args), stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=slave, env=env)
        os.close(slave)
        # Drain the pty while the child runs: macOS discards unread pty
        # output once the last slave descriptor closes.
        chunks = []

        def pump():
            while True:
                try:
                    chunk = os.read(master, 65536)
                except OSError:      # EIO: the child closed the pty
                    return
                if not chunk:
                    return
                chunks.append(chunk)

        reader = threading.Thread(target=pump, daemon=True)
        reader.start()
        out = p.communicate(stdin_data, timeout=timeout * TIME_SCALE)[0]
        reader.join(5.0 * TIME_SCALE)
        os.close(master)
        return p.returncode, out, b"".join(chunks)
    p = subprocess.run(argv(args), input=stdin_data, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, env=env,
                       timeout=timeout * TIME_SCALE)
    return p.returncode, p.stdout, p.stderr


def run_file_tty(script_path, env, timeout=60.0):
    """`-file` with stdin and stdout on a pty, stderr a pipe."""
    master, slave = pty.openpty()
    fcntl.ioctl(slave, termios.TIOCSWINSZ, struct.pack("HHHH", ROWS, COLS, 0, 0))
    p = subprocess.Popen(argv(["-file", script_path]), stdin=slave, stdout=slave,
                         stderr=subprocess.PIPE, env=env)
    os.close(slave)
    out = b""
    deadline = time.monotonic() + timeout * TIME_SCALE
    while True:
        r, _, _ = select.select([master], [], [], 0.2)
        if r:
            try:
                chunk = os.read(master, 65536)
            except OSError:
                chunk = b""
            if not chunk:
                break
            out += chunk
        elif p.poll() is not None:
            break
        if time.monotonic() > deadline:
            p.kill()
            raise CaseFailure("-file under a pty did not finish")
    err = p.stderr.read()
    p.wait()
    os.close(master)
    return p.returncode, out, err


def golden_runs(tmp):
    """The invariance runs: {golden file name: bytes}, plus exit codes."""
    home = os.path.join(tmp, "home")
    os.makedirs(home, exist_ok=True)
    env = base_env(home)
    res = {}
    rc, out, err = run_capture([], env, pipe_request_text().encode())
    res["pipe_basic.ndjson"] = normalise(out)
    res["pipe_basic.stderr"] = normalise(err)
    res["pipe_basic.rc"] = str(rc).encode()
    for name, script in (("file_basic", FILE_SCRIPT),
                         ("file_syntax", FILE_SYNTAX_SCRIPT)):
        path = os.path.join(tmp, name + ".m")
        with open(path, "w") as f:
            f.write(script)
        rc, out, err = run_capture(["-file", path], env)
        res[name + ".stdout"] = normalise(out)
        res[name + ".stderr"] = normalise(err.replace(path.encode(), b"SCRIPT"))
        res[name + ".rc"] = str(rc).encode()
    return res


GOLDEN_DIR = "tests/pty/golden"


def golden(name):
    with open(os.path.join(GOLDEN_DIR, name), "rb") as f:
        return f.read()


def compare(name, got):
    want = golden(name)
    if got != want:
        raise CaseFailure("%s differs from the golden\n--- want ---\n%s\n--- got ---\n%s"
                          % (name, want.decode("utf-8", "replace")[-1200:],
                             got.decode("utf-8", "replace")[-1200:]))


# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

CASES = []


def case(fn):
    CASES.append(fn)
    return fn


@case
def invariance_pipe(tmp):
    res = golden_runs(tmp)
    for name in ("pipe_basic.ndjson", "pipe_basic.stderr", "pipe_basic.rc"):
        compare(name, res[name])
    check(b"\x1b" not in res["pipe_basic.ndjson"] + res["pipe_basic.stderr"],
          "escape sequence in pipe-mode output")


@case
def invariance_pipe_tty_stderr(tmp):
    home = os.path.join(tmp, "home")
    os.makedirs(home)
    rc, out, err = run_capture([], base_env(home), pipe_request_text().encode(),
                               stderr_pty=True)
    check(b"\x1b" not in err, "escape sequence on a tty stderr in pipe mode: %r" % err)
    compare("pipe_basic.ndjson", normalise(out))
    compare("pipe_basic.stderr", normalise(err))


@case
def invariance_file(tmp):
    res = golden_runs(tmp)
    for stem in ("file_basic", "file_syntax"):
        for ext in (".stdout", ".stderr", ".rc"):
            compare(stem + ext, res[stem + ext])
    check(b"\x1b" not in res["file_basic.stdout"] + res["file_basic.stderr"],
          "escape sequence in -file output")


@case
def invariance_file_tty(tmp):
    home = os.path.join(tmp, "home")
    os.makedirs(home)
    path = os.path.join(tmp, "file_basic.m")
    with open(path, "w") as f:
        f.write(FILE_SCRIPT)
    rc, out, err = run_file_tty(path, base_env(home))
    check(b"\x1b" not in out + err, "escape sequence in -file output on a tty")
    compare("file_basic.stdout", normalise(out))
    compare("file_basic.stderr", normalise(err.replace(path.encode(), b"SCRIPT")))
    check(os.listdir(home) == [], "-file created files in $HOME: %r" % os.listdir(home))


@case
def smoke_interactive(tmp):
    s = new_session(tmp)
    try:
        wait_prompt(s, 1)
        s.send("1+1\r")
        s.expect(r"Out\[1\]= 2")
        wait_prompt(s, 2)
        s.send("\x04")
        check(s.exit_code() == 0, "Ctrl-D did not exit 0")
    finally:
        s.close()


@case
def smart_return_continuation(tmp):
    """An open bracket continues the line; closing it evaluates the whole."""
    s = new_session(tmp)
    try:
        wait_prompt(s, 1)
        s.send("f[\r")
        s.drain(0.3)
        check(b"Out[1]" not in s.buf, "f[ was submitted")
        s.send("1]\r")
        s.expect(r"Out\[1\]= f\[1\]")
        wait_prompt(s, 2)
        s.send("\x04")
        check(s.exit_code() == 0, "Ctrl-D did not exit 0")
    finally:
        s.close()


def interactive(tmp, **kw):
    s = new_session(tmp, **kw)
    wait_prompt(s, 1)
    return s


def ctrlc_after(s, input_text, delay, n=1, timeout=2.0):
    """Send input, wait `delay`, press Ctrl-C; expect Out[n]= $Aborted.
    Returns the seconds from the key to the result."""
    s.send(input_text + "\r")
    s.drain(delay)
    check(b"Out[%d]" % n not in s.buf[s.pos:], "%s finished before Ctrl-C" % input_text)
    t0 = time.monotonic()
    s.send("\x03")
    s.expect(r"Out\[%d\]= \$Aborted" % n, timeout)
    return time.monotonic() - t0


def then_still_works(s, n):
    wait_prompt(s, n)
    s.send("1+1\r")
    s.expect(r"Out\[%d\]= 2" % n)
    wait_prompt(s, n + 1)
    s.send("\x04")
    check(s.exit_code() == 0, "Ctrl-D did not exit 0")


@case
def banner(tmp):
    s = interactive(tmp)
    try:
        text = s.text()
        check("Ctrl-C interrupts" in text, "banner lacks the Ctrl-C line")
        check("or CONTROL-C" not in text, "banner still says Ctrl-C exits")
        s.send("\x04")
        s.exit_code()
    finally:
        s.close()


@case
def ctrlc_infinite_loop(tmp):
    s = interactive(tmp)
    try:
        ctrlc_after(s, "Do[i, {i, 10^9}]", 0.5)
        then_still_works(s, 2)
    finally:
        s.close()


@case
def ctrlc_while_true(tmp):
    s = interactive(tmp)
    try:
        ctrlc_after(s, "While[True]", 0.5)
        then_still_works(s, 2)
    finally:
        s.close()


@case
def ctrlc_pause(tmp):
    s = interactive(tmp)
    try:
        dt = ctrlc_after(s, "Pause[30]", 0.3, timeout=1.0)
        check(dt < 1.0, "Pause took %.2fs to abort" % dt)
        then_still_works(s, 2)
    finally:
        s.close()


@case
def ctrlc_kernel(tmp):
    s = interactive(tmp)
    try:
        ctrlc_after(s, "FactorInteger[NextPrime[10^35] NextPrime[10^36]]", 0.5, timeout=3.0)
        then_still_works(s, 2)
    finally:
        s.close()


@case
def ctrlc_checkabort(tmp):
    s = interactive(tmp)
    try:
        s.send("CheckAbort[Do[i, {i, 10^9}], caught]\r")
        s.drain(0.5)
        s.send("\x03")
        s.expect(r"Out\[1\]= caught", 2.0)
        then_still_works(s, 2)
    finally:
        s.close()


@case
def ctrlc_abortprotect(tmp):
    s = interactive(tmp)
    try:
        s.send("AbortProtect[Pause[1]; x = 5]\r")
        s.drain(0.2)
        t0 = time.monotonic()
        s.send("\x03")
        s.expect(r"Out\[1\]= \$Aborted", 3.0)
        dt = time.monotonic() - t0
        check(dt >= 0.7, "AbortProtect body cut short (%.2fs)" % dt)
        wait_prompt(s, 2)
        s.send("x\r")
        s.expect(r"Out\[2\]= 5")
        then_still_works(s, 3)
    finally:
        s.close()


@case
def ctrlc_timeconstrained(tmp):
    s = interactive(tmp)
    try:
        s.send("TimeConstrained[Do[i, {i, 10^9}], 100, fail]\r")
        s.drain(0.5)
        s.send("\x03")
        s.expect(r"Out\[1\]= \$Aborted", 2.0)
        then_still_works(s, 2)
    finally:
        s.close()


@case
def ctrlc_prompt(tmp):
    s = interactive(tmp)
    try:
        s.send("1+")
        s.drain(0.2)
        s.send("\x03")
        s.expect(r"\^C")
        wait_prompt(s, 1)
        s.send("2+2\r")
        s.expect(r"Out\[1\]= 4")
        wait_prompt(s, 2)
        s.send("\x04")
        check(s.exit_code() == 0, "Ctrl-D did not exit 0")
    finally:
        s.close()


@case
def ctrlc_prompt_multiline(tmp):
    s = interactive(tmp)
    try:
        s.send("f[\r")
        s.drain(0.2)
        s.send("x")
        s.drain(0.2)
        s.send("\x03")
        wait_prompt(s, 1)
        s.send("3\r")
        s.expect(r"Out\[1\]= 3")
        wait_prompt(s, 2)
        s.send("\x04")
        check(s.exit_code() == 0, "Ctrl-D did not exit 0")
    finally:
        s.close()


@case
def ctrlc_dialog(tmp):
    s = interactive(tmp)
    try:
        s.send("Dialog[]\r")
        s.expect(r"In\[\d+\]:= ")
        s.send("Do[i, {i, 10^9}]\r")
        s.drain(0.5)
        s.send("\x03")
        s.expect(r"\$Aborted", 2.0)
        s.expect(r"In\[\d+\]:= ")
        s.send("Return[7]\r")
        s.expect(r"Out\[1\]= 7", 5.0)
        then_still_works(s, 2)
    finally:
        s.close()


@case
def double_ctrlc_exit(tmp):
    s = interactive(tmp)
    try:
        s.send("AbortProtect[Pause[20]]\r")
        s.drain(0.3)
        s.send("\x03")
        s.drain(0.3)
        s.send("\x03")
        code = s.exit_code(timeout=2.0)
        check(code == 130, "double Ctrl-C exit status %d, want 130" % code)
    finally:
        s.close()


HIST_HEADER = b"#mathilda-history v1\n"


def quit_session(s):
    s.send("\x04")
    check(s.exit_code() == 0, "Ctrl-D did not exit 0")


@case
def history_persist(tmp):
    hist = os.path.join(tmp, "h")
    extra = {"MATHILDA_HISTFILE": hist}
    s = interactive(tmp, extra=extra)
    try:
        s.send("a1 = 11\r")
        s.expect(r"Out\[1\]= 11")
        wait_prompt(s, 2)
        s.send("f[\r")
        s.drain(0.2)
        s.send("1]\r")
        s.expect(r"Out\[2\]= f\[1\]")
        wait_prompt(s, 3)
        s.send("a1 = 11\r")          # not consecutive: recorded again
        s.expect(r"Out\[3\]= 11")
        wait_prompt(s, 4)
        quit_session(s)
    finally:
        s.close()
    with open(hist, "rb") as f:
        data = f.read()
    want = HIST_HEADER + b"a1 = 11\nf[\\n1]\na1 = 11\n"
    check(data == want, "history file is %r, want %r" % (data, want))

    s = interactive(tmp, extra=extra)
    try:
        s.send("\x1b[A")             # Up: a1 = 11
        s.expect(r"a1 = 11")
        s.send("\x1b[A")             # Up again: the multi-line f[ / 1]
        s.expect(r"1\]")
        s.send("\r")
        s.expect(r"Out\[1\]= f\[1\]")
        wait_prompt(s, 2)
        quit_session(s)
    finally:
        s.close()


@case
def history_disabled(tmp):
    s = interactive(tmp, extra={"MATHILDA_HISTFILE": ""})
    try:
        s.send("1+1\r")
        s.expect(r"Out\[1\]= 2")
        wait_prompt(s, 2)
        s.send("\x1b[A")             # in-session history still works
        s.expect(r"1\+1")
        s.send("\r")
        s.expect(r"Out\[2\]= 2")
        wait_prompt(s, 3)
        quit_session(s)
    finally:
        s.close()
    check(os.listdir(tmp) == [], "files created with history disabled: %r" % os.listdir(tmp))


@case
def history_default_path(tmp):
    s = interactive(tmp)
    try:
        s.send("2+2\r")
        s.expect(r"Out\[1\]= 4")
        wait_prompt(s, 2)
        quit_session(s)
    finally:
        s.close()
    path = os.path.join(tmp, ".mathilda_history")
    with open(path, "rb") as f:
        check(f.read() == HIST_HEADER + b"2+2\n", "default history file content")
    check((os.stat(path).st_mode & 0o777) == 0o600, "history file is not 0600")


@case
def history_dialog(tmp):
    hist = os.path.join(tmp, "h")
    s = interactive(tmp, extra={"MATHILDA_HISTFILE": hist})
    try:
        s.send("Dialog[]\r")
        s.expect(r"In\[\d+\]:= ")
        s.send("inDialog = 3\r")
        s.expect(r"3")
        s.expect(r"In\[\d+\]:= ")
        s.send("Return[1]\r")
        s.expect(r"Out\[1\]= 1")
        wait_prompt(s, 2)
        quit_session(s)
    finally:
        s.close()
    with open(hist, "rb") as f:
        data = f.read()
    check(b"inDialog = 3\n" in data, "dialog line missing from history: %r" % data)


@case
def complete_symbol(tmp):
    s = interactive(tmp, extra={"MATHILDA_HISTFILE": ""})
    try:
        s.send("Integr\t")
        s.expect(r"Integrate")
        s.send("[x, x]\r")
        s.expect(r"Out\[1\]= 1/2 x\^2")
        wait_prompt(s, 2)
        quit_session(s)
    finally:
        s.close()


@case
def complete_option(tmp):
    s = interactive(tmp, extra={"MATHILDA_HISTFILE": ""})
    try:
        s.send("ReplaceRepeated[x, {}, MaxI\t")
        s.expect(r"MaxIterations -> ")
        s.send("5]\r")
        s.expect(r"Out\[1\]= x")
        check(b"ReplaceRepeated[x, {}, MaxIterations -> 5]" in s.buf,
              "the option was not inserted whole")
        wait_prompt(s, 2)
        quit_session(s)
    finally:
        s.close()


@case
def complete_user_option(tmp):
    s = interactive(tmp, extra={"MATHILDA_HISTFILE": ""})
    try:
        s.send("Options[myf] = {Alpha -> 1};\r")
        wait_prompt(s, 2)
        s.send("myf[Al\t\t")           # ambiguous (All, ...): the second Tab lists
        s.expect(r"Alpha -> ")
        s.send("pha -> 2]\r")
        s.expect(r"Out\[2\]= myf\[Alpha -> 2\]")
        wait_prompt(s, 3)
        quit_session(s)
    finally:
        s.close()


@case
def complete_no_eval(tmp):
    s = interactive(tmp, extra={"MATHILDA_HISTFILE": ""})
    try:
        s.send("f[x_] := (Print[\"SIDE\"]; x)\r")
        wait_prompt(s, 2)
        start = len(s.buf)
        s.send("f[1][Op\t")
        s.drain(0.5)
        check(b"SIDE" not in s.buf[start:], "completion evaluated f[1]")
        s.send("\x03")               # abandon the line
        wait_prompt(s, 2)
        quit_session(s)
    finally:
        s.close()


# Colour is SGR (ESC [ ... m). The in-tree editor also emits cursor movement,
# erase and position-report sequences, which are not colour.
SGR = re.compile(rb"\x1b\[[0-9;]*m")


def no_sgr(data):
    return SGR.search(data) is None


def colour_session(tmp, **kw):
    s = interactive(tmp, **kw)
    s.send("1+1\r")
    s.expect(r"Out\[1\]=.*2")
    s.send("1/0\r")
    s.expect(r"ComplexInfinity")
    s.drain(0.3)
    quit_session(s)
    return s


@case
def color_on(tmp):
    s = colour_session(tmp, extra={"MATHILDA_HISTFILE": ""}, drop=("NO_COLOR", "CLICOLOR_FORCE"))
    try:
        b = s.buf
        check(b"\x1b[1;34m" in b, "no coloured In prompt")
        check(b"\x1b[1;35mOut[1]=\x1b[0m" in b, "no coloured Out")
        check(b"\x1b[1;31mPower::infy" in b, "no coloured message head")
    finally:
        s.close()


@case
def color_no_color(tmp):
    s = colour_session(tmp, extra={"MATHILDA_HISTFILE": "", "NO_COLOR": "1"})
    try:
        check(no_sgr(s.buf), "SGR with NO_COLOR=1")
    finally:
        s.close()


@case
def color_dumb(tmp):
    s = colour_session(tmp, term="dumb", extra={"MATHILDA_HISTFILE": ""}, drop=("NO_COLOR",))
    try:
        check(no_sgr(s.buf), "SGR with TERM=dumb")
    finally:
        s.close()


def run_stdout_pipe(tmp, extra):
    """stdin and stderr on a pty, stdout a pipe."""
    master, slave = pty.openpty()
    fcntl.ioctl(slave, termios.TIOCSWINSZ, struct.pack("HHHH", ROWS, COLS, 0, 0))
    env = base_env(tmp, extra=dict({"MATHILDA_HISTFILE": ""}, **extra))
    p = subprocess.Popen(argv([]), stdin=slave, stdout=subprocess.PIPE,
                         stderr=slave, env=env)
    os.close(slave)
    chunks = []

    def pump():
        while True:
            try:
                c = os.read(master, 65536)
            except OSError:
                return
            if not c:
                return
            chunks.append(c)

    t = threading.Thread(target=pump, daemon=True)
    t.start()
    time.sleep(1.0 * TIME_SCALE)
    os.write(master, b"1+1\r")
    time.sleep(0.5 * TIME_SCALE)
    os.write(master, b"\x04")
    out = p.communicate(timeout=10 * TIME_SCALE)[0]
    t.join(2)
    os.close(master)
    return out


@case
def color_stdout_pipe(tmp):
    out = run_stdout_pipe(tmp, {})
    check(b"Out[1]= 2" in out, "no result on the stdout pipe: %r" % out[-300:])
    check(no_sgr(out), "SGR on a piped stdout")


@case
def color_force(tmp):
    out = run_stdout_pipe(tmp, {"CLICOLOR_FORCE": "1"})
    check(b"\x1b[1;35mOut[1]=\x1b[0m 2" in out, "CLICOLOR_FORCE did not colour: %r" % out[-300:])


@case
def took_line(tmp):
    s = interactive(tmp, extra={"MATHILDA_HISTFILE": ""})
    try:
        s.send("Pause[1.2]\r")
        s.expect(r"Out\[1\]=")
        s.expect(r"took 1\.2\d s")
        wait_prompt(s, 2)
        start = len(s.buf)
        s.send("1+1\r")
        s.expect(r"Out\[2\]=.*2")
        wait_prompt(s, 3)
        check(b"took" not in s.buf[start:], "took shown for a fast result")
        quit_session(s)
    finally:
        s.close()


# ---------------------------------------------------------------------------
# The in-tree editor (REPL PR 2). It is the default on a terminal; the cases
# above already run through it. These cover what only the editor does.
# ---------------------------------------------------------------------------

NOHIST = {"MATHILDA_HISTFILE": ""}


@case
def editor_is_default(tmp):
    s = interactive(tmp, extra=NOHIST)
    try:
        check("Tab completes," in s.text() and "Type / for commands" in s.text(),
              "the editor banner lines are missing")
        check(b"\x1b[?2004h" in s.buf, "bracketed paste not enabled")
        quit_session(s)
    finally:
        s.close()


@case
def editor_multiline_updown(tmp):
    """Up moves between the buffer's lines; editing there keeps the rest."""
    s = interactive(tmp, extra=NOHIST)
    try:
        s.send("f[\r")
        s.drain(0.2)
        s.send("2]")
        s.drain(0.2)
        s.send("\x1b[A")            # Up: into the first line, not history
        s.send("\x05")              # Ctrl-E
        s.send("1,")
        s.drain(0.2)
        s.send("\r")
        s.expect(r"Out\[1\]= f\[1, 2\]")
        wait_prompt(s, 2)
        quit_session(s)
    finally:
        s.close()


@case
def editor_force_return(tmp):
    """Esc then Return evaluates an incomplete buffer."""
    s = interactive(tmp, extra=NOHIST)
    try:
        s.send("f[")
        s.drain(0.2)
        s.send("\x1b\r")
        wait_prompt(s, 2, timeout=10)
        check(b"Out[1]" not in s.buf, "an incomplete input produced a result")
        s.send("1+1\r")
        s.expect(r"Out\[2\]= 2")
        quit_session(s)
    finally:
        s.close()


@case
def editor_emacs_keys(tmp):
    s = interactive(tmp, extra=NOHIST)
    try:
        s.send("+1\x01")             # Ctrl-A: to the start
        s.send("1")
        s.send("\r")
        s.expect(r"Out\[1\]= 2")
        wait_prompt(s, 2)
        s.send("x + junk\x17\x17")  # Ctrl-W twice: "x "
        s.send("40\r")
        s.expect(r"Out\[2\]= 40 x")
        wait_prompt(s, 3)
        s.send("2+2\x1f")           # Ctrl-_ undoes the typing
        s.send("3+3\r")
        s.expect(r"Out\[3\]= 6")
        quit_session(s)
    finally:
        s.close()


@case
def editor_paste(tmp):
    """A pasted Return does not submit; the pasted block is one input."""
    s = interactive(tmp, extra=NOHIST)
    try:
        s.send("\x1b[200~f[\r1]\x1b[201~")
        s.drain(0.4)
        check(b"Out[1]" not in s.buf, "the paste was submitted")
        s.send("\r")
        s.expect(r"Out\[1\]= f\[1\]")
        quit_session(s)
    finally:
        s.close()


@case
def editor_typeahead(tmp):
    """Keys typed during an evaluation become the next input."""
    s = interactive(tmp, extra=NOHIST)
    try:
        s.send("Pause[0.6]\r")
        s.drain(0.1)
        s.send("1+1\r")
        s.expect(r"Out\[2\]= 2", 10)
        quit_session(s)
    finally:
        s.close()


@case
def editor_utf8(tmp):
    s = interactive(tmp, extra=NOHIST)
    try:
        s.send("\"\u4e2d\" <> \"a\"\r".encode("utf-8"))
        s.expect('Out\\[1\\]= "\u4e2da"'.encode("utf-8"))
        quit_session(s)
    finally:
        s.close()


@case
def editor_resize(tmp):
    """A resize redraws the buffer at the new width."""
    s = interactive(tmp, extra=NOHIST)
    try:
        s.send("abc")
        s.drain(0.2)
        mark = len(s.buf)
        fcntl.ioctl(s.fd, termios.TIOCSWINSZ, struct.pack("HHHH", 20, 60, 0, 0))
        try:
            os.kill(s.pid, signal.SIGWINCH)
        except OSError:
            pass
        s.drain(0.4)
        check(b"In[1]:= abc" in s.buf[mark:], "no redraw after the resize: %r" % s.buf[mark:])
        s.send("\r")
        s.expect(r"Out\[1\]= abc")
        quit_session(s)
    finally:
        s.close()


@case
def popup_select_accept(tmp):
    """Tab opens the list with a description; Down and Enter accept."""
    s = interactive(tmp, extra=NOHIST)
    try:
        s.send("Sin\t")
        s.expect(r"SinIntegral")
        s.expect(r"Sin\[z\]")                 # the description row
        s.send("\x1b[B")                       # Down: SinIntegral
        s.drain(0.2)
        s.send("\r")                           # Enter accepts, does not submit
        s.drain(0.3)
        check(b"Out[1]" not in s.buf, "Enter in the popup submitted the line")
        s.send("[1.]\r")
        s.expect(r"Out\[1\]= 0\.946")       # SinIntegral[1.], not Sin[1.]
        quit_session(s)
    finally:
        s.close()


@case
def popup_options(tmp):
    """Inside a call, option names come first."""
    s = interactive(tmp, extra=NOHIST)
    try:
        s.send("NIntegrate[x, {x, 0, 1}, \t")
        s.expect(r"AccuracyGoal -> ")
        s.send("\x1b")
        s.drain(0.2)
        s.send("WorkingPrec\t")               # unique: inserted at once
        s.expect(r"WorkingPrecision -> ")
        s.send("20]\r")
        s.expect(r"Out\[1\]= 0\.5")
        quit_session(s)
    finally:
        s.close()


@case
def ctrl_r_search(tmp):
    """Ctrl-R lists matching inputs, multi-line ones whole; Enter edits it."""
    hist = os.path.join(tmp, "h")
    with open(hist, "wb") as f:
        f.write(HIST_HEADER + b"alpha = 1\nf[\\n  2]\nbeta = 3\n")
    s = interactive(tmp, extra={"MATHILDA_HISTFILE": hist})
    try:
        s.send("\x12")
        s.expect(r"search: ")
        s.expect(r"beta = 3")
        s.expect(r"      2\]")                # a multi-line entry, shown whole
        s.send("f[")                           # filter: only the f[ entry is left
        s.drain(0.3)
        check(b"(no match)" not in s.buf, "the filter lost the entry")
        s.send("\r")
        s.drain(0.3)
        check(b"Out[1]" not in s.buf, "Ctrl-R's Enter submitted the entry")
        s.send("\r")
        s.expect(r"Out\[1\]= f\[2\]")
        wait_prompt(s, 2)
        s.send("\x12zzz")
        s.expect(r"\(no match\)")
        s.send("\x07")                         # Ctrl-G: cancel
        s.send("7\r")
        s.expect(r"Out\[2\]= 7")
        quit_session(s)
    finally:
        s.close()


RESULT_INPUTS = ["Sin[x] + 2/3 y^2", "{\"str (* x\", 1.5, f[_, x__]}", "Expand[(a+b)^3]",
                 "Format[fmt[x_]] := \"F\"[x]; fmt[1]", "x -> 2"]


def result_lines(tmp, colour):
    extra = dict(NOHIST)
    drop = ("NO_COLOR", "CLICOLOR_FORCE") if colour else ()
    s = new_session(tmp, extra=extra, drop=drop) if colour else new_session(tmp, extra=extra)
    try:
        wait_prompt(s, 1)
        outs = []
        for i, inp in enumerate(RESULT_INPUTS):
            s.send(inp + "\r")
            m = s.expect(r"Out\[%d\]=(?:\x1b\[0m)? (.*?)\r?\n" % (i + 1))
            outs.append(m.group(1))
            wait_prompt(s, i + 2)
        quit_session(s)
        return outs
    finally:
        s.close()


@case
def result_colour(tmp):
    """Results are coloured in the editor session, with the same characters."""
    coloured = result_lines(tmp, True)
    plain = result_lines(tmp, False)
    check(b"\x1b[36mSin" in coloured[0], "no builtin colour in %r" % coloured[0])
    check(b"\x1b[35m" in coloured[0], "no number colour in %r" % coloured[0])
    for c, p in zip(coloured, plain):
        check(no_sgr(p), "SGR in a NO_COLOR result: %r" % p)
        check(SGR.sub(b"", c) == p, "coloured %r differs from plain %r" % (c, p))


@case
def editor_narrow(tmp):
    """On a 14-column terminal the popup gives way to inline cycling."""
    s = new_session(tmp, extra=NOHIST, cols=14)
    try:
        wait_prompt(s, 1)
        s.send("Sin\t\t")                     # Sin, then SinIntegral
        s.drain(0.3)
        s.send("[1.]\r")
        s.expect(r"0\.946")
        quit_session(s)
    finally:
        s.close()


@case
def editor_latency_e2e(tmp):
    """Key written to the pty -> the editor's redraw read back, with a
    200-line buffer. Reported, not gated: the scheduler of a loaded machine
    dominates it (tests/test_term_editor_latency.c is the gate)."""
    s = interactive(tmp, extra=NOHIST)
    try:
        lines = ",\n".join("f%d[x_] := Sin[x] + %d x^2 (* line %d *)" % (i, i, i) for i in range(200))
        data = ("\x1b[200~{" + lines + "\x1b[201~").encode()
        for i in range(0, len(data), 512):     # in chunks: the redraw must be read
            s.send(data[i:i + 512])            # meanwhile, or both sides block
            s.drain(0.02)
        s.drain(0.5)
        for _ in range(100):
            s.send("\x1b[A")
        s.drain(0.5)
        times = []
        for i in range(200):
            key = b"a" if i % 2 == 0 else b"\x7f"
            mark = len(s.buf)
            t0 = time.monotonic()
            s.send(key)
            deadline = t0 + 2.0
            while len(s.buf) == mark and time.monotonic() < deadline:
                s._read_some(0.05)
            times.append((time.monotonic() - t0) * 1000.0)
            s.drain(0.002)
        times.sort()
        print("    e2e keystroke->redraw (ms): p50 %.2f  p95 %.2f  max %.2f"
              % (times[len(times) // 2], times[len(times) * 95 // 100], times[-1]))
        s.send("\x03")
        wait_prompt(s, 1)
        quit_session(s)
    finally:
        s.close()


def tty_flags(fd):
    attrs = termios.tcgetattr(fd)
    return attrs[3]               # lflag


@case
def editor_restores_terminal(tmp):
    """Raw mode while reading; cooked mode while evaluating and after a SIGTERM."""
    s = interactive(tmp, extra=NOHIST)
    try:
        lflag = tty_flags(s.fd)
        check(not (lflag & termios.ICANON) and not (lflag & termios.ISIG),
              "not in raw mode at the prompt")
        s.send("Pause[1]\r")
        s.drain(0.4)
        lflag = tty_flags(s.fd)
        # cooked during the evaluation; ECHO is off because the status line
        # owns the terminal then (it echoes ^C itself)
        check((lflag & termios.ICANON) and (lflag & termios.ISIG),
              "cooked mode not restored during the evaluation")
        check(not (lflag & termios.ECHO), "ECHO left on under the status line")
        wait_prompt(s, 2)
        lflag = tty_flags(s.fd)
        check(not (lflag & termios.ICANON) and not (lflag & termios.ISIG),
              "not back in raw mode at the next prompt")
        os.kill(s.pid, signal.SIGTERM)
        st = s.finish(5)
        check(os.WIFSIGNALED(st) and os.WTERMSIG(st) == signal.SIGTERM,
              "SIGTERM did not end the process by the signal: %r" % st)
        lflag = tty_flags(s.fd)
        check((lflag & termios.ICANON) and (lflag & termios.ECHO),
              "terminal left in raw mode after SIGTERM")
        check(s.buf.endswith(b"\x1b[?2004l") or b"\x1b[?2004l" in s.buf[-64:],
              "bracketed paste left on after SIGTERM")
    finally:
        s.close()


@case
def editor_readline_fallback(tmp):
    """MATHILDA_EDITOR=readline gives PR 1's readline session."""
    s = interactive(tmp, extra=dict(NOHIST, MATHILDA_EDITOR="readline"))
    try:
        check("Type / for commands" not in s.text(), "editor banner with the readline fallback")
        check(b"\x1b[6n" not in s.buf, "the editor ran despite MATHILDA_EDITOR=readline")
        s.send("f[\r")
        s.drain(0.2)
        s.send("1]\r")
        s.expect(r"Out\[1\]= f\[1\]")
        wait_prompt(s, 2)
        s.send("Integr\t")
        s.expect(r"Integrate")
        s.send("[x, x]\r")
        s.expect(r"Out\[2\]= 1/2 x\^2")
        wait_prompt(s, 3)
        s.send("1+")
        s.drain(0.2)
        s.send("\x03")
        s.expect(r"\^C")
        wait_prompt(s, 3)
        quit_session(s)
    finally:
        s.close()


# ---------------------------------------------------------------------------
# REPL PR 3: the live status line (spinner), progress, took note, status line
# ---------------------------------------------------------------------------

BRAILLE = re.compile(rb"\xe2\xa0[\x80-\xbf]")
EL = b"\x1b[2K"


def spinner_session(tmp, colour=True, **kw):
    extra = dict({"MATHILDA_HISTFILE": ""}, **(kw.pop("extra", None) or {}))
    drop = tuple(kw.pop("drop", ())) + (("NO_COLOR",) if colour else ())
    if not colour:
        extra.setdefault("NO_COLOR", "1")
    return interactive(tmp, extra=extra, drop=drop, **kw)


def segment(s, start):
    """The transcript from `start` to the current match position."""
    return s.buf[start:s.pos]


def frames(data):
    """The spinner frames drawn in `data` (colour mode): text after EL up to
    the next \r."""
    out = []
    for chunk in data.split(EL)[1:]:
        f = chunk.split(b"\r")[0]
        if BRAILLE.search(f):
            out.append(f)
    return out


def strip_sgr(b):
    return re.sub(rb"\x1b\[[0-9;]*m", b"", b)


@case
def spinner_shows(tmp):
    s = spinner_session(tmp)
    try:
        start = len(s.buf)
        s.send("Pause[1]\r")
        s.expect(r"Out\[1\]=")
        seg = segment(s, start)
        fs = frames(seg)
        check(len(fs) >= 5, "expected ~8 frames in 1 s, got %d: %r" % (len(fs), seg[-400:]))
        plain = strip_sgr(fs[-1])
        check(b"Pause" in plain and b"cpu" in plain and b"mem" in plain,
              "frame lacks label/cpu/mem: %r" % plain)
        # the line is erased before Out[1]
        check(seg.rstrip(b"\x1b[1;35mOut[1]=").endswith(b"\r" + EL) or
              (b"\r" + EL + b"\x1b[1;35mOut[1]=") in s.buf, "spinner not erased before Out")
        s.expect(r"took 1\.\d\d s")
        wait_prompt(s, 2)
        quit_session(s)
    finally:
        s.close()


@case
def spinner_quick_no_flicker(tmp):
    s = spinner_session(tmp)
    try:
        start = len(s.buf)
        for n, inp in enumerate(("1+1", 'Print["hi"]; 2', "Table[i^2, {i, 1000}];",
                                 "Expand[(1 + x)^20];"), 1):
            s.send(inp + "\r")
            s.expect(r"Out\[%d\]=" % n)
            wait_prompt(s, n + 1)
        seg = s.buf[start:]
        check(not BRAILLE.search(seg) and EL not in seg,
              "a quick result drew the spinner: %r" % seg[-500:])
        check(b"hi\r\n" in seg, "Print output missing")
        quit_session(s)
    finally:
        s.close()


@case
def spinner_print_coexist(tmp):
    s = spinner_session(tmp)
    try:
        s.send('Do[Print["line", i]; Pause[0.15], {i, 8}]\r')
        s.expect(r"Out\[1\]=")
        seg = s.buf
        # every printed line arrives whole, at the start of a line, after an erase
        # of the spinner (or right after the input)
        for i in range(1, 9):
            m = re.search(rb"(\r\n|\r\x1b\[2K|\r|\x1b\[\?2004l)line%d\r\n" % i, seg)
            check(m, "line%d not printed cleanly: %r" % (i, seg[-600:]))
        # no frame is glued to printed text
        check(not re.search(rb"line\d[^\r]*\xe2\xa0", seg), "spinner overwrote output")
        check(len(frames(seg)) >= 3, "no spinner between prints")
        quit_session(s)
    finally:
        s.close()


@case
def spinner_stderr_order(tmp):
    s = spinner_session(tmp, colour=False)
    try:
        s.send('f::m = "err `1`"; Do[Print["out", i]; Message[f::m, i]; Pause[0.1], {i, 5}]\r')
        s.expect(r"Out\[1\]=")
        text = s.buf.decode("utf-8", "replace")
        order = re.findall(r"out\d|err \d", text)
        want = [x for i in range(1, 6) for x in ("out%d" % i, "err %d" % i)]
        check(order == want, "stdout/stderr order %r, want %r" % (order, want))
        quit_session(s)
    finally:
        s.close()


@case
def spinner_progress(tmp):
    s = spinner_session(tmp, colour=False)
    try:
        s.send("Do[Pause[0.04], {i, 40}]\r")
        s.expect(r"Out\[1\]=", 10.0)
        m = re.findall(rb"Do (\d+)/40 \xe2\x96\x95[^\r]*?\xe2\x96\x8f +(\d+)%", s.buf)
        check(len(m) >= 5, "no progress bar frames: %r" % s.buf[-500:])
        done = [int(a) for a, _ in m]
        check(done == sorted(done) and done[-1] > done[0], "progress not advancing: %r" % done)
        # nested: the outer loop and the one inside it
        s.send("Table[Pause[0.02]; i j, {i, 3}, {j, 10}];\r")
        s.expect(r"Out\[2\]=", 10.0)
        check(re.search(rb"Table \d/3 \xe2\x96\xb8 Table \d+/10", s.buf), "no nested progress")
        quit_session(s)
    finally:
        s.close()


@case
def spinner_ctrlc(tmp):
    s = spinner_session(tmp)
    try:
        s.send("Do[i, {i, 10^9}]\r")
        s.drain(0.7)
        check(frames(s.buf), "spinner not running before Ctrl-C")
        t0 = time.monotonic()
        s.send("\x03")
        s.expect(r"Out\[1\]=\S* (\x1b\[[0-9;]*m)*\$Aborted", 2.0)
        check(time.monotonic() - t0 < 1.0, "abort took too long")
        # the transcript is PR 1's: ^C, newline, then the result -- the spinner
        # line erased first and no stale frame left behind
        check(re.search(rb"\r\x1b\[2K\^C\r\n\x1b\[1;35mOut\[1\]=", s.buf),
              "unexpected Ctrl-C transcript: %r" % s.buf[-300:])
        s.send("1+1\r")
        s.expect(r"Out\[2\]=\S* (\x1b\[[0-9;]*m)*2")
        wait_prompt(s, 3)
        quit_session(s)
    finally:
        s.close()


@case
def spinner_double_ctrlc(tmp):
    s = spinner_session(tmp, colour=False)
    try:
        s.send("AbortProtect[Pause[20]]\r")
        s.drain(0.5)
        s.send("\x03")
        s.drain(0.3)
        s.send("\x03")
        code = s.exit_code(timeout=2.0)
        check(code == 130, "double Ctrl-C exit status %d, want 130" % code)
        tail_b = s.buf[-200:]
        check(tail_b.endswith(b"^C\r\n"), "did not end with ^C newline: %r" % tail_b)
        # each Ctrl-C shows as "^C" + newline (the tty echo is off under the
        # spinner), and no frame is left on screen: anything drawn after the
        # first one was blanked (\r spaces \r) before the second
        first = s.buf.find(b"^C\r\n")
        check(first > 0 and s.buf.count(b"^C\r\n") == 2, "want two ^C lines: %r" % tail_b)
        after = s.buf[first + 4:-4]
        if BRAILLE.search(after):
            check(re.search(rb"\r +\r$", after), "stale frame left: %r" % tail_b)
    finally:
        s.close()


@case
def spinner_no_color_plain(tmp):
    s = spinner_session(tmp, colour=False)
    try:
        s.send("Pause[0.6]\r")
        s.expect(r"Out\[1\]=")
        wait_prompt(s, 2)
        quit_session(s)
        check(BRAILLE.search(s.buf), "no plain spinner with NO_COLOR")
        check(no_sgr(s.buf), "escape sequence with NO_COLOR")
    finally:
        s.close()


@case
def spinner_dumb_off(tmp):
    s = spinner_session(tmp, colour=False, term="dumb")
    try:
        s.send("Pause[0.6]\r")
        s.expect(r"Out\[1\]= Null")
        wait_prompt(s, 2)
        quit_session(s)
        check(not BRAILLE.search(s.buf) and b"\r\r" not in s.buf.split(b"Pause[0.6]")[1].split(b"Out[1]")[0],
              "spinner with TERM=dumb: %r" % s.buf[-300:])
    finally:
        s.close()


@case
def spinner_disabled_env(tmp):
    s = spinner_session(tmp, extra={"MATHILDA_SPINNER": "0"})
    try:
        s.send("Pause[0.6]\r")
        s.expect(r"Out\[1\]=")
        wait_prompt(s, 2)
        quit_session(s)
        check(not BRAILLE.search(s.buf) and EL not in s.buf, "spinner with MATHILDA_SPINNER=0")
    finally:
        s.close()


@case
def spinner_ascii_locale(tmp):
    s = spinner_session(tmp, colour=False, extra={"LANG": "C", "LC_ALL": "C", "LC_CTYPE": "C"})
    try:
        s.send("Pause[0.6]\r")
        s.expect(r"Out\[1\]=")
        wait_prompt(s, 2)
        quit_session(s)
        check(not BRAILLE.search(s.buf), "braille in a C locale")
        check(re.search(rb"\r[|/\\-] Pause \| ", s.buf), "no ASCII spinner: %r" % s.buf[-300:])
    finally:
        s.close()


@case
def spinner_narrow(tmp):
    s = spinner_session(tmp, colour=False, cols=30)
    try:
        s.send("Do[Pause[0.05], {i, 20}]\r")
        s.expect(r"Out\[1\]=", 10.0)
        lines = [l for l in s.buf.split(b"\r") if BRAILLE.search(l)]
        check(lines, "no frames at 30 columns")
        for l in lines:
            w = len(l.decode("utf-8").rstrip(" "))
            check(w <= 29, "frame %d cells wide at 30 columns: %r" % (w, l))
        quit_session(s)
    finally:
        s.close()


@case
def spinner_dialog(tmp):
    s = spinner_session(tmp, colour=False)
    try:
        s.send("Dialog[]\r")
        s.expect(r"In\[\d+\]:= ")
        start = len(s.buf)
        s.send("Pause[0.5]; 5\r")
        s.expect(r"Out\[\d+\]= 5")
        check(BRAILLE.search(s.buf[start:]), "no spinner inside Dialog")
        s.expect(r"In\[\d+\]:= ")
        s.send("Return[7]\r")
        s.expect(r"Out\[1\]= 7", 5.0)
        then_still_works(s, 2)
    finally:
        s.close()


@case
def statusline_env(tmp):
    s = spinner_session(tmp, colour=False)
    try:
        s.send("1+1\r")
        s.expect(r"Out\[1\]= 2")
        wait_prompt(s, 2)
        check(b"last " not in s.buf, "status line shown by default")
        quit_session(s)
    finally:
        s.close()
    s = spinner_session(tmp, colour=False, extra={"MATHILDA_STATUSLINE": "1"})
    try:
        s.send("1+1\r")
        s.expect(r"Out\[1\]= 2")
        s.expect(r"last \d+\.\d\d s \xc2\xb7 mem \d+ [KMG]B \xc2\xb7 Global`\r\n")
        wait_prompt(s, 2)
        quit_session(s)
    finally:
        s.close()


@case
def took_line_cpu(tmp):
    s = spinner_session(tmp, colour=False)
    try:
        s.send("Do[Null, {2*10^6}]; Pause[1]\r")
        s.expect(r"Out\[1\]=")
        s.expect(r"took 1\.\d\d s \xc2\xb7 cpu 0\.\d\d s( \xc2\xb7 peak \d+ [KMG]B)?\r\n")
        quit_session(s)
    finally:
        s.close()


# ---------------------------------------------------------------------------
# The editor and the status line together: the editor owns the terminal while
# reading, the status thread while evaluating, and each hands it back cleanly.
# ---------------------------------------------------------------------------

def is_raw(fd):
    lf = tty_flags(fd)
    return not (lf & termios.ICANON) and not (lf & termios.ISIG)


@case
def combo_multiline_then_spinner(tmp):
    """A multi-line edit, then a slow evaluation: frames drawn, erased before
    Out, and the next prompt is the editor again (raw, history recall)."""
    s = spinner_session(tmp)
    try:
        s.send("Do[Pause[0.1],\r")
        s.drain(0.2)
        s.send("{i, 6}]")
        s.drain(0.1)
        start = len(s.buf)
        s.send("\r")
        s.expect(r"Out\[1\]=")
        seg = segment(s, start)
        check(len(frames(seg)) >= 3, "no spinner after a multi-line edit: %r" % seg[-300:])
        check(b"\r\x1b[2K\x1b[1;35mOut[1]=" in seg, "spinner not erased before Out")
        wait_prompt(s, 2)
        s.drain(0.1)
        check(is_raw(s.fd), "the editor did not get the terminal back")
        s.send("\x1b[A")                 # Up on an empty buffer: history
        s.drain(0.2)
        check(b"{i, 6}]" in strip_sgr(s.buf[s.pos:]), "history recall after the spinner failed")
        s.send("\x03")                   # Ctrl-C at the prompt clears the buffer
        s.send("1+1\r")
        s.expect(r"Out\[2\]=\S* (\x1b\[[0-9;]*m)*2")
        quit_session(s)
    finally:
        s.close()


@case
def combo_print_after_edit(tmp):
    """Print under the spinner, after an edited input, arrives whole and in
    order with the messages."""
    s = spinner_session(tmp, colour=False)
    try:
        s.send('Do[Print["p", i]; Message[General::ovfl]; Pause[0.12],\r')
        s.drain(0.2)
        s.send("{i, 4}]")
        s.drain(0.2)
        start = len(s.buf)
        s.send("\r")
        s.expect(r"Out\[1\]=", 10.0)
        text = s.buf[start:].decode("utf-8", "replace")
        order = re.findall(r"p\d(?=\r\n)|General::ovfl(?=:)", text)
        want = [x for i in range(1, 5) for x in ("p%d" % i, "General::ovfl")]
        check(order == want, "order %r, want %r" % (order, want))
        for i in range(1, 5):
            check(re.search(rb"(\r\n|\r|\x1b\[\?2004l)p%d\r\n" % i, s.buf), "p%d not on its own line" % i)
        wait_prompt(s, 2)
        quit_session(s)
    finally:
        s.close()


@case
def combo_dialog_editor(tmp):
    """Dialog[] under the spinner reads through the editor: the spinner is
    suspended (no frames while the dialog prompt waits, raw mode there) and
    resumed for the dialog's own evaluations."""
    s = spinner_session(tmp, colour=False, extra={"MATHILDA_SPINNER_DELAY_MS": "0"})
    try:
        s.send("Pause[0.3]; Dialog[]\r")
        s.expect(r"In\[\d+\]:= ", 5.0)
        quiet = len(s.buf)
        s.drain(0.5)
        check(not BRAILLE.search(s.buf[quiet:]), "spinner drew over the dialog prompt")
        check(is_raw(s.fd), "the dialog prompt is not the editor")
        s.send("f[\r")                  # the dialog's input is multi-line too
        s.drain(0.2)
        s.send("Pause[0.4]; 5]\r")
        s.expect(r"Out\[\d+\]= f\[5\]", 5.0)
        check(BRAILLE.search(s.buf[quiet:]), "no spinner for the dialog's evaluation")
        s.expect(r"In\[\d+\]:= ")
        s.send("Return[7]\r")
        s.expect(r"Out\[1\]= 7", 5.0)
        then_still_works(s, 2)
    finally:
        s.close()


@case
def combo_ctrlc_multiline(tmp):
    """Ctrl-C under the spinner after a multi-line input: $Aborted, the
    PR 1 transcript, the editor back in raw mode with an empty buffer."""
    s = spinner_session(tmp)
    try:
        s.send("Do[i,\r")
        s.drain(0.2)
        s.send("{i, 10^9}]\r")
        s.drain(0.7)
        check(frames(s.buf), "spinner not running before Ctrl-C")
        s.send("\x03")
        s.expect(r"Out\[1\]=\S* (\x1b\[[0-9;]*m)*\$Aborted", 2.0)
        check(re.search(rb"\r\x1b\[2K\^C\r\n\x1b\[1;35mOut\[1\]=", s.buf),
              "unexpected Ctrl-C transcript: %r" % s.buf[-300:])
        wait_prompt(s, 2)
        s.drain(0.1)
        check(is_raw(s.fd), "editor not in raw mode after the abort")
        s.send("2+2\r")
        s.expect(r"Out\[2\]=\S* (\x1b\[[0-9;]*m)*4")
        quit_session(s)
    finally:
        s.close()


@case
def combo_double_ctrlc_restores(tmp):
    """A double Ctrl-C exit under the spinner leaves a cooked, echoing
    terminal with bracketed paste off."""
    s = spinner_session(tmp, colour=False)
    try:
        s.send("AbortProtect[Pause[20]]\r")
        s.drain(0.5)
        s.send("\x03")
        s.drain(0.3)
        s.send("\x03")
        code = s.exit_code(timeout=2.0)
        check(code == 130, "exit status %d, want 130" % code)
        lf = tty_flags(s.fd)
        check((lf & termios.ICANON) and (lf & termios.ECHO) and (lf & termios.ISIG),
              "terminal not restored after the double Ctrl-C")
    finally:
        s.close()


@case
def combo_sigterm_under_spinner(tmp):
    """SIGTERM during a spinner evaluation: terminal restored (ECHO back)."""
    s = spinner_session(tmp, colour=False)
    try:
        s.send("Pause[5]\r")
        s.drain(0.6)
        check(not (tty_flags(s.fd) & termios.ECHO), "ECHO not off under the spinner")
        os.kill(s.pid, signal.SIGTERM)
        st = s.finish(5)
        check(os.WIFSIGNALED(st) and os.WTERMSIG(st) == signal.SIGTERM,
              "SIGTERM did not end the process by the signal: %r" % st)
        lf = tty_flags(s.fd)
        check((lf & termios.ICANON) and (lf & termios.ECHO), "terminal not restored")
    finally:
        s.close()


# ---------------------------------------------------------------------------
# REPL PR 4: slash commands, themes, mods, Monitor
# ---------------------------------------------------------------------------

def strip_esc(b):
    """The transcript without escape sequences."""
    return re.sub(rb"\x1b\[[0-9;?]*[A-Za-z]", b"", b)


@case
def slash_menu(tmp):
    """'/' on an empty buffer opens the menu; Return runs the selection; no
    In[n] is used and nothing is evaluated."""
    s = interactive(tmp, extra=NOHIST)
    try:
        s.send("/")
        s.expect(r">/help +list the commands")
        s.expect(r"/theme +choose a colour theme")
        s.send("hi")                              # filter: /history
        s.drain(0.2)
        s.send("\r")
        s.expect(r" 1  /history\r\n")              # commands are history entries
        wait_prompt(s, 1)                         # the same In[1] again
        s.send("1+1\r")
        s.expect(r"Out\[1\]= 2")
        check(b"Out[0]" not in s.buf and b"$Failed" not in s.buf, "a command was evaluated")
        s.send("/help\r")
        s.expect(r"Commands\r\n  /help ")
        s.expect(r"Keys\r\n")
        wait_prompt(s, 2)
        s.send("/nope\r")
        s.expect(r"Unknown command /nope")
        wait_prompt(s, 2)
        s.send("%\r")                             # Out[1], not a command result
        s.expect(r"Out\[2\]= 2")
        quit_session(s)
    finally:
        s.close()


@case
def slash_ordinary_operators(tmp):
    """'/' mid-expression, '/.', '//.' and '//' stay operators; no menu."""
    s = interactive(tmp, extra=NOHIST)
    try:
        start = len(s.buf)
        s.send("x^2 /. x -> 3\r")
        s.expect(r"Out\[1\]= 9")
        s.send("{a, b} //. a -> b\r")
        s.expect(r"Out\[2\]= \{b, b\}")
        s.send("16 // Sqrt\r")
        s.expect(r"Out\[3\]= 4")
        s.send("6/3\r")
        s.expect(r"Out\[4\]= 2")
        check(b">/help" not in s.buf[start:], "the menu opened inside an expression")
        quit_session(s)
    finally:
        s.close()


@case
def slash_quit_and_readline(tmp):
    """/quit ends the session with status 0; the readline session runs
    commands too (no menu)."""
    s = interactive(tmp, extra=NOHIST)
    try:
        s.send("/quit\r")
        check(s.exit_code() == 0, "/quit did not exit 0")
    finally:
        s.close()
    s = interactive(tmp, extra=dict(NOHIST, MATHILDA_EDITOR="readline"))
    try:
        s.send("/timing off\r")
        s.expect(r"Timing notes off\.")
        wait_prompt(s, 1)
        s.send("/theme\r")                        # no picker: the list
        s.expect(r"Themes\r\n")
        s.expect(r"/theme NAME chooses one")
        wait_prompt(s, 1)
        s.send("2+2\r")
        s.expect(r"Out\[1\]= 4")
        quit_session(s)
    finally:
        s.close()


def theme_session(tmp, depth="24", **kw):
    extra = dict(NOHIST, MATHILDA_COLOR_DEPTH=depth, **(kw.pop("extra", None) or {}))
    return interactive(tmp, extra=extra, drop=("NO_COLOR",), **kw)


DARK_OUT = b"\x1b[1;38;2;198;120;221mOut["
LIGHT_OUT = b"\x1b[1;38;2;130;80;223mOut["


@case
def theme_picker_apply(tmp):
    """The picker previews live, Return applies and saves the theme, and the
    next session starts with it."""
    s = theme_session(tmp)
    try:
        s.send("1+1\r")
        s.expect(r"Out\[1\]=")
        check(DARK_OUT in s.buf, "the default theme is not dark at 24-bit")
        s.send("/theme\r")
        s.expect(r"Theme\x1b\[0m")
        s.expect(r"(?s)Power::infy.*took 3\.21 s")   # the preview
        s.drain(0.2)
        start = len(s.buf)
        s.send("\x1b[B")                            # Down: light, previewed
        s.drain(0.3)
        check(b"38;2;130;80;223m" in s.buf[start:], "no live preview of light")
        s.send("\r")
        s.expect(r"Theme \x1b\[1mlight\x1b\[0m.*saved in ~/\.mathilda/config")
        check(b"\x1b[?25h" in s.buf[start:], "the cursor was not shown again")
        wait_prompt(s, 2)
        s.send("2+2\r")
        s.expect(r"Out\[2\]=")
        check(LIGHT_OUT in s.buf, "the result is not in the light theme")
        quit_session(s)
    finally:
        s.close()
    with open(os.path.join(tmp, ".mathilda", "config")) as f:
        check(f.read().strip() == "theme = light", "config not written")
    s = theme_session(tmp)
    try:
        s.send("1+1\r")
        s.expect(r"Out\[1\]=")
        check(LIGHT_OUT in s.buf, "the saved theme was not loaded")
        quit_session(s)
    finally:
        s.close()
    s = theme_session(tmp, extra={"MATHILDA_THEME": "mono"})
    try:
        s.send("1+1\r")
        s.expect(r"Out\[1\]=")
        check(b"\x1b[1mOut[1]=" in s.buf, "MATHILDA_THEME did not override the file")
        quit_session(s)
    finally:
        s.close()


@case
def theme_picker_cancel(tmp):
    """Esc reverts the preview and writes nothing."""
    s = theme_session(tmp)
    try:
        s.send("/theme\r")
        s.expect(r"took 3\.21 s")
        s.send("\x1b[B\x1b[B")
        s.drain(0.3)
        s.send("\x1b")
        s.expect(r"\x1b\[\?25h")                     # the picker is gone
        s.drain(0.2)
        s.send("1+1\r")
        m = s.expect(rb"Out\[1\]=\x1b\[0m \x1b\[[0-9;]*m2")   # the result, not the preview
        i = m.start()
        check(DARK_OUT in s.buf[i - 40:i + 8], "Esc did not revert to dark: %r" % s.buf[i - 40:i + 8])
        check(not os.path.exists(os.path.join(tmp, ".mathilda", "config")), "Esc saved the theme")
        quit_session(s)
    finally:
        s.close()


@case
def theme_depths(tmp):
    """256 colours, 16 colours, and NO_COLOR (no escape sequence at all)."""
    s = theme_session(tmp, depth="256")
    try:
        s.send("1+1\r")
        s.expect(r"Out\[1\]=")
        check(re.search(rb"\x1b\[1;38;5;\d+mOut\[1\]=", s.buf), "no 256-colour Out")
        quit_session(s)
    finally:
        s.close()
    s = theme_session(tmp, depth="16")
    try:
        s.send("1+1\r")
        s.expect(r"Out\[1\]=")
        check(b"\x1b[1;35mOut[1]=" in s.buf, "16-colour Out is not the original colour")
        quit_session(s)
    finally:
        s.close()
    s = interactive(tmp, extra=dict(NOHIST, MATHILDA_COLOR_DEPTH="24", MATHILDA_THEME="light"))
    try:                                          # NO_COLOR is set by default here
        s.send("/theme high-contrast\r")
        s.expect(r"Colour is off in this session")
        s.send("1+1\r")
        s.expect(r"Out\[1\]= 2")
        check(no_sgr(s.buf), "SGR with NO_COLOR")
        quit_session(s)
    finally:
        s.close()


MOD_INIT = r'''
greet[s_] := "Hello, " <> If[s == "", "world", s] <> "!";
SlashCommand["greet", greet, "say hello: /greet [name]"];
$Post = (# &);
'''
MOD_BAD = "SlashCommand[\"help\", Print];\n1/0;\n"
MOD_QUIT = "SlashCommand[\"boom\", (Message[General::ovfl]; Abort[]) &, \"always fails\"];\nQuit[];\n"


def write_mods(tmp):
    d = os.path.join(tmp, ".mathilda")
    os.makedirs(os.path.join(d, "mods"), exist_ok=True)
    with open(os.path.join(d, "init.m"), "w") as f:
        f.write(MOD_INIT)
    with open(os.path.join(d, "mods", "bad.m"), "w") as f:
        f.write(MOD_BAD)
    with open(os.path.join(d, "mods", "quitter.m"), "w") as f:
        f.write(MOD_QUIT)


@case
def mods_load_and_list(tmp):
    """init.m and mods/*.m load after the banner; failures are reported and
    the session carries on; /mods lists files, commands and hooks."""
    write_mods(tmp)
    s = new_session(tmp, extra=NOHIST)
    try:
        s.expect(r"Type / for commands")
        s.expect(r"SlashCommand::builtin: /help is a built-in command")
        s.expect(r"mods/bad\.m reported messages; the session continues")
        s.expect(r"mods/quitter\.m called Quit\[\], ignored")
        wait_prompt(s, 1)
        s.send("/greet Ann\r")
        s.expect(r"Hello, Ann!\r\n")             # text, unquoted
        wait_prompt(s, 1)
        s.send("/mods\r")
        s.expect(r"init\.m +loaded \xc2\xb7 /greet \xc2\xb7 sets \$Post")
        s.expect(r"mods/bad\.m +2 messages")
        s.expect(r"mods/quitter\.m +Quit\[\] ignored \xc2\xb7 /boom")
        s.send("/boom\r")
        s.expect(r"General::ovfl")
        s.expect(r"/boom: aborted")
        wait_prompt(s, 1)
        s.send("/\x1b[A")                         # the menu lists mod commands too
        s.expect(r"/greet +say hello")
        s.send("\x03")
        s.send("SlashCommand[\"greet\", None]; 1+1\r")
        s.expect(r"Out\[1\]= 2")
        s.send("/greet\r")
        s.expect(r"Unknown command /greet")
        quit_session(s)
    finally:
        s.close()
    s = interactive(tmp, extra=dict(NOHIST, MATHILDA_MODS="0"))
    try:
        check(b"reported messages" not in s.buf, "MATHILDA_MODS=0 still loaded mods")
        s.send("/mods\r")
        s.expect(r"MATHILDA_MODS=0: no mods were loaded")
        quit_session(s)
    finally:
        s.close()
    s = interactive(tmp, extra=dict(NOHIST, MATHILDA_CONFIG_DIR=""))
    try:
        s.send("/mods\r")
        s.expect(r"Mods are off: MATHILDA_CONFIG_DIR is empty")
        quit_session(s)
    finally:
        s.close()


@case
def monitor_status_line(tmp):
    """Monitor[expr, mon] shows mon's value in the status line."""
    s = spinner_session(tmp, colour=False)
    try:
        s.send("Monitor[Do[Pause[0.05], {i, 30}], {i, \"of 30\"}]\r")
        s.expect(r"Out\[1\]= Null", 10.0)
        vals = [int(v) for v in re.findall(rb"Do \d+/30 [^\r]*?\{(\d+), \"of 30\"\}", s.buf)]
        check(len(vals) >= 5, "no monitor text in the frames: %r" % s.buf[-600:])
        check(vals == sorted(vals) and vals[-1] > vals[0], "monitor not advancing: %r" % vals)
        s.send("Monitor[Total[Range[10]], Print[\"x\"]]\r")   # quick: mon never evaluated
        s.expect(r"Out\[2\]= 55")
        check(b"\r\nx\r\n" not in s.buf, "mon was evaluated for a quick input")
        quit_session(s)
    finally:
        s.close()


@case
def timing_and_statusline_commands(tmp):
    s = interactive(tmp, extra=NOHIST)
    try:
        s.send("/timing off\r")
        s.expect(r"Timing notes off\.")
        s.send("Pause[1.1]\r")
        s.expect(r"Out\[1\]= Null")
        wait_prompt(s, 2)
        check(b"took " not in s.buf, "took note with /timing off")
        s.send("/timing always\r")
        s.expect(r"every result")
        s.send("1+1\r")
        s.expect(r"Out\[2\]= 2\r\ntook \d\.\d\d s")
        s.send("/statusline\r")
        s.expect(r"Status line on")
        s.send("2+2\r")
        s.expect(r"last \d+\.\d\d s")
        s.send("/statusline off\r")
        s.expect(r"Status line off")
        s.send("/clear\r")
        s.expect(r"\x1b\[H\x1b\[2J")
        wait_prompt(s, 4)
        quit_session(s)
    finally:
        s.close()


# Ctrl-Z needs job control: pty.fork() makes the child a session leader whose
# process group is orphaned, and the kernel discards stop signals there. A
# tiny "shell" runs Mathilda in its own foreground group, reports a stop, and
# brings it back (tcsetpgrp, then SIGCONT) on SIGUSR1.
JOB_SHELL = r'''
import os, signal, sys
signal.signal(signal.SIGTTOU, signal.SIG_IGN)
pid = os.fork()
if pid == 0:
    os.setpgid(0, 0)
    os.tcsetpgrp(0, os.getpgrp())
    os.execv(sys.argv[1], sys.argv[1:])
try:
    os.setpgid(pid, pid)
except OSError:
    pass
try:
    os.tcsetpgrp(0, pid)
except OSError:
    pass
def fg(sig, frame):
    os.tcsetpgrp(0, pid)
    os.kill(pid, signal.SIGCONT)
signal.signal(signal.SIGUSR1, fg)
os.write(1, b"[job %d]\r\n" % pid)
while True:
    p, st = os.waitpid(pid, os.WUNTRACED)
    if os.WIFSTOPPED(st):
        os.tcsetpgrp(0, os.getpgrp())
        os.write(1, b"[stopped]\r\n")
        continue
    os._exit(os.WEXITSTATUS(st) if os.WIFEXITED(st) else 128 + os.WTERMSIG(st))
'''


def job_session(tmp, **kw):
    global WRAPPER
    path = os.path.join(tmp, "jobshell.py")
    with open(path, "w") as f:
        f.write(JOB_SHELL)
    saved = WRAPPER
    WRAPPER = [sys.executable, path] + saved
    try:
        s = new_session(tmp, **kw)
    finally:
        WRAPPER = saved
    s.expect(r"\[job (\d+)\]")
    wait_prompt(s, 1)
    return s


@case
def ctrl_z_at_prompt(tmp):
    """Ctrl-Z at the prompt: cooked terminal and paste mode off while stopped;
    on fg raw mode again and the buffer redrawn."""
    s = job_session(tmp, extra=NOHIST)
    try:
        s.send("abc")
        s.drain(0.2)
        s.send("\x1a")
        s.expect(r"\[stopped\]")
        lf = tty_flags(s.fd)
        check((lf & termios.ICANON) and (lf & termios.ECHO), "terminal not cooked while stopped")
        check(b"\x1b[?2004l" in s.buf[-80:], "bracketed paste left on while stopped")
        start = len(s.buf)
        os.kill(s.pid, signal.SIGUSR1)
        s.expect(r"In\[1\]:= (\x1b\[K)?abc")
        s.drain(0.2)
        check(not (tty_flags(s.fd) & termios.ICANON), "not raw after fg")
        check(b"\x1b[?2004h" in s.buf[start:], "paste mode not back after fg")
        s.send("+1\r")
        s.expect(r"Out\[1\]= 1 \+ abc")
        quit_session(s)
    finally:
        s.close()


@case
def ctrl_z_under_spinner(tmp):
    """Ctrl-Z during an evaluation: the line erased and ECHO back while
    stopped; on fg ECHO off again, the spinner resumes, the result arrives."""
    s = job_session(tmp, extra={"MATHILDA_HISTFILE": "", "NO_COLOR": "1"})
    try:
        s.send("Pause[2]\r")
        s.drain(0.7)
        check(BRAILLE.search(s.buf), "spinner not running")
        check(not (tty_flags(s.fd) & termios.ECHO), "ECHO on under the spinner")
        s.send("\x1a")
        s.expect(r"\[stopped\]")
        check((tty_flags(s.fd) & termios.ECHO), "ECHO not restored while stopped")
        before = s.buf[:s.buf.rfind(b"[stopped]")]
        last_frame = before.rfind(b"\xe2\xa0")
        check(re.search(rb"\r +\r", before[last_frame:]), "spinner line not erased before the stop")
        quiet = len(s.buf)
        s.drain(0.5)
        check(not BRAILLE.search(s.buf[quiet:]), "drew while stopped")
        os.kill(s.pid, signal.SIGUSR1)
        s.drain(0.6)
        check(BRAILLE.search(s.buf[quiet:]), "spinner did not resume")
        check(not (tty_flags(s.fd) & termios.ECHO), "ECHO not cleared again after fg")
        s.expect(r"Out\[1\]= Null", 10.0)
        wait_prompt(s, 2)
        s.send("1+1\r")
        s.expect(r"Out\[2\]= 2")
        quit_session(s)
    finally:
        s.close()


@case
def cursor_query_no_wait(tmp):
    """A terminal that never answers ESC [ 6 n costs nothing: the query is
    not waited for, and after two unanswered ones it is not sent again."""
    s = interactive(tmp, extra=NOHIST)
    try:
        gaps = []
        for n in range(1, 6):
            s.send("%d\r" % n)
            s.expect(r"Out\[%d\]= %d" % (n, n))
            t0 = time.monotonic()
            wait_prompt(s, n + 1)
            gaps.append(time.monotonic() - t0)
        check(s.buf.count(b"\x1b[6n") <= 2, "queried %d times" % s.buf.count(b"\x1b[6n"))
        gaps.sort()
        check(gaps[2] < 0.06, "prompts took %r s (a blocking query waits 80 ms)" % gaps)
        quit_session(s)
    finally:
        s.close()


@case
def cursor_query_answered(tmp):
    """A terminal that answers: the reply is consumed, never typed."""
    s = new_session(tmp, extra=NOHIST)
    try:
        s.expect(r"\x1b\[6n")
        s.send("\x1b[24;1R")                     # the reply, after the prompt
        wait_prompt(s, 1)
        s.send("7\r")
        s.expect(r"Out\[1\]= 7")
        s.expect(r"\x1b\[6n")                    # answered once: asked again
        s.send("\x1b[24;1R")
        wait_prompt(s, 2)
        s.send("8\r")
        s.expect(r"Out\[2\]= 8")
        check(b"24;1R" not in strip_esc(s.buf).replace(b"\x1b", b""), "the reply was typed")
        quit_session(s)
    finally:
        s.close()


@case
def redo_key(tmp):
    s = interactive(tmp, extra=NOHIST)
    try:
        s.send("1+2")
        s.send("\x1f")                            # undo
        s.send("\x1b_")                           # Alt-_: redo
        s.send("\r")
        s.expect(r"Out\[1\]= 3")
        quit_session(s)
    finally:
        s.close()


def harness_can_fail(tmp):
    """Deliberately wrong expectation: must FAIL (MATHILDA_PTY_SELFTEST=1)."""
    s = new_session(tmp)
    try:
        wait_prompt(s, 1)
        s.send("1+1\r")
        s.expect(r"Out\[1\]= 3", timeout=3)
    finally:
        s.close()


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def record(directory):
    os.makedirs(directory, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="mth-pty-") as tmp:
        res = golden_runs(tmp)
    for name, data in sorted(res.items()):
        with open(os.path.join(directory, name), "wb") as f:
            f.write(data)
        print("recorded %s (%d bytes)" % (os.path.join(directory, name), len(data)))


def main():
    global BIN, WRAPPER, TIME_SCALE, GOLDEN_DIR
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bin", default="./Mathilda")
    ap.add_argument("--wrapper", default="")
    ap.add_argument("--record", metavar="DIR")
    ap.add_argument("--golden", default="tests/pty/golden")
    ap.add_argument("-k", dest="pattern")
    a = ap.parse_args()

    BIN = a.bin
    GOLDEN_DIR = a.golden
    if not (os.path.isfile(BIN) and os.access(BIN, os.X_OK)):
        print("repl_pty_test: %s not built; skipping" % BIN)
        return 77
    if not BIN.startswith("/"):
        BIN = os.path.abspath(BIN)
    if a.wrapper:
        WRAPPER = shlex.split(a.wrapper)
        TIME_SCALE = 10.0
    if a.record:
        record(a.record)
        return 0

    cases = list(CASES)
    if os.environ.get("MATHILDA_PTY_SELFTEST") == "1":
        cases.append(harness_can_fail)
    if a.pattern:
        rx = re.compile(a.pattern)
        cases = [c for c in cases if rx.search(c.__name__)]
    if not cases:
        print("repl_pty_test: no case matches %r" % a.pattern)
        return 1

    failed = []
    for fn in cases:
        t0 = time.monotonic()
        with tempfile.TemporaryDirectory(prefix="mth-pty-") as tmp:
            try:
                fn(tmp)
                print("PASS %-28s (%.1fs)" % (fn.__name__, time.monotonic() - t0))
            except Exception as e:  # noqa: BLE001 -- report every failure kind
                failed.append(fn.__name__)
                print("FAIL %-28s (%.1fs)\n%s" % (fn.__name__, time.monotonic() - t0, e))
        sys.stdout.flush()
    print("%d/%d cases passed" % (len(cases) - len(failed), len(cases)))
    if failed:
        print("failed: " + " ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
