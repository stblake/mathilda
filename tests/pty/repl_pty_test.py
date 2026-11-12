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


# Readline's own bracketed-paste toggles, already emitted before this work.
BRACKETED_PASTE = re.compile(rb"\x1b\[\?2004[hl]")


def no_sgr(data):
    return b"\x1b[" not in BRACKETED_PASTE.sub(b"", data)


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
            m = re.search(rb"(\r\n|\r\x1b\[2K|\r)line%d\r\n" % i, seg)
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
        s.expect(r"Out\[1\]=\S* \$Aborted", 2.0)
        check(time.monotonic() - t0 < 1.0, "abort took too long")
        # the transcript is PR 1's: ^C, newline, then the result -- the spinner
        # line erased first and no stale frame left behind
        check(re.search(rb"\r\x1b\[2K\^C\r\n\x1b\[1;35mOut\[1\]=", s.buf),
              "unexpected Ctrl-C transcript: %r" % s.buf[-300:])
        s.send("1+1\r")
        s.expect(r"Out\[2\]=\S* 2")
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
