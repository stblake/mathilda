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
    extra = kw.pop("extra", None)
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
        out = p.communicate(stdin_data, timeout=timeout * TIME_SCALE)[0]
        err = b""
        while True:
            r, _, _ = select.select([master], [], [], 0.2)
            if not r:
                break
            try:
                chunk = os.read(master, 65536)
            except OSError:
                break
            if not chunk:
                break
            err += chunk
        os.close(master)
        return p.returncode, out, err
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
