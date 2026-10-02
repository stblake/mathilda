# Issue #85 — the xeus kernel only builds against xeus 6

Reported by mkoeppe packaging Mathilda 0.252 for Sage: `kernel/` fails to
compile because `mathilda_interpreter.{hpp,cpp}` targets the xeus **6**
`xinterpreter` ABI while his toolchain ships xeus **5.2.x**.
`kernel/README.md` already claims "xeus ≥ 5" — the code never honoured it.

## The three xeus 5 → 6 deltas that touch us (5.2.8 vs 6.0.6 headers diffed)

| | xeus 5.0 – 5.2.8 | xeus 6.0+ |
|---|---|---|
| shutdown | `void shutdown_request_impl()` | `nl::json shutdown_request_impl(bool restart)` |
| interrupt | *absent* | `nl::json interrupt_request_impl()` |
| info reply | `create_info_reply(protocol_version, implementation, …, banner, debugger, help_links)` | `create_info_reply(implementation, …, banner, help_links, supported_features)` |

The first two are the reported compile errors. The third is worse: all ten of our
arguments are `std::string`, so against xeus 5 the call binds with every argument
**shifted one slot left** — `"xmathilda"` becomes the protocol version and the
banner is dropped. It compiles and lies.

Everything else is source-identical across 5.x/6.x, so `main.cpp` needs no
change. xeus 4 is out of scope (`execute_request_impl` took an `xrequest_context`).

## Second, independent bug (found while reading xeus 6's core)

xeus 6 `xkernel_core::shutdown_request` does `std::string reply_status =
reply["status"];` and only calls `p_server->stop()` when it is `"ok"`. Our two
handlers returned a bare `nl::json::object()`, so nlohmann threw 302 and no reply
was ever sent. Confirmed against the installed pre-fix kernel:

```
ERROR: received bad message: [json.exception.type_error.302] type must be string, but is null
Message type: interrupt_request        -> NO REPLY
Message type: shutdown_request         -> NO REPLY
process still alive after shutdown_request: True
```

Unnoticed because `kernel/test/kernel_test.py` tore the kernel down with
`shutdown_kernel(now=True)` — a SIGKILL that never sends the request.

## Plan

- [x] Diff every xeus header the kernel touches, 5.2.8 vs 6.0.6
- [x] Reproduce the control-channel bug against the installed kernel
- [x] `mathilda_interpreter.hpp` — `MATHILDA_XEUS_6` macro + guarded declarations
- [x] `mathilda_interpreter.cpp` — guarded definitions; `create_shutdown_reply` /
      `create_interrupt_reply` on xeus 6 (fixes bug 2); guarded `create_info_reply`
- [x] `CMakeLists.txt` — `find_package(xeus 5.0 REQUIRED)`, report the API branch
- [x] `kernel_test.py` — graceful interrupt + shutdown checks (regression test
      for bug 2); must FAIL on the pre-fix binary
- [x] `README.md` — supported xeus versions
- [x] Verify: build + test under xeus 6.0.6
- [x] Verify: build + test under xeus 5.2.8 (fresh env) — this is what closes #85
- [x] Version bump 0.253 → 0.254, changelog, commit, tag
- [x] Reply on the issue (v0.254, `824721c2`, tagged)

## Review

**What changed.** `kernel/src/mathilda_interpreter.hpp` now derives one macro,
`MATHILDA_XEUS_6`, from `XEUS_VERSION_MAJOR` and declares the two divergent
virtuals under it; `.cpp` carries the matching pair of definitions plus a guarded
`create_info_reply` (the xeus 5 branch passes `XEUS_KERNEL_PROTOCOL_VERSION`
first and `/*debugger*/ false` in its own slot). On xeus 6 the handlers now
return `xeus::create_shutdown_reply(restart)` / `create_interrupt_reply()`.
`CMakeLists.txt` floors `xeus` at 5.0 and `xeus-zmq` at 3.0 and prints which API
branch it selected. `kernel_test.py` gained the two control-channel checks and
now shuts the kernel down the way Jupyter does.

**Verification.** 19/19 checks pass under xeus 6.0.6 AND under a fresh xeus
5.2.8 env (`~/micromamba/envs/mathilda-xeus5`, xeus-zmq 3.1.0) — the first time
the kernel has ever been built against xeus 5. The two new checks fail on the
pre-fix binary (`NO REPLY` / still alive), so they earn their place. Root `make`
and `make check-c99` clean.

**Loose end worth knowing.** A true mid-computation interrupt still needs an
abort flag in the evaluator loop; `interrupt_request` now *acknowledges*
correctly rather than throwing, but it does not stop a running `Integrate[]`.
That was already on the kernel's follow-up list and is unchanged by this fix.
