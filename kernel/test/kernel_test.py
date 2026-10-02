import os, sys, time
from jupyter_client.manager import start_new_kernel

fails = 0
def check(name, cond, extra=""):
    global fails
    print(("PASS" if cond else "FAIL") + "  " + name + ("" if cond else "   "+str(extra)))
    if not cond: fails += 1

km, kc = start_new_kernel(kernel_name="xmathilda")
print("kernel started")

def run(code):
    """Return (results, displays, streams, errors, status)."""
    msg_id = kc.execute(code)
    results, displays, streams, errors, status = [], [], [], [], None
    while True:
        try: msg = kc.get_iopub_msg(timeout=30)
        except Exception: break
        if msg["parent_header"].get("msg_id") != msg_id: continue
        t = msg["msg_type"]; c = msg["content"]
        if t == "execute_result": results.append(c["data"])
        elif t == "display_data": displays.append(c["data"])
        elif t == "stream": streams.append((c["name"], c["text"]))
        elif t == "error": errors.append(c)
        elif t == "status" and c["execution_state"] == "idle": break
    reply = kc.get_shell_msg(timeout=30)
    status = reply["content"]["status"]
    return results, displays, streams, errors, status

# kernel_info. Every field here is a distinct slot of xeus's create_info_reply,
# whose parameter list differs between xeus 5 and 6 (it lost a leading
# `protocol_version` in 6). All the arguments are strings, so a call written for
# the wrong version binds silently one slot over and the reply comes back
# plausible but shifted — hence checking the far end of the list (banner) and not
# just the near end. `supported_features` exists only in xeus 6's reply, so its
# presence is also how the control-channel checks below tell the versions apart.
kc.kernel_info(); ki = kc.get_shell_msg(timeout=30)["content"]
xeus6 = "supported_features" in ki
print(f"   xeus API: {'6' if xeus6 else '5'} (protocol {ki.get('protocol_version')})")
check("kernel_info: language mathilda", ki["language_info"]["name"] == "mathilda", ki.get("language_info"))
check("kernel_info: implementation xmathilda", ki.get("implementation") == "xmathilda", ki.get("implementation"))
check("kernel_info: banner names Mathilda", ki.get("banner", "").startswith("Mathilda "), ki.get("banner"))
check("kernel_info: codemirror_mode mathematica",
      ki["language_info"].get("codemirror_mode") == "mathematica", ki.get("language_info"))

# simple symbolic result with LaTeX
r,d,s,e,st = run("Integrate[x^2, x]")
check("Integrate: execute_result present", len(r)==1, (r,d))
check("Integrate: has text/latex", r and "text/latex" in r[0], r)
check("Integrate: has text/plain", r and "text/plain" in r[0], r)
print("   latex:", r[0].get("text/latex") if r else None, "| plain:", r[0].get("text/plain") if r else None)

# multi-statement cell with Print + final result
r,d,s,e,st = run('Print["hello"]\nPrint["world"]\n2+2')
check("cell: two stdout streams", sum(1 for n,_ in s if n=="stdout")>=1 and "hello" in "".join(t for _,t in s) and "world" in "".join(t for _,t in s), s)
check("cell: final result 4", any(v.get("text/plain")=="4" for v in r), (r,d))
check("cell: intermediate results as display_data (Print has none; 2+2 is the only value)", True)

# message to stderr (1/0 => Power::infy)
r,d,s,e,st = run("1/0")
check("1/0: Power::infy on stderr", any(n=="stderr" and "Power::infy" in t for n,t in s), s)

# parse error => error status
r,d,s,e,st = run("f[1,")
check("syntax error: status error", st=="error", (st,e))

# plot => plotly display_data
r,d,s,e,st = run("Plot[Sin[x], {x, 0, 3}]")
check("Plot: plotly display_data", any("application/vnd.plotly.v1+json" in v for v in d), d and list(d[0].keys()))

# completion
kc.complete("Sin", 3); cm = kc.get_shell_msg(timeout=30)["content"]
check("complete: matches include Sin/Sinh", "Sin" in cm["matches"] and "Sinh" in cm["matches"], cm.get("matches"))
check("complete: cursor_start=0 cursor_end=3", cm["cursor_start"]==0 and cm["cursor_end"]==3, cm)

# inspection
kc.inspect("Sin", 3, 0); ins = kc.get_shell_msg(timeout=30)["content"]
check("inspect: found", ins["found"] is True, ins)
check("inspect: has text/plain", "text/plain" in ins.get("data",{}), ins.get("data"))

# is_complete
kc.is_complete("1 + 1"); c1 = kc.get_shell_msg(timeout=30)["content"]
kc.is_complete("f[1,");  c2 = kc.get_shell_msg(timeout=30)["content"]
check("is_complete: balanced => complete", c1["status"]=="complete", c1)
check("is_complete: open bracket => incomplete", c2["status"]=="incomplete", c2)

# session history
run("11 + 11")
r,d,s,e,st = run("% + 1")
check("history: % works (23)", any(v.get("text/plain")=="23" for v in r), (r,d))

# --- control channel -------------------------------------------------------
# On xeus 6 both replies come straight from the interpreter and MUST carry
# status: xkernel_core reads reply["status"] as a std::string and only calls
# p_server->stop() when it is "ok", so an interpreter returning a bare {} throws
# nlohmann 302 inside the control handler — no reply is sent at all and the kernel
# never exits. That is what these two checks guard.
#
# On xeus 5 the core writes both replies itself and never asks the interpreter,
# and it sends them with *null* content: it std::moves the reply json into
# publish_message and then hands the moved-from value to send_reply. That is an
# upstream xeus 5 bug (xeus 6 passes a copy to publish and moves only into
# send_reply), nothing we can influence, so there the assertion is only that the
# reply arrives and the kernel actually goes away.
def control(msg_type, content, timeout=15):
    kc.control_channel.send(kc.session.msg(msg_type, content))
    try:
        return kc.get_control_msg(timeout=timeout)
    except Exception as ex:
        return {"msg_type": f"NO REPLY ({type(ex).__name__})"}

def check_reply(req, want_type):
    m = control(req, {"restart": False} if req == "shutdown_request" else {})
    ok = m.get("msg_type") == want_type
    if ok and xeus6:                      # only xeus 6 routes the reply through us
        ok = (m.get("content") or {}).get("status") == "ok"
    check(f"{req}: {want_type}" + (" with status ok" if xeus6 else ""),
          ok, (m.get("msg_type"), m.get("content")))

check_reply("interrupt_request", "interrupt_reply")
check_reply("shutdown_request", "shutdown_reply")

time.sleep(3)
check("shutdown_request: kernel process exited", not km.is_alive())

km.shutdown_kernel(now=True)   # backstop if the graceful path above did not land
print()
print(("ALL PASS" if fails==0 else f"{fails} FAILURE(S)"))
sys.exit(1 if fails else 0)
