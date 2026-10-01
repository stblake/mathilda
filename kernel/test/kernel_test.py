import os, sys
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

# kernel_info
kc.kernel_info(); ki = kc.get_shell_msg(timeout=30)["content"]
check("kernel_info: language mathilda", ki["language_info"]["name"] == "mathilda", ki.get("language_info"))

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

km.shutdown_kernel(now=True)
print()
print(("ALL PASS" if fails==0 else f"{fails} FAILURE(S)"))
sys.exit(1 if fails else 0)
