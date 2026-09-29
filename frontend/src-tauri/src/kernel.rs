// kernel.rs — manages the Mathilda CAS child process

use serde_json::{json, Value};
use std::sync::Arc;
use std::time::Duration;
use tauri::async_runtime::Receiver;
use tauri::ipc::Channel;
use tauri::Emitter;
use tauri_plugin_shell::process::{CommandChild, CommandEvent};
use tauri_plugin_shell::ShellExt;
use tokio::sync::Mutex;

// ---------------------------------------------------------------------------

struct KernelState {
    child: CommandChild,
    rx: Receiver<CommandEvent>,
    next_id: u32,
}

/// Handle stored in Tauri State — wraps the Mathilda child process.
pub struct MathildaKernel {
    state: Arc<Mutex<Option<KernelState>>>,
    app: tauri::AppHandle,
}

impl MathildaKernel {
    /// Create a kernel handle synchronously — state starts empty.
    /// Call `start()` asynchronously after managing.
    pub fn empty(app: tauri::AppHandle) -> Self {
        Self {
            state: Arc::new(Mutex::new(None)),
            app,
        }
    }

    /// Spawn Mathilda sidecar and wait for ping/pong handshake.
    pub async fn new(app: tauri::AppHandle) -> Result<Self, String> {
        let kernel = Self::empty(app);
        kernel.spawn_inner().await?;
        Ok(kernel)
    }

    /// Start the kernel asynchronously (call after `empty()` + manage).
    pub async fn start(&self) -> Result<(), String> {
        self.spawn_inner().await
    }

    async fn spawn_inner(&self) -> Result<(), String> {
        let (rx, child) = self
            .app
            .shell()
            .sidecar("mathilda")
            .map_err(|e| format!("sidecar lookup: {e}"))?
            .spawn()
            .map_err(|e| format!("spawn: {e}"))?;

        {
            let mut guard = self.state.lock().await;
            *guard = Some(KernelState {
                child,
                rx,
                next_id: 1,
            });
        }

        self.ping_inner().await
    }

    async fn ping_inner(&self) -> Result<(), String> {
        let mut guard = self.state.lock().await;
        let state = guard.as_mut().ok_or("no kernel")?;

        state
            .child
            .write(b"{\"type\":\"ping\"}\n")
            .map_err(|e| format!("write: {e}"))?;

        let deadline = tokio::time::Instant::now() + Duration::from_secs(10);
        loop {
            let remaining = deadline.saturating_duration_since(tokio::time::Instant::now());
            if remaining.is_zero() {
                return Err("Kernel ping timed out".into());
            }
            match tokio::time::timeout(remaining, state.rx.recv()).await {
                Ok(Some(CommandEvent::Stdout(bytes))) => {
                    if String::from_utf8_lossy(&bytes).contains("\"pong\"") {
                        return Ok(());
                    }
                }
                Ok(Some(CommandEvent::Terminated(_))) => {
                    return Err("Kernel terminated before pong".into());
                }
                Ok(None) | Err(_) => {
                    return Err("Kernel ping timed out".into());
                }
                _ => {}
            }
        }
    }

    /// Ping the kernel publicly (for frontend health check).
    pub async fn ping(&self) -> Result<(), String> {
        self.ping_inner().await
    }

    /// Evaluate `expr`, forwarding output messages to `channel` until "done".
    pub async fn evaluate(&self, expr: String, channel: Channel<Value>) -> Result<(), String> {
        let mut guard = self.state.lock().await;
        let state = guard.as_mut().ok_or("kernel not running")?;

        let id = state.next_id;
        state.next_id = state.next_id.wrapping_add(1);

        let request = cell_request(id, &expr)?;
        state
            .child
            .write(request.as_bytes())
            .map_err(|e| format!("write: {e}"))?;

        loop {
            match state.rx.recv().await {
                Some(CommandEvent::Stdout(bytes)) => {
                    match classify_stdout_line(&String::from_utf8_lossy(&bytes), id) {
                        LineAction::Forward(msg) => {
                            channel.send(msg).map_err(|e| format!("channel: {e}"))?;
                        }
                        LineAction::Done { memory } => {
                            /* `done` carries the kernel's resident memory, which the status bar
                            shows. Forwarded as its own event rather than through `channel`:
                            the channel is the CELL's output stream and a memory reading is
                            not output -- putting it there would make every consumer of cell
                            output filter it out. It is attached to `done` on the kernel side
                            because memory can only change when something was evaluated, so a
                            poll would either lag or add traffic for a number already in
                            hand. */
                            if let Some(bytes) = memory {
                                let _ = self.app.emit("kernel-memory", bytes);
                            }
                            break;
                        }
                        LineAction::Ignore => {}
                    }
                }
                Some(CommandEvent::Terminated(p)) => {
                    return Err(format!("Kernel died (code {:?})", p.code));
                }
                Some(CommandEvent::Stderr(bytes)) => {
                    /* The kernel captures its own diagnostics during an evaluation and sends
                    them as "message" lines, so stderr here is whatever escaped that capture.
                    It arrives while this cell is evaluating, so it is shown in the cell as a
                    warning rather than only logged where nobody looks. */
                    let text = String::from_utf8_lossy(&bytes);
                    log::debug!("[kernel] {text}");
                    if let Some(msg) = stderr_message(&text, id) {
                        channel.send(msg).map_err(|e| format!("channel: {e}"))?;
                    }
                }
                _ => {}
            }
        }
        Ok(())
    }

    /// Kill the current kernel (state → None). Frontend should call
    /// restart_kernel afterwards if it wants a fresh session.
    pub async fn kill(&self) -> Result<(), String> {
        let mut guard = self.state.lock().await;
        if let Some(old) = guard.take() {
            // Graceful quit first (best-effort, kernel may be busy).
            let mut child = old.child;
            let _ = child.write(b"{\"type\":\"quit\"}\n");
            drop(old.rx);
            tokio::time::sleep(Duration::from_millis(300)).await;
            let _ = child.kill(); // consumes child
        }
        Ok(())
    }

    /// Restart: kill current kernel and spawn a fresh one with backoff.
    pub async fn restart(&self) -> Result<(), String> {
        self.kill().await?;

        let mut delay = Duration::from_millis(500);
        for attempt in 1..=3u32 {
            log::info!("Restart attempt {attempt}/3");
            match self.spawn_inner().await {
                Ok(()) => return Ok(()),
                Err(e) => {
                    log::warn!("Attempt {attempt} failed: {e}");
                    if attempt < 3 {
                        tokio::time::sleep(delay).await;
                        delay *= 2;
                    } else {
                        return Err(format!("Kernel failed to start after 3 attempts: {e}"));
                    }
                }
            }
        }
        Ok(())
    }

    /// Interrupt: kill the child so the current evaluation stops.
    /// The state is set to None; caller should restart_kernel.
    pub async fn interrupt(&self) -> Result<(), String> {
        self.kill().await
    }

    pub async fn is_running(&self) -> bool {
        self.state.lock().await.is_some()
    }
}

// ---------------------------------------------------------------------------
// Protocol: what one line from the kernel means for the evaluation `id`.
//
// Pure functions, so the routing is unit tested without spawning a process.

/// The request line for evaluating `expr` as notebook cell `id`.
///
/// `"cell": true` asks the kernel for notebook semantics (see src/repl.c): several
/// statements per cell, Print output and messages captured and sent as "stream" /
/// "message" lines before each result, and no result for `;` or Null. A kernel
/// that predates the flag ignores it and answers as before, with non-JSON Print
/// lines and stderr then caught by the fallbacks below.
pub(crate) fn cell_request(id: u32, expr: &str) -> Result<String, String> {
    let line = serde_json::to_string(&json!({ "id": id, "expr": expr, "cell": true }))
        .map_err(|e| e.to_string())?;
    Ok(line + "\n")
}

/// What to do with one stdout line while evaluation `id` is in flight.
#[derive(Debug, PartialEq)]
pub(crate) enum LineAction {
    /// A message for the evaluating cell (result, error, Print text, warning).
    Forward(Value),
    /// The kernel finished request `id`; `memory` is its resident bytes, if sent.
    Done { memory: Option<u64> },
    /// Not ours (another request id, a pong) or blank.
    Ignore,
}

/// Route one stdout line (trailing newline optional).
///
/// A protocol line for `id` is forwarded (or ends the evaluation when it is
/// `done`). A line that is NOT JSON is text the kernel printed outside the
/// protocol -- `Print` output from a kernel that predates the "stream" message,
/// or anything writing to stdout directly -- and is forwarded as a `stream`
/// message so it appears in the cell, where it used to be logged and dropped.
pub(crate) fn classify_stdout_line(line: &str, id: u32) -> LineAction {
    let line = line.trim_end_matches(['\n', '\r']);
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return LineAction::Ignore;
    }
    match serde_json::from_str::<Value>(trimmed) {
        Ok(msg) if msg.is_object() => {
            let for_us = msg.get("id").and_then(|v| v.as_u64()) == Some(id as u64);
            if !for_us {
                return LineAction::Ignore;
            }
            if msg["type"] == "done" {
                return LineAction::Done {
                    memory: msg.get("memory").and_then(|v| v.as_u64()),
                };
            }
            LineAction::Forward(msg)
        }
        _ => {
            LineAction::Forward(json!({ "id": id, "type": "stream", "text": format!("{line}\n") }))
        }
    }
}

/// A stderr chunk seen during evaluation `id`, as a `message` for the cell.
/// None when it is only whitespace.
pub(crate) fn stderr_message(text: &str, id: u32) -> Option<Value> {
    let text = text.trim_end();
    if text.trim().is_empty() {
        return None;
    }
    Some(json!({ "id": id, "type": "message", "text": text }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn request_is_one_line_with_the_cell_flag() {
        let req = cell_request(4, "a = 1\nPrint[\"x\"]").unwrap();
        assert!(req.ends_with('\n'));
        assert_eq!(
            req.matches('\n').count(),
            1,
            "the request must be ONE line: {req}"
        );
        let v: Value = serde_json::from_str(req.trim_end()).unwrap();
        assert_eq!(
            v,
            json!({ "id": 4, "expr": "a = 1\nPrint[\"x\"]", "cell": true })
        );
    }

    #[test]
    fn result_for_this_request_is_forwarded() {
        let line = r#"{"id":3,"type":"expr","payload":"2","latex":"2"}"#;
        assert_eq!(
            classify_stdout_line(line, 3),
            LineAction::Forward(json!({"id":3,"type":"expr","payload":"2","latex":"2"}))
        );
    }

    #[test]
    fn stream_and_message_kinds_are_forwarded_in_order_before_done() {
        let lines = [
            r#"{"id":7,"type":"stream","text":"hello\n"}"#,
            r#"{"id":7,"type":"message","text":"Power::infy: Infinite expression 1/0 encountered."}"#,
            r#"{"id":7,"type":"expr","payload":"ComplexInfinity"}"#,
            r#"{"id":7,"type":"done","memory":1234}"#,
        ];
        let got: Vec<LineAction> = lines.iter().map(|l| classify_stdout_line(l, 7)).collect();
        assert!(matches!(&got[0], LineAction::Forward(m) if m["type"] == "stream"));
        assert!(matches!(&got[1], LineAction::Forward(m) if m["type"] == "message"));
        assert!(matches!(&got[2], LineAction::Forward(m) if m["type"] == "expr"));
        assert_eq!(got[3], LineAction::Done { memory: Some(1234) });
    }

    #[test]
    fn done_without_memory() {
        assert_eq!(
            classify_stdout_line(r#"{"id":1,"type":"done"}"#, 1),
            LineAction::Done { memory: None }
        );
    }

    #[test]
    fn other_ids_pongs_and_blank_lines_are_ignored() {
        assert_eq!(
            classify_stdout_line(r#"{"id":2,"type":"expr","payload":"x"}"#, 1),
            LineAction::Ignore
        );
        assert_eq!(
            classify_stdout_line(r#"{"id":2,"type":"done"}"#, 1),
            LineAction::Ignore
        );
        assert_eq!(
            classify_stdout_line(r#"{"type":"pong"}"#, 1),
            LineAction::Ignore
        );
        assert_eq!(classify_stdout_line("   \n", 1), LineAction::Ignore);
    }

    #[test]
    fn non_json_stdout_becomes_stream_text() {
        // Print output from a kernel that predates the "stream" message.
        assert_eq!(
            classify_stdout_line("hello world\n", 5),
            LineAction::Forward(json!({"id":5,"type":"stream","text":"hello world\n"}))
        );
        // A bare JSON scalar is not a protocol message either.
        assert_eq!(
            classify_stdout_line("42", 5),
            LineAction::Forward(json!({"id":5,"type":"stream","text":"42\n"}))
        );
    }

    /// End to end against the real kernel, when it is built: the request this
    /// module writes and the routing it applies, over the actual pipe protocol.
    /// Skipped (passes) when ../../Mathilda does not exist.
    #[test]
    fn real_kernel_cell_round_trip() {
        use std::io::Write;
        use std::process::{Command, Stdio};
        let bin = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../Mathilda");
        if !bin.exists() {
            eprintln!("skipping: {} not built", bin.display());
            return;
        }
        let mut child = Command::new(&bin)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .env("MATHILDA_NO_WINDOW", "1")
            .spawn()
            .expect("spawn kernel");
        {
            let stdin = child.stdin.as_mut().unwrap();
            let cell = "x = 2;\nPrint[\"hi\"]\n1/0\nx + 1";
            stdin
                .write_all(cell_request(1, cell).unwrap().as_bytes())
                .unwrap();
            stdin.write_all(b"{\"type\":\"quit\"}\n").unwrap();
        }
        let out = child.wait_with_output().expect("kernel output");
        let mut seen = Vec::new();
        for line in String::from_utf8_lossy(&out.stdout).lines() {
            match classify_stdout_line(line, 1) {
                LineAction::Forward(m) => seen.push(m),
                LineAction::Done { .. } => break,
                LineAction::Ignore => {}
            }
        }
        let kinds: Vec<&str> = seen
            .iter()
            .map(|m| m["type"].as_str().unwrap_or(""))
            .collect();
        // A "line" line precedes each of the four statements: it carries the
        // kernel's $Line, which is what `%` / In[n] / Out[n] resolve against and
        // what the front end labels the cell with (see src/repl.c).
        assert_eq!(
            kinds,
            [
                "line", // x = 2;   -- suppressed, but still a numbered line
                "line", "stream", // Print["hi"]
                "line", "message", "expr", // 1/0
                "line", "expr",   // x + 1
            ],
            "{seen:?}"
        );
        let of_type = |t: &str| -> Vec<&Value> {
            seen.iter().filter(|m| m["type"] == t).collect()
        };
        assert_eq!(of_type("stream")[0]["text"], "hi\n");
        assert!(of_type("message")[0]["text"]
            .as_str()
            .unwrap()
            .starts_with("Power::infy"));
        assert_eq!(of_type("expr")[1]["payload"], "3");
        /* The lines are consecutive and start at 1: one per statement, which is
         * how `%` means "the previous result" inside a cell as well as across
         * cells. */
        let lines: Vec<u64> = of_type("line")
            .iter()
            .map(|m| m["line"].as_u64().unwrap())
            .collect();
        assert_eq!(lines, [1, 2, 3, 4], "{seen:?}");
    }

    #[test]
    fn stderr_becomes_a_message() {
        assert_eq!(
            stderr_message("General::stop: Further output suppressed.\n", 9),
            Some(
                json!({"id":9,"type":"message","text":"General::stop: Further output suppressed."})
            )
        );
        assert_eq!(stderr_message(" \n", 9), None);
    }
}
