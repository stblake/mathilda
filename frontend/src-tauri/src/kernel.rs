// kernel.rs — manages the Mathilda CAS child process

use serde_json::Value;
use std::sync::Arc;
use std::time::Duration;
use tauri::async_runtime::Receiver;
use tauri::ipc::Channel;
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
            *guard = Some(KernelState { child, rx, next_id: 1 });
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

        let expr_json = serde_json::to_string(&expr).map_err(|e| e.to_string())?;
        let request = format!("{{\"id\":{id},\"expr\":{expr_json}}}\n");
        state
            .child
            .write(request.as_bytes())
            .map_err(|e| format!("write: {e}"))?;

        loop {
            match state.rx.recv().await {
                Some(CommandEvent::Stdout(bytes)) => {
                    let line = String::from_utf8_lossy(&bytes);
                    let trimmed = line.trim();
                    if trimmed.is_empty() {
                        continue;
                    }
                    match serde_json::from_str::<Value>(trimmed) {
                        Ok(msg) => {
                            let for_us = msg
                                .get("id")
                                .and_then(|v| v.as_u64())
                                .map(|v| v == id as u64)
                                .unwrap_or(false);
                            if !for_us {
                                continue;
                            }
                            if msg["type"] == "done" {
                                break;
                            }
                            channel.send(msg).map_err(|e| format!("channel: {e}"))?;
                        }
                        Err(e) => {
                            log::warn!("Parse error on kernel output: {e} — {trimmed}");
                        }
                    }
                }
                Some(CommandEvent::Terminated(p)) => {
                    return Err(format!("Kernel died (code {:?})", p.code));
                }
                Some(CommandEvent::Stderr(bytes)) => {
                    log::debug!("[kernel] {}", String::from_utf8_lossy(&bytes));
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
