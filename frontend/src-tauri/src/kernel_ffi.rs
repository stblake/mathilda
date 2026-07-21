//! kernel_ffi.rs — in-process Mathilda kernel for mobile targets.
//!
//! Drop-in replacement for the sidecar-based `kernel.rs` used on desktop. It
//! exposes the identical public API (`empty`/`start`/`evaluate`/`ping`/
//! `restart`/`interrupt`/`is_running`) so `commands.rs` and `lib.rs` compile
//! unchanged, but runs the kernel via FFI (`crate::ffi`) instead of spawning a
//! child process — mandatory on iOS/Android where process spawning is blocked.
//!
//! The C kernel is single-threaded and non-reentrant, so every evaluation is
//! serialized behind a mutex and executed on a blocking thread (the eval can be
//! long-running; we must not block the async runtime's worker).

use crate::ffi;
use serde_json::{json, Value};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use tauri::ipc::Channel;
use tauri::Manager;
use tokio::sync::Mutex;

// The internal module tree (init.m + derivative/integral tables) embedded at
// compile time — extracted to the app data dir at startup on Android, where
// APK-bundled assets are not reachable via the C loader's fopen(). Path is
// relative to this crate (src-tauri): ../../src/internal is the repo tree.
#[cfg(target_os = "android")]
use include_dir::{include_dir, Dir};
#[cfg(target_os = "android")]
static INTERNAL_TREE: Dir<'_> = include_dir!("$CARGO_MANIFEST_DIR/../../src/internal");

/// Handle stored in Tauri State — owns the in-process kernel lifecycle.
pub struct MathildaKernel {
    /// Serializes access to the non-reentrant C kernel.
    lock: Arc<Mutex<()>>,
    ready: Arc<AtomicBool>,
    app: tauri::AppHandle,
}

impl MathildaKernel {
    /// Create a handle synchronously; the kernel is initialized in `start()`.
    pub fn empty(app: tauri::AppHandle) -> Self {
        Self {
            lock: Arc::new(Mutex::new(())),
            ready: Arc::new(AtomicBool::new(false)),
            app,
        }
    }

    /// Convenience: build and start in one step.
    pub async fn new(app: tauri::AppHandle) -> Result<Self, String> {
        let k = Self::empty(app);
        k.start().await?;
        Ok(k)
    }

    /// Locate the bundled `internal/` module tree (init.m + tables) inside the
    /// app bundle's resources, so the kernel finds its bootstrap on a device
    /// where there is no source tree and no writable CWD. Tauri places declared
    /// resources under `<resources>/assets/…`, but the exact prefix has varied
    /// across versions/platforms, so probe the known candidates and return the
    /// first that actually contains `init.m`.
    #[cfg(not(target_os = "android"))]
    fn home_dir(&self) -> Option<String> {
        let res = self.app.path().resource_dir().ok()?;
        for sub in ["assets/internal", "internal", "_up_/src/internal"] {
            let cand = res.join(sub);
            if cand.join("init.m").is_file() {
                return Some(cand.to_string_lossy().into_owned());
            }
        }
        // Fall back to the conventional location even if the probe missed, so
        // the loader emits a diagnostic path rather than nothing.
        Some(res.join("assets/internal").to_string_lossy().into_owned())
    }

    /// Android bundles app resources inside the APK, which the C loader cannot
    /// `fopen()`. So the `internal/` module tree (init.m + tables, ~160 KB) is
    /// embedded in the binary at compile time and extracted once to the app's
    /// writable data dir; MATHILDA_HOME then points at that real path.
    #[cfg(target_os = "android")]
    fn home_dir(&self) -> Option<String> {
        let base = self.app.path().app_data_dir().ok()?;
        let dst = base.join("mathilda-internal");
        if !dst.join("init.m").is_file() {
            let _ = std::fs::create_dir_all(&dst);
            if let Err(e) = INTERNAL_TREE.extract(&dst) {
                log::error!("failed to extract internal tree to {dst:?}: {e}");
            }
        }
        Some(dst.to_string_lossy().into_owned())
    }

    /// Initialize the kernel on a blocking thread (idempotent). Sets
    /// MATHILDA_HOME from the bundle resources first so init.m resolves.
    pub async fn start(&self) -> Result<(), String> {
        if self.ready.load(Ordering::Acquire) {
            return Ok(());
        }
        let home = self.home_dir();
        let ready = self.ready.clone();
        let lock = self.lock.clone();
        tauri::async_runtime::spawn_blocking(move || {
            // Hold the lock across init so no eval races the symbol table.
            let _guard = lock.blocking_lock();
            if let Some(h) = home {
                ffi::set_home(&h);
            }
            ffi::init();
            ready.store(true, Ordering::Release);
        })
        .await
        .map_err(|e| format!("kernel init join error: {e}"))?;
        log::info!("Mathilda in-process kernel ready (v{})", ffi::version());
        Ok(())
    }

    /// Readiness probe (kept for API parity with the sidecar kernel).
    pub async fn ping(&self) -> Result<(), String> {
        if self.ready.load(Ordering::Acquire) {
            Ok(())
        } else {
            Err("kernel not initialized".into())
        }
    }

    /// Evaluate `expr` in-process and stream one message to `channel`. The
    /// desktop kernel streams many NDJSON lines then a "done"; here the C API
    /// returns a single structured result, so we emit one message and return
    /// (the frontend treats the resolved promise as "done").
    pub async fn evaluate(&self, expr: String, channel: Channel<Value>) -> Result<(), String> {
        if !self.ready.load(Ordering::Acquire) {
            // Lazily initialize if evaluate arrives before start() completed.
            self.start().await?;
        }
        let lock = self.lock.clone();
        // One structured eval (not eval + eval_latex): a single evaluation keeps
        // side effects from running twice AND carries the Plotly payload for
        // Graphics/Graphics3D, so Plot[...] renders as a chart instead of raw
        // red Graphics[...] text.
        let raw = tauri::async_runtime::spawn_blocking(move || {
            let _guard = lock.blocking_lock();
            ffi::eval_json(&expr)
        })
        .await
        .map_err(|e| format!("eval join error: {e}"))?;

        // The kernel already emitted {"type":...,"payload"/"latex"/"message":...};
        // stamp the request id and forward it verbatim.
        let mut msg = serde_json::from_str::<Value>(&raw).unwrap_or_else(|_| {
            json!({ "type": "error", "message": format!("bad kernel output: {raw}") })
        });
        if let Some(obj) = msg.as_object_mut() {
            obj.insert("id".into(), json!(0));
        }
        channel.send(msg).map_err(|e| format!("channel: {e}"))?;
        Ok(())
    }

    /// No separate process to kill; state is process-global. No-op for parity.
    pub async fn kill(&self) -> Result<(), String> {
        Ok(())
    }

    /// The in-process kernel has no resettable session yet (the C symbol table
    /// persists for the app's lifetime). Report success so "Run All" proceeds;
    /// a true reset would need a `mathilda_ffi_reset()` entry point.
    pub async fn restart(&self) -> Result<(), String> {
        log::warn!("restart requested: in-process kernel keeps global state (no-op)");
        Ok(())
    }

    /// Synchronous C evaluations cannot be interrupted from another thread.
    pub async fn interrupt(&self) -> Result<(), String> {
        log::warn!("interrupt requested: not supported for in-process kernel (no-op)");
        Ok(())
    }

    pub async fn is_running(&self) -> bool {
        self.ready.load(Ordering::Acquire)
    }
}
