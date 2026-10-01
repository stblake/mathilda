// commands.rs — Tauri commands exposed to the Svelte frontend

use crate::kernel::MathildaKernel;
use crate::notebook_format::{parse_stanzas, serialize_stanzas};
use serde_json::Value;
use tauri::ipc::Channel;
use tauri::State;

/// Evaluate a Mathilda expression, streaming output messages through
/// `channel` until the kernel emits "done".
#[tauri::command]
pub async fn evaluate_cell(
    expr: String,
    channel: Channel<Value>,
    kernel: State<'_, MathildaKernel>,
) -> Result<(), String> {
    kernel.evaluate(expr, channel).await
}

/// Restart the Mathilda kernel process (used for "Run All" fresh context).
#[tauri::command]
pub async fn restart_kernel(kernel: State<'_, MathildaKernel>) -> Result<(), String> {
    kernel.restart().await
}

/// Interrupt the current computation by sending SIGINT to the kernel.
#[tauri::command]
pub async fn interrupt_kernel(kernel: State<'_, MathildaKernel>) -> Result<(), String> {
    kernel.interrupt().await
}

/// Ping the kernel — returns Ok(()) if alive.
#[tauri::command]
pub async fn ping_kernel(kernel: State<'_, MathildaKernel>) -> Result<(), String> {
    kernel.ping().await
}

/// Evaluate one expression quietly (no cell history) and return {payload, latex,
/// error}. Used by the output "Convert To" menu (FullForm / TeXForm).
#[tauri::command]
pub async fn eval_once(
    expr: String,
    kernel: State<'_, MathildaKernel>,
) -> Result<Value, String> {
    kernel.eval_once(expr).await
}

/// Parse an expression for structural (bottom-up) selection; returns [start, end]
/// byte-span pairs for every subexpression. Does not evaluate.
#[tauri::command]
pub async fn syntax_spans(
    expr: String,
    kernel: State<'_, MathildaKernel>,
) -> Result<Vec<[i64; 2]>, String> {
    kernel.fetch_spans(expr).await
}

/// Save one notebook to a `.mnb` file (see `notebook_format.rs`).
/// `cells` is a JSON array of objects: [{type, source}, ...].
/// Only type and source are written; outputs are ephemeral.
#[tauri::command]
pub async fn save_notebook(path: String, cells: Vec<Value>) -> Result<(), String> {
    std::fs::write(&path, serialize_stanzas(&cells)).map_err(|e| format!("save: {e}"))
}

/// Load a notebook from a `.mnb` file.
/// Returns a JSON array of cell objects: [{type, source}, ...].
#[tauri::command]
pub async fn load_notebook(path: String) -> Result<Vec<Value>, String> {
    let content = std::fs::read_to_string(&path).map_err(|e| format!("load: {e}"))?;
    Ok(parse_stanzas(&content))
}

/// Save a library JSON blob to a .lb file.
/// `json` is the full serialized library string (produced by `serializeLibrary` in canvas.ts).
#[tauri::command]
pub async fn save_library(path: String, json: String) -> Result<(), String> {
    std::fs::write(&path, &json).map_err(|e| format!("save_library: {e}"))
}

/// Load a library from a .lb file.
/// Returns the raw JSON string for parsing on the frontend.
#[tauri::command]
pub async fn load_library(path: String) -> Result<String, String> {
    std::fs::read_to_string(&path).map_err(|e| format!("load_library: {e}"))
}

/// Set the native OS window title from the frontend.
/// More reliable than the JS getCurrentWindow().setTitle() in WebKit.
#[tauri::command]
pub async fn set_window_title(
    app: tauri::AppHandle,
    title: String,
) -> Result<(), String> {
    use tauri::Manager;
    if let Some(win) = app.get_webview_window("main") {
        win.set_title(&title).map_err(|e| e.to_string())
    } else {
        Err("Window 'main' not found".into())
    }
}

// ---------------------------------------------------------------------------
// File > Open Recent (see recent.rs)
//
// These take an `AppHandle` and look the state up, rather than declaring
// `State<'_, RecentFiles>`: `tauri::menu` — and so `RecentFiles` — is desktop-only, and a
// desktop-only parameter type would force `generate_handler!` into two divergent copies of
// the whole command list, one per platform, which is exactly the kind of list that drifts.
// With the gate INSIDE each body the signature is platform-independent, there is one handler
// list, and on mobile (no menu bar to serve) each call is a no-op returning nothing.
//
// `try_state` rather than `state`: the latter panics when unmanaged, and a panic in the file
// path would be a worse bug than a recent list that quietly does not update.

/// The recent-files list, newest first, as absolute paths.
///
/// The front end indexes into this with the `recent-<i>` menu id it was clicked with, so the
/// order here IS the order of the menu items.
#[tauri::command]
pub async fn recent_files(app: tauri::AppHandle) -> Result<Vec<String>, String> {
    #[cfg(desktop)]
    {
        use tauri::Manager;
        return Ok(app
            .try_state::<crate::recent::RecentFiles>()
            .map(|r| r.list())
            .unwrap_or_default());
    }
    #[cfg(not(desktop))]
    {
        let _ = app;
        Ok(Vec::new())
    }
}

/// Record a file that was just opened or saved.
#[tauri::command]
pub async fn push_recent_file(app: tauri::AppHandle, path: String) -> Result<(), String> {
    #[cfg(desktop)]
    {
        use tauri::Manager;
        if let Some(r) = app.try_state::<crate::recent::RecentFiles>() {
            return r.push(std::path::Path::new(&path));
        }
        return Ok(());
    }
    #[cfg(not(desktop))]
    {
        let _ = (app, path);
        Ok(())
    }
}

/// Drop a file from the list — the front end calls this when opening one failed, so a moved
/// or deleted file leaves the menu instead of sitting there failing.
#[tauri::command]
pub async fn forget_recent_file(app: tauri::AppHandle, path: String) -> Result<(), String> {
    #[cfg(desktop)]
    {
        use tauri::Manager;
        if let Some(r) = app.try_state::<crate::recent::RecentFiles>() {
            return r.forget(std::path::Path::new(&path));
        }
        return Ok(());
    }
    #[cfg(not(desktop))]
    {
        let _ = (app, path);
        Ok(())
    }
}

/// Empty the list — File > Open Recent > Clear Menu.
#[tauri::command]
pub async fn clear_recent_files(app: tauri::AppHandle) -> Result<(), String> {
    #[cfg(desktop)]
    {
        use tauri::Manager;
        if let Some(r) = app.try_state::<crate::recent::RecentFiles>() {
            return r.clear();
        }
        return Ok(());
    }
    #[cfg(not(desktop))]
    {
        let _ = app;
        Ok(())
    }
}
