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
