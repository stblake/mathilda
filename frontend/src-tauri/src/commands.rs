// commands.rs — Tauri commands exposed to the Svelte frontend

use crate::kernel::MathildaKernel;
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

/// Save notebook source to a file (plain-text .mathilda format).
/// `cells` is a JSON array of objects: [{type, source}, ...].
/// Only the source field is written; outputs are ephemeral.
#[tauri::command]
pub async fn save_notebook(path: String, cells: Vec<Value>) -> Result<(), String> {
    let mut out = String::new();
    for (i, cell) in cells.iter().enumerate() {
        let cell_type = cell["type"].as_str().unwrap_or("code");
        let source = cell["source"].as_str().unwrap_or("");
        // Stanza format:
        //   (* cell: code *)
        //   <source>
        //
        // Blank line between cells; compatible with Mathilda's .m format.
        if i > 0 {
            out.push('\n');
        }
        out.push_str(&format!("(* cell: {cell_type} *)\n"));
        out.push_str(source);
        if !source.ends_with('\n') {
            out.push('\n');
        }
    }
    std::fs::write(&path, &out).map_err(|e| format!("save: {e}"))
}

/// Load a notebook from a .mathilda file.
/// Returns a JSON array of cell objects: [{type, source}, ...].
#[tauri::command]
pub async fn load_notebook(path: String) -> Result<Vec<Value>, String> {
    let content = std::fs::read_to_string(&path).map_err(|e| format!("load: {e}"))?;
    let mut cells: Vec<Value> = Vec::new();
    let mut current_type = "code".to_string();
    let mut current_source = String::new();
    let mut in_cell = false;

    for line in content.lines() {
        if let Some(rest) = line.strip_prefix("(* cell:") {
            // Save the previous cell if any.
            if in_cell {
                let src = current_source.trim_end_matches('\n').to_string();
                cells.push(serde_json::json!({
                    "type": current_type,
                    "source": src,
                }));
                current_source.clear();
            }
            // Parse cell type from "(* cell: code *)" etc.
            current_type = rest
                .trim()
                .trim_end_matches("*)")
                .trim()
                .to_string();
            in_cell = true;
        } else if in_cell {
            current_source.push_str(line);
            current_source.push('\n');
        }
    }
    // Flush last cell.
    if in_cell {
        let src = current_source.trim_end_matches('\n').to_string();
        cells.push(serde_json::json!({
            "type": current_type,
            "source": src,
        }));
    }

    // If file has no stanza markers, treat entire file as a single code cell.
    if cells.is_empty() && !content.trim().is_empty() {
        cells.push(serde_json::json!({
            "type": "code",
            "source": content.trim_end_matches('\n'),
        }));
    }

    Ok(cells)
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
