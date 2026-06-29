mod commands;
mod kernel;

use commands::{evaluate_cell, interrupt_kernel, load_library, load_notebook, ping_kernel, restart_kernel, save_library, save_notebook, set_window_title};
use kernel::MathildaKernel;
use tauri::menu::{Menu, MenuItem, PredefinedMenuItem, Submenu};
use tauri::{Emitter, Manager};

fn build_menu(app: &tauri::App) -> tauri::Result<Menu<tauri::Wry>> {
    let file = Submenu::with_items(
        app,
        "File",
        true,
        &[
            &MenuItem::with_id(app, "open",    "Open…",    true, Some("CmdOrCtrl+O"))?,
            &MenuItem::with_id(app, "save",    "Save",     true, Some("CmdOrCtrl+S"))?,
            &MenuItem::with_id(app, "save-as", "Save As…", true, Some("CmdOrCtrl+Shift+S"))?,
            &PredefinedMenuItem::separator(app)?,
            &PredefinedMenuItem::close_window(app, None)?,
        ],
    )?;

    let edit = Submenu::with_items(
        app,
        "Edit",
        true,
        &[
            &MenuItem::with_id(app, "add-cell",     "Add Cell Below",   true, Some("CmdOrCtrl+B"))?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "copy-cells",   "Copy Cell(s)",     true, Some("CmdOrCtrl+C"))?,
            &MenuItem::with_id(app, "paste-cells",  "Paste Cell(s)",    true, Some("CmdOrCtrl+V"))?,
            &MenuItem::with_id(app, "delete-cells", "Delete Cell(s)",   true, Some("Backspace"))?,
        ],
    )?;

    let kernel_menu = Submenu::with_items(
        app,
        "Kernel",
        true,
        &[
            &MenuItem::with_id(app, "run-all",   "Run All",         true, Some("CmdOrCtrl+Shift+Return"))?,
            &MenuItem::with_id(app, "restart",   "Restart Kernel",  true, Some("CmdOrCtrl+Shift+R"))?,
            &MenuItem::with_id(app, "interrupt", "Interrupt",       true, Some("CmdOrCtrl+Period"))?,
        ],
    )?;

    let view = Submenu::with_items(
        app,
        "View",
        true,
        &[
            &MenuItem::with_id(app, "toggle-dark", "Toggle Dark Mode", true, Some("CmdOrCtrl+Shift+D"))?,
        ],
    )?;

    Menu::with_items(app, &[
        &Submenu::with_items(app, "Mathilda", true, &[
            &PredefinedMenuItem::about(app, None, None)?,
            &PredefinedMenuItem::separator(app)?,
            &PredefinedMenuItem::services(app, None)?,
            &PredefinedMenuItem::separator(app)?,
            &PredefinedMenuItem::hide(app, None)?,
            &PredefinedMenuItem::hide_others(app, None)?,
            &PredefinedMenuItem::separator(app)?,
            &PredefinedMenuItem::quit(app, None)?,
        ])?,
        &file,
        &edit,
        &kernel_menu,
        &view,
    ])
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_log::Builder::default().build())
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_fs::init())
        .setup(|app| {
            // Native menu
            let menu = build_menu(app)?;
            app.set_menu(menu)?;
            app.on_menu_event(|app, event| {
                let id = event.id().as_ref().to_string();
                let _ = app.emit(&format!("menu:{id}"), ());
            });

            // Kernel — managed synchronously, spawned async
            let kernel = MathildaKernel::empty(app.handle().clone());
            app.manage(kernel);
            let handle = app.handle().clone();
            tauri::async_runtime::spawn(async move {
                let kernel = handle.state::<MathildaKernel>();
                match kernel.start().await {
                    Ok(()) => log::info!("Mathilda kernel ready"),
                    Err(e) => log::error!("Kernel failed to start: {e}"),
                }
            });

            Ok(())
        })
        .invoke_handler(tauri::generate_handler![
            evaluate_cell,
            restart_kernel,
            interrupt_kernel,
            ping_kernel,
            save_notebook,
            load_notebook,
            save_library,
            load_library,
            set_window_title,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
