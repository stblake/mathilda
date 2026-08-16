mod commands;

// Kernel backend is chosen at compile time:
//   * desktop  -> `kernel.rs`: spawns the `mathilda` sidecar over stdio.
//   * mobile   -> `kernel_ffi.rs`: runs the kernel in-process via FFI, because
//                 iOS/Android sandboxes forbid spawning child processes.
// Both expose an identical `MathildaKernel` API, so the rest of the app is
// backend-agnostic.
#[cfg(not(mobile))]
mod kernel;
#[cfg(mobile)]
mod ffi;
#[cfg(mobile)]
#[path = "kernel_ffi.rs"]
mod kernel;

use commands::{evaluate_cell, interrupt_kernel, load_library, load_notebook, ping_kernel, restart_kernel, save_library, save_notebook, set_window_title};
use kernel::MathildaKernel;
#[cfg(desktop)]
use tauri::menu::{Menu, MenuItem, PredefinedMenuItem, Submenu};
#[cfg(desktop)]
use tauri::Emitter;
use tauri::Manager;

// Native menu bar is a desktop-only concept; iOS/Android have no app menu.
#[cfg(desktop)]
fn build_menu(app: &tauri::App) -> tauri::Result<Menu<tauri::Wry>> {
    // Item ids are a CONTRACT with the webview: each one is emitted as `menu:<id>` and handled by
    // runMenuCommand in src/lib/menuCommands.ts, whose MENU_IDS list is what App.svelte subscribes
    // to. An id added here without a case there is a menu item that does nothing, which is why the
    // dispatcher warns on an unknown id rather than ignoring it.
    let file = Submenu::with_items(
        app,
        "File",
        true,
        &[
            &MenuItem::with_id(app, "file.new",  "New Notebook", true, Some("CmdOrCtrl+N"))?,
            &MenuItem::with_id(app, "open",      "Open…",        true, Some("CmdOrCtrl+O"))?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "save",      "Save",         true, Some("CmdOrCtrl+S"))?,
            &MenuItem::with_id(app, "save-as",   "Save As…",     true, Some("CmdOrCtrl+Shift+S"))?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "file.print", "Print…",      true, Some("CmdOrCtrl+P"))?,
            &PredefinedMenuItem::separator(app)?,
            // Two different closes, deliberately both present: one closes the focused NOTEBOOK,
            // the other the window. Cmd+W goes to the notebook because that is the one a reader
            // reaches for repeatedly.
            &MenuItem::with_id(app, "file.close", "Close Notebook", true, Some("CmdOrCtrl+W"))?,
            &PredefinedMenuItem::close_window(app, Some("Close Window"))?,
        ],
    )?;

    // Use the PREDEFINED clipboard/undo items so the standard shortcuts (Cmd+C/X/V/A, Cmd+Z) route
    // through the macOS responder chain to the focused text editor. Binding custom items to those
    // accelerators would hijack the keys app-wide and break copy/paste inside cells. The items
    // below are the ones with no native equivalent, so they carry ids of their own.
    let edit = Submenu::with_items(
        app,
        "Edit",
        true,
        &[
            &PredefinedMenuItem::undo(app, None)?,
            &PredefinedMenuItem::redo(app, None)?,
            &PredefinedMenuItem::separator(app)?,
            &PredefinedMenuItem::cut(app, None)?,
            &PredefinedMenuItem::copy(app, None)?,
            &PredefinedMenuItem::paste(app, None)?,
            &PredefinedMenuItem::select_all(app, None)?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "edit.comment", "Un/Comment Selection", true,
                               Some("CmdOrCtrl+/"))?,
            &MenuItem::with_id(app, "edit.indent",  "Indent Selected Lines",  true, None::<&str>)?,
            &MenuItem::with_id(app, "edit.outdent", "Outdent Selected Lines", true, None::<&str>)?,
            &MenuItem::with_id(app, "edit.dupLine", "Duplicate Line", true,
                               Some("CmdOrCtrl+Shift+L"))?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "edit.findDoc", "Documentation for Selection", true,
                               Some("CmdOrCtrl+Shift+F"))?,
        ],
    )?;

    let insert = Submenu::with_items(
        app,
        "Insert",
        true,
        &[
            &MenuItem::with_id(app, "insert.code",    "Input Cell",   true, Some("CmdOrCtrl+B"))?,
            &MenuItem::with_id(app, "insert.text",    "Text Cell",    true, None::<&str>)?,
            &MenuItem::with_id(app, "insert.section", "Section Cell", true, None::<&str>)?,
        ],
    )?;

    // Convert To is flat rather than a nested submenu: three items read better than a submenu
    // holding three, and macOS submenus cost a second gesture to reach.
    let cell = Submenu::with_items(
        app,
        "Cell",
        true,
        &[
            &MenuItem::with_id(app, "cell.toInput",   "Convert to Input",   true, None::<&str>)?,
            &MenuItem::with_id(app, "cell.toText",    "Convert to Text",    true, None::<&str>)?,
            &MenuItem::with_id(app, "cell.toSection", "Convert to Section", true, None::<&str>)?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "cell.divide",    "Divide Cell", true,
                               Some("CmdOrCtrl+Shift+D"))?,
            &MenuItem::with_id(app, "cell.merge",     "Merge Cells", true,
                               Some("CmdOrCtrl+Shift+M"))?,
            &MenuItem::with_id(app, "cell.duplicate", "Duplicate Cell", true, None::<&str>)?,
            &MenuItem::with_id(app, "cell.delete",    "Delete Cell",    true, None::<&str>)?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "cell.clearOutput",    "Delete Output",     true, None::<&str>)?,
            &MenuItem::with_id(app, "cell.clearAllOutput", "Delete All Output", true, None::<&str>)?,
        ],
    )?;

    let evaluation = Submenu::with_items(
        app,
        "Evaluation",
        true,
        &[
            &MenuItem::with_id(app, "eval.cell", "Evaluate Cell", true,
                               Some("CmdOrCtrl+Return"))?,
            &MenuItem::with_id(app, "run-all",   "Evaluate Notebook", true,
                               Some("CmdOrCtrl+Shift+Return"))?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "interrupt", "Abort Evaluation", true,
                               Some("CmdOrCtrl+Period"))?,
            &MenuItem::with_id(app, "restart",   "Restart Kernel",   true,
                               Some("CmdOrCtrl+Shift+R"))?,
        ],
    )?;

    // Graphics is a documentation entry point rather than a set of editing commands: the renderer
    // has no interactive object model to act on, so greyed drawing tools would be noise.
    let graphics = Submenu::with_items(
        app,
        "Graphics",
        true,
        &[
            &MenuItem::with_id(app, "gfx.plot",     "Plot Documentation",     true, None::<&str>)?,
            &MenuItem::with_id(app, "gfx.image",    "Image Documentation",    true, None::<&str>)?,
            &MenuItem::with_id(app, "gfx.image3d",  "Image3D Documentation",  true, None::<&str>)?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "gfx.graphics", "Graphics Documentation", true, None::<&str>)?,
        ],
    )?;

    let view = Submenu::with_items(
        app,
        "View",
        true,
        &[
            &MenuItem::with_id(app, "toggle-dark", "Toggle Dark Mode", true,
                               Some("CmdOrCtrl+Shift+T"))?,
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
        &insert,
        &cell,
        &evaluation,
        &graphics,
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
            // Native menu (desktop only — mobile has no app menu bar).
            #[cfg(desktop)]
            {
                let menu = build_menu(app)?;
                app.set_menu(menu)?;
                app.on_menu_event(|app, event| {
                    let id = event.id().as_ref().to_string();
                    let _ = app.emit(&format!("menu:{id}"), ());
                });
            }

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
