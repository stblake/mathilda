mod commands;
mod notebook_format;
// File > Open Recent: the list, its store, and the native submenu. Desktop-only, because it
// exists to serve a menu bar and mobile has none.
#[cfg(desktop)]
mod recent;

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

use commands::{evaluate_cell, interrupt_kernel, load_library, load_notebook, ping_kernel, restart_kernel, save_library, save_notebook, set_window_title, syntax_spans, eval_once};
use commands::{clear_recent_files, forget_recent_file, push_recent_file, recent_files};
use commands::{cancel_quit, confirm_and_quit, sync_view_menu, QuitGuard};
use std::sync::atomic::Ordering;
use kernel::MathildaKernel;
#[cfg(desktop)]
use recent::RecentFiles;
#[cfg(desktop)]
use tauri::menu::{CheckMenuItem, Menu, MenuItem, PredefinedMenuItem, Submenu};
use tauri::Emitter;
use tauri::Manager;

/// Handles for the View-menu items whose checkmark mirrors a webview store, so
/// `sync_view_menu` (commands.rs) can update them when the store changes. Managed
/// in app state; desktop-only, like the native menu itself.
#[cfg(desktop)]
pub struct ViewMenu {
    pub autocomplete: CheckMenuItem<tauri::Wry>,
    /// (colour-scheme id, its menu item); exactly one is checked at a time.
    pub schemes: Vec<(String, CheckMenuItem<tauri::Wry>)>,
}

// Native menu bar is a desktop-only concept; iOS/Android have no app menu.
//
// Returns the menu together with the Open Recent state, because that submenu's items have to
// be reachable after startup to redraw the list and there is no way to look a `MenuItem` back
// up out of an installed `Menu`. `setup` installs the menu, draws the list once, and hands the
// state to `app.manage`.
#[cfg(desktop)]
fn build_menu(app: &tauri::App) -> tauri::Result<(Menu<tauri::Wry>, RecentFiles, ViewMenu)> {
    // Item ids are a CONTRACT with the webview: each one is emitted as `menu:<id>` and handled by
    // runMenuCommand in src/lib/menuCommands.ts, whose MENU_IDS list is what App.svelte subscribes
    // to. An id added here without a case there is a menu item that does nothing, which is why the
    // dispatcher warns on an unknown id rather than ignoring it.

    // Open Recent. The ten slots are spelled out one per line, and deliberately so: the ids are
    // read out of THIS FILE as source text by tools/check_menu_ids.py and diffed against the
    // handler cases in menuCommands.ts, so a loop building `format!("recent-{i}")` would leave the
    // Rust half of that check blind and the TypeScript half looking like dead wiring. Slot i shows
    // recent-list entry i; RecentFiles derives its cap from how many slots it is given, so the
    // number ten is not written down anywhere else.
    let recent_slots = vec![
        MenuItem::with_id(app, "recent-0", "", true, None::<&str>)?,
        MenuItem::with_id(app, "recent-1", "", true, None::<&str>)?,
        MenuItem::with_id(app, "recent-2", "", true, None::<&str>)?,
        MenuItem::with_id(app, "recent-3", "", true, None::<&str>)?,
        MenuItem::with_id(app, "recent-4", "", true, None::<&str>)?,
        MenuItem::with_id(app, "recent-5", "", true, None::<&str>)?,
        MenuItem::with_id(app, "recent-6", "", true, None::<&str>)?,
        MenuItem::with_id(app, "recent-7", "", true, None::<&str>)?,
        MenuItem::with_id(app, "recent-8", "", true, None::<&str>)?,
        MenuItem::with_id(app, "recent-9", "", true, None::<&str>)?,
    ];
    let recent_clear = MenuItem::with_id(app, "recent-clear", "Clear Menu", true, None::<&str>)?;
    // Disabled, so it never emits and needs no handler; a submenu that opens onto nothing reads
    // as broken, one that says it is empty does not.
    let recent_empty = MenuItem::new(app, "No Recent Files", false, None::<&str>)?;
    // Built empty and filled by RecentFiles::rebuild once the menu is installed.
    let recent = Submenu::with_items(app, "Open Recent", true, &[])?;

    let file = Submenu::with_items(
        app,
        "File",
        true,
        &[
            &MenuItem::with_id(app, "file-new",  "New Notebook", true, Some("CmdOrCtrl+N"))?,
            &MenuItem::with_id(app, "open",      "Open…",        true, Some("CmdOrCtrl+O"))?,
            &recent,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "save",      "Save",         true, Some("CmdOrCtrl+S"))?,
            &MenuItem::with_id(app, "save-as",   "Save As…",     true, Some("CmdOrCtrl+Shift+S"))?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "file-print", "Print…",      true, Some("CmdOrCtrl+P"))?,
            &PredefinedMenuItem::separator(app)?,
            // Two different closes, deliberately both present: one closes the focused NOTEBOOK,
            // the other the window. Cmd+W goes to the notebook because that is the one a reader
            // reaches for repeatedly.
            &MenuItem::with_id(app, "file-close", "Close Notebook", true, Some("CmdOrCtrl+W"))?,
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
            &MenuItem::with_id(app, "edit-comment", "Un/Comment Selection", true,
                               Some("CmdOrCtrl+/"))?,
            &MenuItem::with_id(app, "edit-indent",  "Indent Selected Lines",  true, None::<&str>)?,
            &MenuItem::with_id(app, "edit-outdent", "Outdent Selected Lines", true, None::<&str>)?,
            &MenuItem::with_id(app, "edit-dupLine", "Duplicate Line", true,
                               Some("CmdOrCtrl+Shift+L"))?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "edit-copyInputAbove", "Copy Input from Above", true,
                               Some("CmdOrCtrl+L"))?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "edit-findDoc", "Documentation for Selection", true,
                               Some("CmdOrCtrl+Shift+F"))?,
        ],
    )?;

    let insert = Submenu::with_items(
        app,
        "Insert",
        true,
        &[
            &MenuItem::with_id(app, "insert-code",    "Input Cell",   true, Some("CmdOrCtrl+B"))?,
            &MenuItem::with_id(app, "insert-text",    "Text Cell",    true, None::<&str>)?,
            &MenuItem::with_id(app, "insert-section", "Section Cell", true, None::<&str>)?,
        ],
    )?;

    // Convert To is flat rather than a nested submenu: a style is a one-gesture choice, and macOS
    // submenus cost a second gesture to reach. The heading levels are grouped after a separator and
    // run in outline order, so the block reads as the document structure it sets.
    let cell = Submenu::with_items(
        app,
        "Cell",
        true,
        &[
            &MenuItem::with_id(app, "cell-toInput",   "Convert to Input",   true, None::<&str>)?,
            &MenuItem::with_id(app, "cell-toText",    "Convert to Text",    true, None::<&str>)?,
            &PredefinedMenuItem::separator(app)?,
            // The heading ladder, in outline order, so the menu reads as the document structure it
            // sets. Still flat rather than a nested submenu: a style is a one-gesture choice, and
            // macOS submenus cost a second gesture to reach.
            &MenuItem::with_id(app, "cell-toTitle",    "Convert to Title",    true, None::<&str>)?,
            &MenuItem::with_id(app, "cell-toSubtitle", "Convert to Subtitle", true, None::<&str>)?,
            &MenuItem::with_id(app, "cell-toChapter",  "Convert to Chapter",  true, None::<&str>)?,
            &MenuItem::with_id(app, "cell-toSection",  "Convert to Section",  true, None::<&str>)?,
            &MenuItem::with_id(app, "cell-toSubsection", "Convert to Subsection", true,
                               None::<&str>)?,
            &MenuItem::with_id(app, "cell-toSubsubsection", "Convert to Subsubsection", true,
                               None::<&str>)?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "cell-divide",    "Divide Cell", true,
                               Some("CmdOrCtrl+Shift+D"))?,
            &MenuItem::with_id(app, "cell-merge",     "Merge Cells", true,
                               Some("CmdOrCtrl+Shift+M"))?,
            &MenuItem::with_id(app, "cell-duplicate", "Duplicate Cell", true, None::<&str>)?,
            &MenuItem::with_id(app, "cell-delete",    "Delete Cell",    true, None::<&str>)?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "cell-clearOutput",    "Delete Output",     true, None::<&str>)?,
            &MenuItem::with_id(app, "cell-clearAllOutput", "Delete All Output", true, None::<&str>)?,
        ],
    )?;

    let evaluation = Submenu::with_items(
        app,
        "Evaluation",
        true,
        &[
            &MenuItem::with_id(app, "eval-cell", "Evaluate Cell", true,
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
            &MenuItem::with_id(app, "gfx-plot",     "Plot Documentation",     true, None::<&str>)?,
            &MenuItem::with_id(app, "gfx-image",    "Image Documentation",    true, None::<&str>)?,
            &MenuItem::with_id(app, "gfx-image3d",  "Image3D Documentation",  true, None::<&str>)?,
            &PredefinedMenuItem::separator(app)?,
            &MenuItem::with_id(app, "gfx-graphics", "Graphics Documentation", true, None::<&str>)?,
        ],
    )?;

    // `toggle-dark` stays a plain item. `toggle-autocomplete` and the scheme
    // items are CheckMenuItems whose checkmark mirrors a webview store, kept in
    // step by sync_view_menu. The ids are a contract with menuCommands.ts and are
    // spelled out longhand (not generated), because check_menu_ids.py reads both
    // sides as source text; the scheme labels mirror COLOR_SCHEMES in schemes.ts.
    let toggle_autocomplete =
        CheckMenuItem::with_id(app, "toggle-autocomplete", "Autocomplete", true, true, None::<&str>)?;
    let sc_default   = CheckMenuItem::with_id(app, "scheme-default",         "Default (adaptive)", true, true,  None::<&str>)?;
    let sc_dracula   = CheckMenuItem::with_id(app, "scheme-dracula",         "Dracula",            true, false, None::<&str>)?;
    let sc_nord      = CheckMenuItem::with_id(app, "scheme-nord",            "Nord",               true, false, None::<&str>)?;
    let sc_monokai   = CheckMenuItem::with_id(app, "scheme-monokai",         "Monokai",            true, false, None::<&str>)?;
    let sc_sol_dark  = CheckMenuItem::with_id(app, "scheme-solarized-dark",  "Solarized Dark",     true, false, None::<&str>)?;
    let sc_sol_light = CheckMenuItem::with_id(app, "scheme-solarized-light", "Solarized Light",    true, false, None::<&str>)?;
    let sc_gruvbox   = CheckMenuItem::with_id(app, "scheme-gruvbox",         "Gruvbox",            true, false, None::<&str>)?;
    let sc_one_dark  = CheckMenuItem::with_id(app, "scheme-one-dark",        "One Dark",           true, false, None::<&str>)?;
    let sc_tokyo     = CheckMenuItem::with_id(app, "scheme-tokyo-night",     "Tokyo Night",        true, false, None::<&str>)?;
    let sc_github    = CheckMenuItem::with_id(app, "scheme-github-light",    "GitHub Light",       true, false, None::<&str>)?;
    let sc_sixties   = CheckMenuItem::with_id(app, "scheme-sixties",         "60s (Teletype)",       true, false, None::<&str>)?;
    let sc_seventies = CheckMenuItem::with_id(app, "scheme-seventies",       "70s (Green Phosphor)", true, false, None::<&str>)?;
    let sc_eighties  = CheckMenuItem::with_id(app, "scheme-eighties",        "80s (Brown/Orange)", true, false, None::<&str>)?;
    let sc_grayscale = CheckMenuItem::with_id(app, "scheme-grayscale",       "Grayscale",          true, false, None::<&str>)?;
    let sc_off       = CheckMenuItem::with_id(app, "scheme-off",             "Off (none)",         true, false, None::<&str>)?;
    let syntax = Submenu::with_items(app, "Syntax Highlighting", true, &[
        &sc_default, &sc_dracula, &sc_nord, &sc_monokai, &sc_sol_dark,
        &sc_sol_light, &sc_gruvbox, &sc_one_dark, &sc_tokyo, &sc_github,
        &sc_sixties, &sc_seventies, &sc_eighties, &sc_grayscale,
        &PredefinedMenuItem::separator(app)?,
        &sc_off,
    ])?;
    let view = Submenu::with_items(
        app,
        "View",
        true,
        &[
            &MenuItem::with_id(app, "toggle-dark", "Toggle Dark Mode", true,
                               Some("CmdOrCtrl+Shift+T"))?,
            &PredefinedMenuItem::separator(app)?,
            &toggle_autocomplete,
            &syntax,
        ],
    )?;
    let view_menu = ViewMenu {
        autocomplete: toggle_autocomplete,
        schemes: vec![
            ("default".into(), sc_default),
            ("dracula".into(), sc_dracula),
            ("nord".into(), sc_nord),
            ("monokai".into(), sc_monokai),
            ("solarized-dark".into(), sc_sol_dark),
            ("solarized-light".into(), sc_sol_light),
            ("gruvbox".into(), sc_gruvbox),
            ("one-dark".into(), sc_one_dark),
            ("tokyo-night".into(), sc_tokyo),
            ("github-light".into(), sc_github),
            ("sixties".into(), sc_sixties),
            ("seventies".into(), sc_seventies),
            ("eighties".into(), sc_eighties),
            ("grayscale".into(), sc_grayscale),
            ("off".into(), sc_off),
        ],
    };

    let menu = Menu::with_items(app, &[
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
    ])?;

    let recent = RecentFiles::new(app, recent, recent_slots, recent_clear, recent_empty)?;
    Ok((menu, recent, view_menu))
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_log::Builder::default().build())
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_fs::init())
        // Save-on-close gate, half one. The red traffic-light button and
        // File > Close Window arrive here as CloseRequested. Hold the close and
        // let the webview's Save/Don't Save/Cancel modal decide; it flips the
        // QuitGuard and calls `confirm_and_quit` (or `cancel_quit`). Cmd+Q and
        // dock > Quit never reach this event — they take the ExitRequested path
        // in the run handler below.
        .on_window_event(|window, event| {
            if let tauri::WindowEvent::CloseRequested { api, .. } = event {
                let guard = window.state::<QuitGuard>();
                if !guard.confirmed.load(Ordering::SeqCst) {
                    api.prevent_close();
                    let _ = window.emit("quit-requested", ());
                }
            }
        })
        .setup(|app| {
            // The gate's shared state, managed before any close/quit event can fire.
            app.manage(QuitGuard::default());

            // Native menu (desktop only — mobile has no app menu bar).
            #[cfg(desktop)]
            {
                let (menu, recent, view_menu) = build_menu(app)?;
                app.set_menu(menu)?;
                // Draw File > Open Recent only once the menu is installed: the submenu has to
                // be attached to a live NSMenu before appending to it reaches the menu bar.
                if let Err(e) = recent.rebuild() {
                    log::warn!("{e}");
                }
                app.manage(recent);
                app.manage(view_menu);
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
            syntax_spans,
            eval_once,
            restart_kernel,
            interrupt_kernel,
            ping_kernel,
            save_notebook,
            load_notebook,
            save_library,
            load_library,
            set_window_title,
            recent_files,
            push_recent_file,
            forget_recent_file,
            clear_recent_files,
            confirm_and_quit,
            cancel_quit,
            sync_view_menu,
        ])
        .build(tauri::generate_context!())
        .expect("error while building tauri application")
        // Save-on-close gate, half two. Cmd+Q and dock > Quit send
        // NSApplication terminate:, which bypasses the window's CloseRequested
        // entirely and arrives as ExitRequested. Same gate, same webview modal.
        .run(|app_handle, event| {
            if let tauri::RunEvent::ExitRequested { api, .. } = event {
                let guard = app_handle.state::<QuitGuard>();
                // The user already confirmed; this is our own app.exit(0). Let it go.
                if guard.confirmed.load(Ordering::SeqCst) {
                    return;
                }
                // A prompt is already pending. A SECOND quit attempt means the
                // user insists (or the webview is wedged): let the exit through
                // rather than ever trapping the app with an unanswerable prompt.
                if guard.prevented.swap(true, Ordering::SeqCst) {
                    return;
                }
                api.prevent_exit();
                let _ = app_handle.emit("quit-requested", ());
            }
        });
}
