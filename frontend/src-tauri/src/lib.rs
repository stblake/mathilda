mod commands;
mod kernel;

use commands::{evaluate_cell, interrupt_kernel, load_notebook, ping_kernel, restart_kernel, save_notebook};
use kernel::MathildaKernel;
use tauri::Manager;

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_log::Builder::default().build())
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_fs::init())
        .setup(|app| {
            // Create and manage the kernel handle synchronously so that
            // State<MathildaKernel> is always available when commands run.
            // The actual subprocess spawn happens in the async task below.
            let kernel = MathildaKernel::empty(app.handle().clone());
            app.manage(kernel);

            let handle = app.handle().clone();
            tauri::async_runtime::spawn(async move {
                let kernel = handle.state::<MathildaKernel>();
                if let Err(e) = kernel.start().await {
                    log::error!("Kernel failed to start: {e}");
                } else {
                    log::info!("Mathilda kernel ready");
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
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
