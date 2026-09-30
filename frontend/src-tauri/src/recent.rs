//! recent.rs — the File > Open Recent list and the native submenu that shows it.
//!
//! Split the way `notebook_format.rs` is: the list logic is pure functions of plain data,
//! unit tested without an app handle, and the Tauri glue (the submenu handles, the managed
//! state, the JSON store) sits on top of it.
//!
//! WHO OWNS WHAT. The submenu is native, so it can only be built and mutated here in Rust;
//! and the list has to be on disk before the first paint, which rules out the webview
//! pushing it up after mount. But the webview owns every open and save path there is, so it
//! is what RECORDS an entry (`push_recent_file`) and what ACTS on a click — Rust emits
//! `menu:recent-3` like any other menu item and the front end maps the index back to a path.
//!
//! THE ID POOL IS FIXED AND LITERAL. `tools/check_menu_ids.py` reads the ids out of
//! `lib.rs` as source text and diffs them against the handler cases in `menuCommands.ts`,
//! so an id built at runtime with `format!("recent-{i}")` would be invisible to the Rust
//! half of that check and the TypeScript half would be reported as dead wiring. The ten
//! slots are therefore spelled out longhand in `lib.rs` and handed to [`RecentFiles::new`]
//! in order; `MAX_RECENT` is derived from how many arrived, so the number 10 lives in
//! exactly one place.
//!
//! What this is NOT: a general preferences store. `src/lib/properties.ts` records the
//! decision that UI preferences do not persist until they all can. A recent-files list is
//! document history rather than a preference — an Open Recent that forgets on restart is
//! not the feature — so it persists, and the rest of that decision stands.

use std::path::{Path, PathBuf};
use std::sync::Mutex;

use tauri::menu::{MenuItem, PredefinedMenuItem, Submenu};
use tauri::{Manager, Wry};

/// Filename of the JSON store inside the app config directory.
const STORE_FILE: &str = "recent.json";

// =============================================================================
// Pure list logic — no Tauri, no filesystem. Everything here is unit tested.
// =============================================================================

/// Put `path` at the front of `list`, removing any earlier occurrence, and cap the result
/// at `max`.
///
/// Move-to-front rather than push: reopening the file you opened an hour ago should lift it
/// back to the top, not add a second line saying the same thing. Comparison is on the path
/// as given — two spellings of one file (a symlink, a trailing slash) would read as two
/// entries, which is the same thing every other app's Open Recent does.
pub fn promote(list: &mut Vec<PathBuf>, path: &Path, max: usize) {
    list.retain(|p| p != path);
    list.insert(0, path.to_path_buf());
    list.truncate(max);
}

/// Drop `path` from `list` if it is there. Used when an entry turns out not to open.
pub fn forget(list: &mut Vec<PathBuf>, path: &Path) {
    list.retain(|p| p != path);
}

/// The menu label for each path, in the same order.
///
/// A file name on its own, as macOS shows it — except where two entries share one, in which
/// case both get their parent folder appended (`notes.mnb` / `notes.mnb — drafts`). Two menu
/// lines reading exactly the same thing is the failure worth spending eight lines to avoid:
/// the reader cannot tell which is which, and picking wrong replaces their canvas.
pub fn labels(list: &[PathBuf]) -> Vec<String> {
    let names: Vec<String> = list
        .iter()
        .map(|p| {
            p.file_name()
                .map(|n| n.to_string_lossy().into_owned())
                // A path ending in `..` or `/` has no file name; showing the whole path is
                // better than showing nothing.
                .unwrap_or_else(|| p.to_string_lossy().into_owned())
        })
        .collect();

    names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let ambiguous = names
                .iter()
                .enumerate()
                .any(|(j, other)| j != i && other == name);
            if !ambiguous {
                return name.clone();
            }
            match list[i].parent().and_then(|d| d.file_name()) {
                Some(dir) => format!("{name} — {}", dir.to_string_lossy()),
                None => name.clone(),
            }
        })
        .collect()
}

/// Parse the store's JSON (an array of path strings), keeping only entries that still exist.
///
/// A recent list is a list of files, not of names: an entry whose file has been deleted or
/// moved cannot be opened, so offering it is offering a dead end. Pruning on load is the
/// cheap half; `App.svelte` prunes the other half, when an open actually fails.
pub fn parse_store(json: &str, max: usize) -> Vec<PathBuf> {
    let parsed: Vec<String> = serde_json::from_str(json).unwrap_or_default();
    let mut out: Vec<PathBuf> = Vec::new();
    for s in parsed {
        let p = PathBuf::from(s);
        if p.is_file() && !out.contains(&p) {
            out.push(p);
        }
    }
    out.truncate(max);
    out
}

/// Serialise the list for the store.
pub fn render_store(list: &[PathBuf]) -> String {
    let strs: Vec<String> = list.iter().map(|p| p.to_string_lossy().into_owned()).collect();
    serde_json::to_string_pretty(&strs).unwrap_or_else(|_| "[]".to_string())
}

// =============================================================================
// Managed state: the list, its store, and the native submenu that mirrors it.
// =============================================================================

/// Everything File > Open Recent needs, held in Tauri managed state.
///
/// The menu handles are all `Arc`-backed and `Send + Sync`, so `app.manage` takes the struct
/// whole and every command can reach the same submenu.
pub struct RecentFiles {
    paths: Mutex<Vec<PathBuf>>,
    /// The "Open Recent" submenu itself, nested inside File.
    menu: Submenu<Wry>,
    /// `recent-0` … `recent-9`, in order. Slot *i* shows list entry *i*.
    slots: Vec<MenuItem<Wry>>,
    /// `recent-clear` — the "Clear Menu" item at the foot, macOS style.
    clear: MenuItem<Wry>,
    /// Shown, disabled, when the list is empty. A submenu that opens onto nothing looks
    /// broken; one that says so does not.
    placeholder: MenuItem<Wry>,
    /// The rule above the "Clear Menu" item. Held rather than made fresh each redraw, so
    /// every handle the submenu ever contains outlives the menu — see `rebuild`.
    separator: PredefinedMenuItem<Wry>,
    store: PathBuf,
}

impl RecentFiles {
    /// Adopt the menu handles `lib.rs` just built, then load the stored list.
    ///
    /// Does not touch the menu — call [`RecentFiles::rebuild`] once the menu is installed.
    pub fn new(
        app: &tauri::App,
        menu: Submenu<Wry>,
        slots: Vec<MenuItem<Wry>>,
        clear: MenuItem<Wry>,
        placeholder: MenuItem<Wry>,
    ) -> tauri::Result<Self> {
        let store = app
            .path()
            .app_config_dir()
            .unwrap_or_else(|_| PathBuf::from("."))
            .join(STORE_FILE);

        let max = slots.len();
        let paths = std::fs::read_to_string(&store)
            .map(|j| parse_store(&j, max))
            .unwrap_or_default();

        Ok(Self {
            paths: Mutex::new(paths),
            separator: PredefinedMenuItem::separator(app)?,
            menu,
            slots,
            clear,
            placeholder,
            store,
        })
    }

    /// How many entries the list holds — the size of the id pool `lib.rs` built, so the
    /// number is not written down twice.
    fn max(&self) -> usize {
        self.slots.len()
    }

    /// The list, newest first, as strings for the webview.
    pub fn list(&self) -> Vec<String> {
        self.paths
            .lock()
            .map(|p| p.iter().map(|p| p.to_string_lossy().into_owned()).collect())
            .unwrap_or_default()
    }

    /// Record a file that was just opened or saved.
    pub fn push(&self, path: &Path) -> Result<(), String> {
        self.mutate(|list, max| promote(list, path, max))
    }

    /// Drop a file — used when opening it failed, so the menu stops offering a dead end.
    pub fn forget(&self, path: &Path) -> Result<(), String> {
        self.mutate(|list, _| forget(list, path))
    }

    /// Empty the list (the "Clear Menu" item).
    pub fn clear(&self) -> Result<(), String> {
        self.mutate(|list, _| list.clear())
    }

    /// Apply `f` to the list, then persist and redraw. One funnel, so no caller can change
    /// the list and forget to do either.
    fn mutate<F: FnOnce(&mut Vec<PathBuf>, usize)>(&self, f: F) -> Result<(), String> {
        let max = self.max();
        {
            let mut list = self.paths.lock().map_err(|e| format!("recent: {e}"))?;
            f(&mut list, max);
            let json = render_store(&list);
            if let Some(dir) = self.store.parent() {
                let _ = std::fs::create_dir_all(dir);
            }
            // A store that cannot be written costs the list its memory across restarts and
            // nothing else, so it is reported and stepped over rather than failing the open
            // or save that triggered it.
            if let Err(e) = std::fs::write(&self.store, json) {
                log::warn!("recent: could not write {}: {e}", self.store.display());
            }
        }
        self.rebuild()
    }

    /// Redraw the submenu from the list.
    ///
    /// Drain and re-append rather than rebuild the menu: `App::set_menu` tears the whole menu
    /// bar down first (`remove_menu` → `setMainMenu(None)` on macOS), so redrawing ten rows
    /// that way would flicker the Mathilda and Edit menus too, and drop every key equivalent
    /// for an instant. muda supports the drain directly — an `NSMenuItem` is minted fresh from
    /// the Rust-side item on each append and merely detached on remove, so a pool handle is
    /// reusable and `set_text` on a detached item is picked up by its next append.
    ///
    /// Three invariants, each load-bearing:
    ///
    /// 1. **The whole redraw happens in ONE main-thread hop.** Every muda call individually
    ///    posts to the main thread and blocks the caller until the event loop pumps, so from
    ///    an async command a naive redraw is ~thirty blocking round trips — and each stalls
    ///    for as long as a menu is being tracked or a modal file dialog is up. Inside
    ///    `run_on_main_thread` the inner calls take muda's already-on-main-thread branch and
    ///    run inline.
    /// 2. **A pool handle is appended at most once per redraw.** muda removes a child by id
    ///    but detaches *every* `NSMenuItem` it has in that menu, so the same handle appended
    ///    twice leaves its item list and the real NSMenu disagreeing, and later `remove_at`
    ///    calls then target the wrong row.
    /// 3. **The lock is released before any menu call.** The main thread takes this lock too
    ///    (a click arrives there), so holding the guard across a call that waits on the main
    ///    thread is a two-party deadlock. The list is snapshotted into labels first.
    pub fn rebuild(&self) -> Result<(), String> {
        let labels = {
            let list = self.paths.lock().map_err(|e| format!("recent: {e}"))?;
            labels(&list)
        };

        // Cloned handles, not `&self`: the closure outlives this call when it is posted from
        // a worker thread. Each is an Arc onto the same native item.
        let menu = self.menu.clone();
        let slots = self.slots.clone();
        let clear = self.clear.clone();
        let placeholder = self.placeholder.clone();
        let separator = self.separator.clone();

        self.menu
            .app_handle()
            .run_on_main_thread(move || {
                let redraw = || -> tauri::Result<()> {
                    while menu.remove_at(0)?.is_some() {}
                    if labels.is_empty() {
                        return menu.append(&placeholder);
                    }
                    for (slot, label) in slots.iter().zip(labels.iter()) {
                        slot.set_text(label)?;
                        menu.append(slot)?;
                    }
                    menu.append(&separator)?;
                    menu.append(&clear)
                };
                // The redraw may run after this call returned, so a failure has nowhere to be
                // returned to. A stale Open Recent submenu is not worth failing the open or
                // save that triggered it.
                if let Err(e) = redraw() {
                    log::warn!("recent menu: could not redraw: {e}");
                }
            })
            .map_err(|e| format!("recent menu: {e}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(s: &str) -> PathBuf {
        PathBuf::from(s)
    }

    fn strs(list: &[PathBuf]) -> Vec<String> {
        list.iter().map(|p| p.to_string_lossy().into_owned()).collect()
    }

    #[test]
    fn promote_puts_the_newest_first() {
        let mut list = vec![p("/a.mnb"), p("/b.mnb")];
        promote(&mut list, &p("/c.mnb"), 10);
        assert_eq!(strs(&list), ["/c.mnb", "/a.mnb", "/b.mnb"]);
    }

    #[test]
    fn reopening_moves_to_front_rather_than_duplicating() {
        let mut list = vec![p("/a.mnb"), p("/b.mnb"), p("/c.mnb")];
        promote(&mut list, &p("/c.mnb"), 10);
        assert_eq!(strs(&list), ["/c.mnb", "/a.mnb", "/b.mnb"]);
    }

    #[test]
    fn the_list_is_capped_and_drops_the_oldest() {
        let mut list: Vec<PathBuf> = (0..3).map(|i| p(&format!("/{i}.mnb"))).collect();
        promote(&mut list, &p("/new.mnb"), 3);
        assert_eq!(strs(&list), ["/new.mnb", "/0.mnb", "/1.mnb"]);
    }

    #[test]
    fn forget_removes_one_entry_and_tolerates_a_miss() {
        let mut list = vec![p("/a.mnb"), p("/b.mnb")];
        forget(&mut list, &p("/a.mnb"));
        assert_eq!(strs(&list), ["/b.mnb"]);
        forget(&mut list, &p("/nope.mnb"));
        assert_eq!(strs(&list), ["/b.mnb"]);
    }

    #[test]
    fn labels_are_file_names() {
        let list = vec![p("/home/u/work/notes.mnb"), p("/home/u/lib.lb")];
        assert_eq!(labels(&list), ["notes.mnb", "lib.lb"]);
    }

    #[test]
    fn a_shared_file_name_gains_its_folder_on_both_entries() {
        // The point of the disambiguation: neither line may read just "notes.mnb", or the
        // reader cannot tell which one they are about to open.
        let list = vec![
            p("/home/u/drafts/notes.mnb"),
            p("/home/u/final/notes.mnb"),
            p("/home/u/other.mnb"),
        ];
        assert_eq!(
            labels(&list),
            ["notes.mnb — drafts", "notes.mnb — final", "other.mnb"]
        );
    }

    #[test]
    fn parse_store_drops_files_that_are_gone() {
        // This source file exists; the other two do not.
        let me = file!();
        let json = format!(r#"["{me}", "/definitely/not/here.mnb", "{me}"]"#);
        let got = parse_store(&json, 10);
        assert_eq!(strs(&got), [me], "only the existing, de-duplicated path survives");
    }

    #[test]
    fn parse_store_survives_a_corrupt_file() {
        assert!(parse_store("not json at all", 10).is_empty());
        assert!(parse_store("", 10).is_empty());
    }

    #[test]
    fn the_store_round_trips() {
        let me = PathBuf::from(file!());
        let json = render_store(&[me.clone()]);
        assert_eq!(parse_store(&json, 10), vec![me]);
    }
}
