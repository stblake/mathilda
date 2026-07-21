//! ffi.rs — safe Rust bindings to the in-process Mathilda kernel.
//!
//! On desktop the notebook talks to a spawned `mathilda` sidecar over stdio
//! (see `kernel.rs`). Mobile OS sandboxes (iOS/Android) forbid spawning child
//! processes, so there the kernel is compiled to a static library
//! (`libmathilda.a`, built from `src/ffi/mathilda_ffi.c`) and linked directly
//! into this app. This module wraps that C ABI in safe Rust.
//!
//! Safety: the C kernel uses a process-global symbol table and is NOT
//! reentrant. Every function here must be called from a single logical thread
//! at a time — callers serialize through a mutex (see `kernel_ffi.rs`).

#![allow(dead_code)] // desktop builds compile this module but call none of it.

use std::ffi::{CStr, CString};
use std::os::raw::c_char;

extern "C" {
    fn mathilda_ffi_init();
    fn mathilda_ffi_set_home(dir: *const c_char);
    fn mathilda_ffi_eval(input: *const c_char) -> *mut c_char;
    fn mathilda_ffi_eval_latex(input: *const c_char) -> *mut c_char;
    fn mathilda_ffi_eval_json(input: *const c_char) -> *mut c_char;
    fn mathilda_ffi_free(s: *mut c_char);
    fn mathilda_ffi_version() -> *const c_char;
}

/// Point the kernel at the bundled `internal/` module tree (contains `init.m`
/// and the derivative/integral tables). Call once before `init`. On iOS this
/// is the app bundle's resource directory.
pub fn set_home(dir: &str) {
    if let Ok(c) = CString::new(dir) {
        // SAFETY: `c` outlives the call; the C side copies the string.
        unsafe { mathilda_ffi_set_home(c.as_ptr()) };
    }
}

/// Initialize the kernel (idempotent). Loads builtins + `init.m`.
pub fn init() {
    // SAFETY: no arguments; C side guards against double-init.
    unsafe { mathilda_ffi_init() };
}

/// Take ownership of a C string returned by the kernel, copying it into a Rust
/// `String` and freeing the C allocation. Returns "" on NULL / bad UTF-8.
unsafe fn take_c_string(ptr: *mut c_char) -> String {
    if ptr.is_null() {
        return String::new();
    }
    let s = CStr::from_ptr(ptr).to_string_lossy().into_owned();
    mathilda_ffi_free(ptr);
    s
}

/// Evaluate `input`, returning the formatted output text (same as the REPL's
/// `Out[n]=` line). Parse errors return "$Failed (parse error)".
pub fn eval(input: &str) -> String {
    let c = match CString::new(input) {
        Ok(c) => c,
        Err(_) => return "$Failed (embedded NUL in input)".into(),
    };
    // SAFETY: `c` outlives the call; the returned pointer is owned by us and
    // freed inside `take_c_string`.
    unsafe { take_c_string(mathilda_ffi_eval(c.as_ptr())) }
}

/// Evaluate `input`, returning a LaTeX rendering of the result (for KaTeX).
pub fn eval_latex(input: &str) -> String {
    let c = match CString::new(input) {
        Ok(c) => c,
        Err(_) => return String::new(),
    };
    // SAFETY: as in `eval`.
    unsafe { take_c_string(mathilda_ffi_eval_latex(c.as_ptr())) }
}

/// Evaluate `input` exactly ONCE and return a JSON object string describing the
/// result — the notebook's preferred entry point. It carries the Plotly plot
/// payload for `Graphics`/`Graphics3D` (which `eval`/`eval_latex` cannot) and,
/// by evaluating a single time, avoids double-running side effects. Shape:
///   {"type":"plot","payload":{...}}
///   {"type":"expr","payload":"...","latex":"..."}   (latex may be absent)
///   {"type":"error","message":"..."}
/// See `mathilda_ffi_eval_json` in `src/ffi/mathilda_ffi.c`.
pub fn eval_json(input: &str) -> String {
    let c = match CString::new(input) {
        Ok(c) => c,
        Err(_) => return r#"{"type":"error","message":"embedded NUL in input"}"#.into(),
    };
    // SAFETY: as in `eval`.
    unsafe { take_c_string(mathilda_ffi_eval_json(c.as_ptr())) }
}

/// Kernel version string (e.g. "0.015").
pub fn version() -> String {
    // SAFETY: the C side returns a pointer to static storage; we only borrow it.
    unsafe {
        let p = mathilda_ffi_version();
        if p.is_null() {
            String::new()
        } else {
            CStr::from_ptr(p).to_string_lossy().into_owned()
        }
    }
}
