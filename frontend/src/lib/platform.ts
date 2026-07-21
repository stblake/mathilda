// platform.ts — runtime device-capability detection.
//
// The notebook canvas offers mouse affordances (drag a card to move it,
// drag on empty space to rubber-band select) that make no sense without a
// pointer you can hover and click precisely. On a touch device those same
// gestures fight the natural expectation that one finger pans and two fingers
// zoom. We branch on the device being touch-primary rather than sniffing the
// OS, so it is correct on Android, iOS, and touch tablets alike with no
// Tauri/OS plugin.
//
// Detection: `navigator.maxTouchPoints > 0` is the reliable signal for "has a
// touchscreen" and, unlike the `(pointer: coarse)` media query, it is true in
// the Android emulator (whose WebView advertises a mouse-like fine pointer
// because the host cursor drives it, so `(pointer: coarse)` reports false
// there). We keep the media query as a fallback for any engine that reports
// touch only that way. A conventional desktop (macOS/Windows Tauri webview)
// has maxTouchPoints === 0 and a fine pointer, so it stays on the mouse path.
//
// Evaluated once at module load; the primary input type does not change during
// a session on real devices.
export const isTouchDevice: boolean =
  typeof window !== 'undefined' &&
  ((typeof navigator !== 'undefined' && navigator.maxTouchPoints > 0) ||
    (typeof window.matchMedia === 'function' &&
      window.matchMedia('(pointer: coarse)').matches));
