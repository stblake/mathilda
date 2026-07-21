// platform.ts — runtime device-capability detection.
//
// The notebook canvas offers mouse affordances (drag a card to move it, drag on
// empty space to rubber-band select) that make no sense on a touch device,
// where one finger should pan and two should zoom. `isTouchDevice` selects
// between the two interaction models.
//
// We key off the OS/WebView, NOT a pointer-capability media query or
// `maxTouchPoints` alone, because neither distinguishes the two cases that must
// behave differently:
//   - Android emulator — advertises a mouse (`(pointer: coarse)` is *false*,
//     `(hover: hover)` is true) yet is a touch device. `maxTouchPoints` catches
//     it; the media query does not.
//   - Desktop touchscreen (a Windows/Linux 2-in-1 running the *desktop* app) —
//     has touch points but is mouse-primary here and must stay on the mouse
//     path. A bare `maxTouchPoints > 0` check would wrongly flip it to touch.
//
// The Tauri WebView user-agent is the reliable discriminator: Android and iOS
// WebViews carry their platform token; desktop WebViews (macOS/Windows/Linux)
// do not. iPadOS is the one exception — its WKWebView reports a desktop
// "Macintosh" UA — so it is caught via touch points (a real Mac has none).
function detectTouchDevice(): boolean {
  if (typeof navigator === 'undefined') return false;
  const ua = navigator.userAgent || '';
  if (/Android|iPhone|iPod|iPad/i.test(ua)) return true;
  // iPadOS in a native WKWebView presents a "Macintosh" UA; desktop Macs have
  // no touchscreen, so touch points > 1 means this is really an iPad.
  if (/Macintosh/.test(ua) && navigator.maxTouchPoints > 1) return true;
  return false;
}

// Evaluated once at module load; the platform does not change during a session.
export const isTouchDevice: boolean =
  typeof window !== 'undefined' && detectTouchDevice();
