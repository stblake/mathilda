/* structuralSelection.ts — a self-contained CodeMirror extension that adds
 * Mathilda's bottom-up (precedence-aware) click selection to ANY EditorView,
 * editable or read-only. Used by the input cells (CellShell) and the read-only
 * output views (CodeView) so both behave identically.
 *
 * The precedence-aware spans come from the kernel (fetchSpans); a synchronous
 * local bracket/token span set is the instant fallback for the first click until
 * the kernel answers, and the gesture recomputes its ladder each click so it
 * upgrades seamlessly. The pure gesture logic lives in structSelect.ts. */

import { EditorView, ViewPlugin, keymap } from '@codemirror/view';
import type { PluginValue, ViewUpdate } from '@codemirror/view';
import { EditorSelection, Prec } from '@codemirror/state';
import { tokenize } from './mathildaLex';
import { buildSpans, pickAtom, clickGesture, emptyGesture, stepGesture } from './structSelect';
import type { Span, GestureState } from './structSelect';
import { fetchSpans } from './ipc';

class StructSel implements PluginValue {
  private localSpans: Span[] | null = null;   // synchronous fallback
  private kernelSpans: Span[] | null = null;  // precedence-aware, from the kernel
  private spansDoc = '';
  private fetching = false;
  private gesture: GestureState = emptyGesture();
  private timer: ReturnType<typeof setTimeout> | null = null;

  constructor(private view: EditorView) {}

  private docOf(): string { return this.view.state.doc.toString(); }

  /** Best spans available right now (kernel if ready for this doc, else local),
   *  kicking off a kernel fetch when needed. */
  private bestSpans(): Span[] {
    const doc = this.docOf();
    if (doc !== this.spansDoc) { this.localSpans = null; this.kernelSpans = null; this.spansDoc = doc; }
    if (this.kernelSpans) return this.kernelSpans;
    if (!this.localSpans) this.localSpans = buildSpans(tokenize(doc), doc);
    void this.request(doc);
    return this.localSpans;
  }

  private async request(doc: string): Promise<void> {
    if (this.fetching || (this.kernelSpans && this.spansDoc === doc)) return;
    this.fetching = true;
    try {
      const raw = await fetchSpans(doc);
      if (this.docOf() === doc) { this.kernelSpans = raw.map(([from, to]) => ({ from, to })); this.spansDoc = doc; }
    } catch {
      /* kernel busy/down: keep the local fallback */
    } finally {
      this.fetching = false;
    }
  }

  private prefetch(): void {
    if (this.timer) clearTimeout(this.timer);
    this.timer = setTimeout(() => void this.request(this.docOf()), 200);
  }

  private select(sp: Span): void {
    this.view.dispatch({ selection: EditorSelection.range(sp.from, sp.to) });
  }

  /** First click selects the token; each further click on it climbs one level. */
  onMousedown(event: MouseEvent): boolean {
    try {
      if (event.button !== 0) return false;
      const pos = this.view.posAtCoords({ x: event.clientX, y: event.clientY });
      if (pos == null) { this.gesture = emptyGesture(); return false; }
      const atom = pickAtom(tokenize(this.docOf()), pos);
      const r = clickGesture(this.gesture, this.bestSpans(), atom);
      this.gesture = r.state;
      if (r.selection) {
        event.preventDefault();
        this.view.focus();   // so the selection renders (esp. in read-only output views)
        this.select(r.selection);
        return true;
      }
      return false;
    } catch {
      this.gesture = emptyGesture();
      return false;
    }
  }

  /** Alt-Up (grow) / Alt-Down (shrink). */
  step(grow: boolean): boolean {
    try {
      const main = this.view.state.selection.main;
      const atom = pickAtom(tokenize(this.docOf()), main.head);
      if (!atom) return false;
      const r = stepGesture(this.gesture, this.bestSpans(), atom, { from: main.from, to: main.to }, grow);
      if (!r.selection) return false;
      this.gesture = r.state;
      this.select(r.selection);
      return true;
    } catch { return false; }
  }

  update(u: ViewUpdate): void {
    if (u.docChanged) {
      this.localSpans = null; this.kernelSpans = null; this.spansDoc = '';
      this.gesture = emptyGesture();
      this.prefetch();
    }
    if (u.focusChanged && u.view.hasFocus) void this.request(this.docOf());
  }

  destroy(): void { if (this.timer) clearTimeout(this.timer); }
}

const structSelPlugin = ViewPlugin.fromClass(StructSel);

/** Drop-in extension: Mathilda bottom-up click selection + Alt-Up/Down. */
export const mathildaStructuralSelection = [
  structSelPlugin,
  // Prec.highest so our mousedown runs before CodeMirror's own selection logic.
  Prec.highest(EditorView.domEventHandlers({
    mousedown(event, view) { return view.plugin(structSelPlugin)?.onMousedown(event) ?? false; },
  })),
  keymap.of([
    { key: 'Alt-ArrowUp',   run: (view) => view.plugin(structSelPlugin)?.step(true) ?? false },
    { key: 'Alt-ArrowDown', run: (view) => view.plugin(structSelPlugin)?.step(false) ?? false },
  ]),
];
