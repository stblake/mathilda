// searchHighlight.ts — show the find bar's current match WITHOUT moving focus.
//
// The find bar used to "reveal" a match by selecting it in the owning editor and
// focusing that editor. Selection is the only highlight a CodeMirror view without
// drawSelection has, and a DOM selection only shows where focus is -- so focus
// had to move, and then the next Enter went to the CELL and replaced the match
// with a newline. The find bar must keep focus for Enter / Shift+Enter to walk
// the matches, so the current match is painted instead of selected:
//
//   * code cells: a CodeMirror mark decoration held in a StateField and set by a
//     StateEffect. It lives in the editor's own state, so it survives redraws
//     and maps through edits.
//   * prose cells: the CSS Custom Highlight API over a DOM Range. Where the
//     WebView lacks it the cell is still scrolled into view; only the paint is
//     missing.
//
// The selection is placed for real when the bar closes (Escape), which is what
// every editor's find does.

import { StateEffect, StateField } from '@codemirror/state';
import { Decoration, EditorView } from '@codemirror/view';
import type { DecorationSet } from '@codemirror/view';

/** A half-open source range, or null to clear. */
export type MarkRange = { start: number; end: number } | null;

/** Set (or clear, with null) the current-match mark in one editor. */
export const setSearchMark = StateEffect.define<MarkRange>({
  map: (v, mapping) => (v ? { start: mapping.mapPos(v.start), end: mapping.mapPos(v.end) } : null),
});

const currentMatch = Decoration.mark({ class: 'cm-search-current' });

const searchMarkField = StateField.define<DecorationSet>({
  create: () => Decoration.none,
  update(deco, tr) {
    deco = deco.map(tr.changes);
    for (const e of tr.effects) {
      if (!e.is(setSearchMark)) continue;
      const len = tr.state.doc.length;
      const r = e.value;
      const from = r ? Math.max(0, Math.min(r.start, len)) : 0;
      const to = r ? Math.max(0, Math.min(r.end, len)) : 0;
      deco = r && to > from ? Decoration.set([currentMatch.range(from, to)]) : Decoration.none;
    }
    return deco;
  },
  provide: f => EditorView.decorations.from(f),
});

/** The extension a code cell installs so the find bar can mark a match in it. */
export const searchMarkExtension = [
  searchMarkField,
  EditorView.baseTheme({
    '.cm-search-current': {
      background: 'rgba(250, 179, 135, 0.45)',
      outline: '1px solid rgba(250, 179, 135, 0.95)',
      borderRadius: '2px',
    },
  }),
];

/** Mark a range in a code cell's editor and scroll it into view, without focus. */
export function markInEditor(view: EditorView, range: MarkRange) {
  const effects: StateEffect<unknown>[] = [setSearchMark.of(range)];
  if (range) effects.push(EditorView.scrollIntoView(Math.min(range.start, view.state.doc.length), { y: 'nearest' }));
  view.dispatch({ effects });
}

/* ---- prose cells --------------------------------------------------------- */

const HIGHLIGHT_NAME = 'mathilda-search';

/** Map source offsets onto the DOM of an element painted with `innerText =
 *  source`. Setting innerText turns each newline into a <br>, so the element is
 *  a run of text nodes and <br>s; a <br> counts as the one '\n' it replaced.
 *  Returns null when the offsets do not land in the DOM (it is not what this
 *  assumes, and marking the wrong characters is worse than marking none). */
export function rangeForOffsets(el: HTMLElement, start: number, end: number): Range | null {
  const doc = el.ownerDocument;
  const walker = doc.createTreeWalker(el, NodeFilter.SHOW_TEXT | NodeFilter.SHOW_ELEMENT);
  let pos = 0;
  let startNode: Node | null = null, startOff = 0;
  let endNode: Node | null = null, endOff = 0;
  for (let n = walker.nextNode(); n; n = walker.nextNode()) {
    if (n.nodeType === Node.TEXT_NODE) {
      const len = n.textContent?.length ?? 0;
      if (!startNode && start <= pos + len) { startNode = n; startOff = start - pos; }
      if (!endNode && end <= pos + len) { endNode = n; endOff = end - pos; break; }
      pos += len;
    } else if ((n as Element).tagName === 'BR') {
      pos += 1;
    }
  }
  if (!startNode || !endNode) return null;
  const range = doc.createRange();
  range.setStart(startNode, startOff);
  range.setEnd(endNode, endOff);
  return range;
}

/** Paint `range` as the current match (null clears). False when the WebView has
 *  no CSS Custom Highlight API, in which case nothing is painted. */
export function paintDomRange(range: Range | null): boolean {
  const g = globalThis as unknown as {
    CSS?: { highlights?: { set(k: string, v: unknown): void; delete(k: string): void } };
    Highlight?: new (...r: Range[]) => unknown;
  };
  const reg = g.CSS?.highlights;
  if (!reg || !g.Highlight) return false;
  if (range) reg.set(HIGHLIGHT_NAME, new g.Highlight(range));
  else reg.delete(HIGHLIGHT_NAME);
  return true;
}
