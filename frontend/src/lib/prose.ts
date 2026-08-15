// prose.ts — rendering and editing for Markdown text cells.
//
// The cell-type picker has described a text cell as "Prose / markdown" since it
// was written, and nothing rendered any markdown: a text cell was a
// contenteditable holding raw characters. This module is the half that makes the
// label true, plus the selection edits the toolbar's Text group needs.
//
// Kept out of CellShell.svelte on purpose. The rendering is a pure function of a
// string, and the selection edits are pure functions of (element, markers), so
// both are readable and reviewable without mounting a component -- the same
// reason cellCommands.ts exists next door.

import { marked } from 'marked';

/** Markdown to HTML for a text cell.
 *
 *  `breaks: true` because a text cell is typed like prose, not like a Markdown
 *  file: a single newline is meant as a line break, and requiring two trailing
 *  spaces to get one would read as the cell ignoring the Return key.
 *
 *  ON TRUST. The result is inserted with `{@html}`, so raw HTML in a cell's
 *  source reaches the DOM. That is not a new capability being granted: a
 *  notebook is an executable document whose code cells already evaluate
 *  arbitrary Mathilda, so anyone in a position to put a `<script>` in the prose
 *  of a notebook you open is already in a position to run code in it. What would
 *  be wrong is claiming otherwise -- a regex pass over the output would look like
 *  sanitisation while missing most of what it appears to cover, and `marked`
 *  dropped its own `sanitize` option precisely because it could not do the job.
 *  RefPage.svelte renders generated documentation through the same `{@html}`. */
export function renderProse(src: string): string {
  const text = src ?? '';
  if (!text.trim()) return '';
  return marked.parse(text, { async: false, breaks: true, gfm: true }) as string;
}

/** Where to leave the caret, counted back from the end of the inserted text. */
export type CaretBack = {
  /** Applied when the selection was empty. Usually `after.length`, so typing
   *  continues INSIDE the markers just inserted. */
  collapsed: number;
  /** Applied when text was selected. Usually 0 -- the caret belongs after the
   *  closing marker, since the wrapped words are finished. */
  selected: number;
};

/** Wrap the current selection in `el` with Markdown markers.
 *
 *  Uses `document.execCommand('insertText')`, which is deprecated and is still
 *  the right call here: it is the only API that edits a contenteditable while
 *  keeping the browser's native undo stack intact, and it fires `input` so the
 *  cell's existing handler writes the new source to the store with no extra
 *  plumbing. Replacing a Range by hand gives an edit that Cmd-Z cannot undo,
 *  which is a worse outcome than using a deprecated call in a WebView whose
 *  engine is known.
 *
 *  Returns false and does nothing when there is no selection inside `el`, so a
 *  toolbar button pressed with the caret elsewhere is inert rather than editing
 *  whichever cell happens to be active. */
export function wrapSelection(el: HTMLElement | null, before: string, after: string,
                              caretBack: CaretBack): boolean {
  if (!el) return false;
  const sel = window.getSelection();
  if (!sel || sel.rangeCount === 0) return false;
  if (!sel.anchorNode || !el.contains(sel.anchorNode)) return false;

  const selected = sel.toString();
  if (!document.execCommand('insertText', false, before + selected + after)) return false;

  const back = selected.length === 0 ? caretBack.collapsed : caretBack.selected;
  if (back > 0) moveCaretBack(back);
  return true;
}

/** Step the caret back `n` characters within its own text node.
 *
 *  Deliberately does not walk into previous nodes: every caller moves back by at
 *  most the few characters it just inserted, which are in the node the caret
 *  landed in. A short node means the insert did something unexpected, and leaving
 *  the caret where the browser put it is better than guessing across a boundary. */
function moveCaretBack(n: number): void {
  const sel = window.getSelection();
  if (!sel || sel.rangeCount === 0) return;
  const range = sel.getRangeAt(0);
  const off = range.startOffset - n;
  if (off < 0) return;
  range.setStart(range.startContainer, off);
  range.collapse(true);
  sel.removeAllRanges();
  sel.addRange(range);
}

/* The four the Text group offers. Emphasis uses `**`/`*` rather than `__`/`_`,
 * which read as literal underscores inside identifiers -- and identifiers are
 * exactly what gets written about in a notebook's prose.
 *
 * There is no bullet-list button. Inserting "- " at the start of the current LINE
 * means finding that line's start in a contenteditable whose line boxes may be
 * text nodes, <div>s or <br>s depending on how the content was entered, and a
 * bullet that lands mid-word is worse than a button that is not there. It wants a
 * real caret to verify against, not a guess. */
export const PROSE_BOLD   = { before: '**', after: '**', caret: { collapsed: 2, selected: 0 } };
export const PROSE_ITALIC = { before: '*',  after: '*',  caret: { collapsed: 1, selected: 0 } };
export const PROSE_CODE   = { before: '`',  after: '`',  caret: { collapsed: 1, selected: 0 } };
/** Caret lands between the parentheses either way: with words selected the text
 *  is done and the URL is what you still have to type. */
export const PROSE_LINK   = { before: '[',  after: ']()', caret: { collapsed: 3, selected: 1 } };
