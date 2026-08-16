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
import katex from 'katex';

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
  /* Math comes OUT before Markdown runs and goes back in after. The order is the
     whole point: `$a_1 * b_2$` handed to Markdown first loses `_1 * b_2` to
     emphasis, and no amount of care in the TeX pass afterwards can recover what
     the emphasis rules already ate. */
  const { text: masked, spans, token } = extractMath(text);
  let html = marked.parse(masked, { async: false, breaks: true, gfm: true }) as string;
  for (let i = 0; i < spans.length; i++) {
    html = html.replace(token(i), renderTex(spans[i]));
  }
  return html;
}

/** One extracted math span. */
export type MathSpan = { display: boolean; tex: string };

function renderTex(span: MathSpan): string {
  try {
    /* throwOnError: false makes KaTeX render its own error marker in place rather
       than take the whole cell down for one bad macro -- the same call
       Output.svelte makes for kernel output. */
    return katex.renderToString(span.tex, { throwOnError: false, displayMode: span.display });
  } catch {
    /* KaTeX can still throw on input it cannot even parse into an error. The raw
       source is then the most useful thing to show: it is what was typed. */
    const delim = span.display ? '$$' : '$';
    return `<code>${delim}${escapeHtml(span.tex)}${delim}</code>`;
  }
}

function escapeHtml(s: string): string {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

/** Pull `$…$` and `$$…$$` out of `src`, replacing each with an opaque token.
 *
 *  A hand-written scan rather than a regex, because the things that must NOT be
 *  math are all context: a `$` inside a code span or fence is a literal dollar,
 *  and `\$` is a literal dollar anywhere. A regex over the whole string cannot see
 *  which of those it is inside.
 *
 *  TWO HEURISTICS, both standard (pandoc and markdown-it use the same pair), and
 *  both there to keep prose about money from turning into equations: an opening
 *  `$` must be followed by a non-space, and a closing `$` must be preceded by one.
 *  So "it cost $5 and $6" stays prose, while `$x + 1$` is math. An unclosed `$`
 *  is likewise left exactly as typed rather than swallowing the rest of the cell.
 *
 *  Exported for the checks in scripts/check-prose.mjs, which is the only reason it
 *  is not private -- the tokenizer's edge cases are the part worth testing. */
export function extractMath(src: string): {
  text: string; spans: MathSpan[]; token: (i: number) => string;
} {
  /* A sentinel the source cannot already contain, so a cell that happens to
     discuss this very placeholder still renders correctly. */
  let stem = '@@MATHILDA-MATH-';
  while (src.includes(stem)) stem += 'X';
  const token = (i: number) => `${stem}${i}@@`;

  const spans: MathSpan[] = [];
  let out = '';
  let i = 0;
  const n = src.length;

  while (i < n) {
    const c = src[i];

    /* An escaped dollar is a literal one. Copied through WITH its backslash, so
       Markdown is the thing that unescapes it -- exactly as it would without any
       math support. */
    if (c === '\\' && i + 1 < n && src[i + 1] === '$') { out += src.slice(i, i + 2); i += 2; continue; }

    /* A fenced code block: copy to the closing fence, or to the end if there is
       none. Nothing inside is math. */
    if (src.startsWith('```', i)) {
      const close = src.indexOf('```', i + 3);
      const end = close < 0 ? n : close + 3;
      out += src.slice(i, end); i = end; continue;
    }

    /* An inline code span, closed by a run of backticks of the SAME length --
       which is how ``a ` b`` holds a literal backtick. */
    if (c === '`') {
      let run = 0;
      while (i + run < n && src[i + run] === '`') run++;
      const fence = '`'.repeat(run);
      const close = src.indexOf(fence, i + run);
      const end = close < 0 ? n : close + run;
      out += src.slice(i, end); i = end; continue;
    }

    if (src.startsWith('$$', i)) {
      const close = src.indexOf('$$', i + 2);
      if (close > i + 2) {
        spans.push({ display: true, tex: src.slice(i + 2, close) });
        out += token(spans.length - 1); i = close + 2; continue;
      }
      out += '$$'; i += 2; continue;          /* unclosed: literal */
    }

    if (c === '$') {
      const next = src[i + 1];
      if (next !== undefined && !/\s/.test(next)) {
        /* Scan for a closing `$` preceded by a non-space, skipping escaped ones. */
        let j = i + 1;
        let found = -1;
        while (j < n) {
          if (src[j] === '\\') { j += 2; continue; }
          if (src[j] === '$' && !/\s/.test(src[j - 1])) { found = j; break; }
          if (src[j] === '\n' && src[j + 1] === '\n') break;   /* math does not cross a blank line */
          j++;
        }
        if (found > i + 1) {
          spans.push({ display: false, tex: src.slice(i + 1, found) });
          out += token(spans.length - 1); i = found + 1; continue;
        }
      }
      out += '$'; i += 1; continue;           /* not math: literal */
    }

    out += c; i += 1;
  }

  return { text: out, spans, token };
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
/** Inline math. The fifth marker, and the one the rendering pass above exists for:
 *  a button that writes `$…$` is how you find out the cell renders it. */
export const PROSE_MATH   = { before: '$',  after: '$',  caret: { collapsed: 1, selected: 0 } };
