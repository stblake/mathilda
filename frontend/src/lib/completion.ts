/*
 * completion.ts — autocomplete for the notebook's code cells.
 *
 * Offers two kinds of candidate in CodeMirror's completion dropdown (a vertical
 * list when there is more than one match):
 *
 *   1. BUILTINS — every documented builtin symbol. The names come from the
 *      reference-page index (`public/refpages/index.json`, ~1000+ entries),
 *      already fetched and cached at startup by `loadRefpageIndex()` in
 *      refpages.ts. The index value is `<category>/<Name>.md`, so the category
 *      doubles as the completion's detail line.
 *
 *   2. SESSION SYMBOLS — identifiers the user has used or defined anywhere on the
 *      canvas, harvested live from every code cell's source. Definitions are
 *      typed, so scanning the sources catches both "used" and "defined" without a
 *      kernel round-trip. Builtins are removed from this set so a name never
 *      appears twice.
 *
 * Wired into the editor in CellShell.svelte's initEditor(), which only runs for
 * `cell.type === 'code'` — so completion is code-cell-only with no extra gate.
 */

import {
  autocompletion,
  acceptCompletion,
  startCompletion,
  type Completion,
  type CompletionContext,
  type CompletionResult,
} from '@codemirror/autocomplete';
import { tooltips, type EditorView, type KeyBinding } from '@codemirror/view';
import type { Extension } from '@codemirror/state';
import { get } from 'svelte/store';
import { canvasState } from './canvas';
import { autocompleteEnabled } from './properties';
import { loadRefpageIndex } from './refpages';

/* A Mathilda identifier: a letter or $ followed by letters, digits or $.
 * Kept in step with `identifierAt` in refpages.ts. The /g twin is for scanning a
 * whole cell; the plain one is anchored at the cursor by context.matchBefore. */
const IDENT = /[A-Za-z$][A-Za-z0-9$]*/;
const IDENT_SCAN = /[A-Za-z$][A-Za-z0-9$]*/g;

/* Builtins change only when the app is rebuilt, so the list is built once from
 * the index and cached for the life of the page. */
let builtinOptions: Completion[] | null = null;
let builtinNames: Set<string> | null = null;

/** `elementary-functions/Sin.md` -> `elementary functions`. */
function categoryOf(path: string): string {
  return (path.split('/')[0] ?? '').replace(/-/g, ' ');
}

/** The builtin completions and the set of their names, built once. */
async function builtins(): Promise<{ options: Completion[]; names: Set<string> }> {
  if (builtinOptions && builtinNames) {
    return { options: builtinOptions, names: builtinNames };
  }
  const index = await loadRefpageIndex();
  const names = new Set<string>();
  const options: Completion[] = [];
  for (const [name, path] of Object.entries(index)) {
    names.add(name);
    options.push({ label: name, type: 'function', detail: categoryOf(path) });
  }
  builtinOptions = options;
  builtinNames = names;
  return { options, names };
}

/** Identifiers used or defined across every code cell on the canvas, minus the
 *  builtins (offered separately) and the word currently being typed. */
function sessionSymbols(builtinSet: Set<string>, current: string): Completion[] {
  const seen = new Set<string>();
  const out: Completion[] = [];
  let notebooks;
  try {
    notebooks = get(canvasState).notebooks;
  } catch {
    return out; /* canvas not ready — just offer builtins */
  }
  for (const nb of notebooks) {
    let cells;
    try {
      cells = nb.store.allCells();
    } catch {
      continue;
    }
    for (const cell of cells) {
      if (cell.type !== 'code' || !cell.source) continue;
      const tokens = cell.source.match(IDENT_SCAN);
      if (!tokens) continue;
      for (const tok of tokens) {
        if (tok === current || seen.has(tok) || builtinSet.has(tok)) continue;
        seen.add(tok);
        out.push({ label: tok, type: 'variable' });
      }
    }
  }
  return out;
}

/** CodeMirror completion source: builtins + session symbols for the identifier
 *  at the cursor. Session symbols come first so the user's own names win ties. */
export async function mathildaCompletionSource(
  context: CompletionContext,
): Promise<CompletionResult | null> {
  if (!get(autocompleteEnabled)) return null; /* feature off (Properties panel) */
  const word = context.matchBefore(IDENT);
  /* Nothing to complete, and not an explicit Ctrl-Space / Tab request. */
  if (!word || (word.from === word.to && !context.explicit)) return null;

  const { options: builtinOpts, names } = await builtins();
  const session = sessionSymbols(names, word.text);

  return {
    from: word.from,
    options: [...session, ...builtinOpts],
    /* Keep the same result set while the user types more identifier characters
     * (CodeMirror re-filters it itself) instead of re-querying each keystroke. */
    validFor: /^[A-Za-z0-9$]*$/,
  };
}

/** The autocomplete extension: the dropdown (typing-triggered, with the default
 *  Up/Down/Enter/Escape keymap) PLUS a tooltip host on document.body.
 *
 *  The host matters. Cells live inside containers with `overflow: hidden`, and
 *  the canvas applies a `transform: scale()` for zoom. A tooltip rendered in the
 *  editor's own DOM (CodeMirror's default) is therefore clipped by the cell
 *  boundary and positioned relative to the transformed ancestor. Rendering it on
 *  <body> lifts it out of both so the list sits above the cells.
 *
 *  A function, not a const, so document.body is read when a cell's editor is
 *  built (CellShell's initEditor, in onMount) rather than at module load, when
 *  body may not exist yet. CodeMirror flattens the returned array. */
export function mathildaAutocomplete(): Extension {
  return [
    autocompletion({
      override: [mathildaCompletionSource],
      activateOnTyping: true,
      icons: true,
    }),
    tooltips({ parent: document.body }),
  ];
}

/** Tab-completion. If the dropdown is open, Tab accepts the highlighted match.
 *  Otherwise it opens the dropdown ONLY when an identifier character sits just
 *  before the cursor — so Tab keeps its normal meaning (indent / move on) when
 *  there is nothing to complete, e.g. at line start or after whitespace. */
function tabComplete(view: EditorView): boolean {
  if (!get(autocompleteEnabled)) return false; /* off: let Tab behave normally */
  if (acceptCompletion(view)) return true;
  const pos = view.state.selection.main.head;
  const before = view.state.doc.sliceString(Math.max(0, pos - 1), pos);
  if (!/[A-Za-z0-9$]/.test(before)) return false;
  return startCompletion(view);
}

/** Add this BEFORE defaultKeymap so Tab is claimed while a completion is active. */
export const mathildaCompletionKeymap: KeyBinding[] = [
  { key: 'Tab', run: tabComplete },
];
