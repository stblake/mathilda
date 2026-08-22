// search.ts — find across every cell of the focused notebook.
//
// NOT @codemirror/search. That package's openSearchPanel searches ONE editor --
// the cell that happens to hold focus -- and a find bar that silently ignores the
// other forty cells while calling itself notebook search is worse than no find
// bar, because you would trust its "no matches". So the matching runs over the
// notebook's own model, `store.allCells()`, and navigation then drives whichever
// editor owns the match it landed on.
//
// The matcher is a pure function of (cells, query, caseSensitive) with no store
// and no DOM, which is what lets `npm run check:search` run the real thing rather
// than a paraphrase of it.

import { writable } from 'svelte/store';

/** One occurrence, as offsets into that cell's own source. */
export type SearchMatch = {
  cellId: string;
  /** Offset of the first character of the match within the cell's source. */
  start: number;
  /** Offset one past the last character. */
  end: number;
};

/** The subset of a cell this module needs, so the matcher is testable without
 *  building a whole Cell. */
export type SearchableCell = { id: string; source: string };

export const searchOpen = writable(false);
export const searchQuery = writable('');
export const searchCaseSensitive = writable(false);
/** Index into the current match list. Clamped by the bar, not here. */
export const searchIndex = writable(0);

/** Every occurrence of `query`, in document order: cells in notebook order, and
 *  within a cell by increasing offset.
 *
 *  Matches do NOT overlap. Searching "aa" in "aaaa" gives two matches at 0 and 2,
 *  not three at 0, 1 and 2 -- overlapping hits would make "next match" step
 *  through positions that look identical to the reader, and every editor's find
 *  behaves this way.
 *
 *  An empty or whitespace-only query finds nothing rather than matching at every
 *  offset. A query of a single space is a legitimate search, so only the empty
 *  string is rejected; "  " finds real double spaces.
 */
export function findMatches(cells: SearchableCell[], query: string,
                            caseSensitive: boolean): SearchMatch[] {
  const out: SearchMatch[] = [];
  if (!query) return out;

  const needle = caseSensitive ? query : query.toLowerCase();
  for (const cell of cells) {
    const src = cell.source ?? '';
    const hay = caseSensitive ? src : src.toLowerCase();
    let from = 0;
    for (;;) {
      const at = hay.indexOf(needle, from);
      if (at < 0) break;
      out.push({ cellId: cell.id, start: at, end: at + query.length });
      from = at + query.length;   /* non-overlapping */
    }
  }
  return out;
}

/** Step `index` by `delta` within `count` matches, wrapping in both directions.
 *
 *  Separate from the component because off-by-one and negative-modulo are exactly
 *  the two things worth checking: in JavaScript `-1 % 3` is `-1`, not `2`, so the
 *  naive form sends Shift+Enter at the first match to a negative index and the
 *  bar shows "0 of 3". */
export function stepIndex(index: number, delta: number, count: number): number {
  if (count <= 0) return 0;
  return ((index + delta) % count + count) % count;
}
