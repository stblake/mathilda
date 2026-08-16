#!/usr/bin/env python3
"""check_menu_ids.py -- does every native menu item actually do something?

WHY THIS EXISTS. The menu bar is built in Rust (`frontend/src-tauri/src/lib.rs`) and its items are
handled in TypeScript (`frontend/src/lib/menuCommands.ts`), joined only by a string id travelling
through a `menu:<id>` event. Nothing in either language checks the other, and both failure
directions are SILENT:

  * an id in the Rust menu with no case in the dispatcher is a menu item that does nothing when
    clicked -- no error, no log, just a dead command;
  * a case in the dispatcher whose id no menu item emits is dead code that reads like wiring, which
    is how three Edit handlers survived after the real work moved to predefined native items.

Both were present the first time this was checked by hand. A menu is exactly the kind of surface
nobody exercises item by item, so the check belongs in a script rather than in a habit.

It also checks MENU_IDS, the list App.svelte subscribes to: an id that is handled but not listed is
never subscribed, so the handler still never runs.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RUST = ROOT / "frontend" / "src-tauri" / "src" / "lib.rs"
TS = ROOT / "frontend" / "src" / "lib" / "menuCommands.ts"

# Ids are alphanumeric with dots and dashes; DIGITS MATTER -- a pattern of [a-zA-Z.-] only silently
# skipped `gfx.image3d` and reported it as missing from a list it was in.
ID_CHARS = r"[A-Za-z0-9.\-]+"


def main():
    if not RUST.exists() or not TS.exists():
        print("menu sources not found; skipping")
        return 0
    rust = RUST.read_text()
    ts = TS.read_text()

    native = set(re.findall(r'with_id\(app,\s*"(' + ID_CHARS + r')"', rust))
    try:
        seg = ts[ts.index("export const MENU_IDS") : ts.index("] as const")]
    except ValueError:
        print("MENU_IDS not found in menuCommands.ts")
        return 1
    listed = set(re.findall(r"'(" + ID_CHARS + r")'", seg))
    cases = set(re.findall(r"case '(" + ID_CHARS + r")':", ts))

    if not native or not listed or not cases:
        # Any of the three coming back empty means a pattern stopped matching, not that a menu is
        # empty. Reporting success on that would be a green light with nothing behind it.
        print(
            "one of the three sets is empty (native=%d listed=%d cases=%d);"
            " the patterns no longer match the sources"
            % (len(native), len(listed), len(cases))
        )
        return 1

    problems = []
    for label, missing, why in [
        ("native item with no handler", native - cases, "clicking it does nothing"),
        (
            "handler for an id nothing emits",
            cases - native,
            "dead code that reads like wiring",
        ),
        (
            "handled but not in MENU_IDS",
            cases - listed,
            "App.svelte never subscribes, so the handler never runs",
        ),
        (
            "in MENU_IDS but not native",
            listed - native,
            "a subscription for an event that is never sent",
        ),
    ]:
        if missing:
            problems.append((label, sorted(missing), why))

    print(
        "%d native menu ids, %d subscribed, %d handled"
        % (len(native), len(listed), len(cases))
    )
    if not problems:
        print("menu wiring is complete in both directions")
        return 0
    for label, ids, why in problems:
        print("\n%s (%s):" % (label, why))
        print("  " + ", ".join(ids))
    return 1


if __name__ == "__main__":
    sys.exit(main())
