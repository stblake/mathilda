// properties.ts — state for the properties sidebar.
//
// Two stores, and they are here rather than in theme.ts or status.ts because
// neither is about a theme or a measurement: one is whether the panel is open,
// the other is a display preference the panel owns.
//
// Plain writables with no persistence, matching darkMode in theme.ts. Nothing in
// this app persists UI state yet, and adding localStorage for one preference
// would make this the only setting that survives a restart -- a surprise, not a
// feature. When persistence arrives it should arrive for all of them at once.

import { writable } from 'svelte/store';

/** Whether the properties sidebar is showing. Toolbar's Sidebar group toggles it. */
export const propertiesOpen = writable(false);

/** Whether a code cell shows its `In[n]` label once it has been evaluated.
 *
 *  Lives here rather than in the toolbar because it is a persistent preference
 *  about how the notebook reads, not an action -- the toolbar is for verbs. */
export const showExecLabels = writable(true);
