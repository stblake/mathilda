// plotly.js-dist-min ships a prebuilt bundle with no type declarations, and
// there is no @types package for the -dist-min variant. Output.svelte imports it
// dynamically and uses only newPlot/purge, so an untyped module is honest.
declare module 'plotly.js-dist-min';
