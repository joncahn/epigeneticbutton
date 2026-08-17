// Behavioral tests for the Genome column's multi-select editor in
// tools/epicc-builder.html.
//
// The builder has no build step and no npm dependencies, so instead of jsdom
// this extracts the functions under test from the shipped HTML and runs them
// against the DOM shim below.
//
// Run via tests/unit/test_builder_genome_picker.py, or directly:
//   node tests/unit/js/genome_picker_test.js
// Exits non-zero on failure.

const fs = require("fs");
const path = require("path");
const vm = require("vm");

const BUILDER = process.env.EPICC_BUILDER_HTML ||
  path.join(__dirname, "..", "..", "..", "tools", "epicc-builder.html");
const src = fs.readFileSync(BUILDER, "utf8");

function extract(name) {
  // Top-level `function name(...) { ... }` up to the closing brace in column 0.
  const m = src.match(new RegExp("^function " + name + "\\([\\s\\S]*?^}", "m"));
  if (!m) throw new Error("could not extract " + name + " from " + BUILDER);
  return m[0];
}

// --------------------------------------------------------------------------
// Minimal DOM shim -- only what the editor actually touches.
// --------------------------------------------------------------------------
class El {
  constructor(tag) {
    this.tagName = tag.toUpperCase();
    this.children = [];
    this.style = {};
    this.attrs = {};
    this.listeners = {};
    this._text = "";
    this.className = "";
    this.checked = false;
    this.value = "";
    this.parent = null;
  }
  set textContent(v) { this._text = v; this.children = []; }
  get textContent() {
    return this._text + this.children.map((c) => c.textContent).join("");
  }
  appendChild(c) { c.parent = this; this.children.push(c); return c; }
  setAttribute(k, v) { this.attrs[k] = v; }
  addEventListener(t, fn) { (this.listeners[t] = this.listeners[t] || []).push(fn); }
  removeEventListener(t, fn) {
    if (this.listeners[t]) this.listeners[t] = this.listeners[t].filter((f) => f !== fn);
  }
  fire(t, ev) {
    (this.listeners[t] || []).slice().forEach((fn) =>
      fn(ev || { preventDefault() {}, stopPropagation() {} }));
  }
  remove() {
    if (this.parent) this.parent.children = this.parent.children.filter((c) => c !== this);
    this.parent = null;
  }
  contains(n) { return n === this || this.children.some((c) => c.contains(n)); }
  focus() {}
  select() {}
  getBoundingClientRect() {
    return { left: 10, top: 10, right: 210, bottom: 40, width: 200, height: 30 };
  }
  find(pred, out = []) {
    if (pred(this)) out.push(this);
    this.children.forEach((c) => c.find(pred, out));
    return out;
  }
}

const docListeners = {};
const document = {
  body: new El("div"),
  createElement: (t) => new El(t),
  addEventListener(t, fn) { (docListeners[t] = docListeners[t] || []).push(fn); },
  removeEventListener(t, fn) {
    if (docListeners[t]) docListeners[t] = docListeners[t].filter((f) => f !== fn);
  },
};
const window = {
  innerWidth: 1400, innerHeight: 900,
  addEventListener() {}, removeEventListener() {},
};

// --------------------------------------------------------------------------
// Load the functions under test
// --------------------------------------------------------------------------
let DEFINED = [];
const sandbox = {
  document, window, console,
  definedGenomes: [],
  getDefinedGenomes: () => DEFINED.slice(),
  module: {},
};
vm.createContext(sandbox);
vm.runInContext(
  extract("parseGenomeList") + "\n" + extract("genomeMultiEditor") +
  "\n;module.exports = { parseGenomeList, genomeMultiEditor };",
  sandbox
);
const { parseGenomeList, genomeMultiEditor } = sandbox.module.exports;

// --------------------------------------------------------------------------
// Harness
// --------------------------------------------------------------------------
let pass = 0;
const failures = [];
function eq(got, want, label) {
  const g = JSON.stringify(got), w = JSON.stringify(want);
  if (g === w) { pass++; console.log("  ok   " + label); }
  else {
    failures.push(label);
    console.log("  FAIL " + label + "\n         got  " + g + "\n         want " + w);
  }
}

function panels() {
  return document.body.children.filter((c) => c.className === "genome-picker");
}
function openEditor(cellValue, defined) {
  DEFINED = defined.slice();
  document.body.children = [];   // isolate tests from each other
  const cellEl = document.body.appendChild(new El("div"));
  const cell = { getValue: () => cellValue, getElement: () => cellEl };
  let committed, cancelled = false, rendered;
  const input = genomeMultiEditor(
    cell,
    (fn) => { rendered = fn; },
    (v) => { committed = v; },
    () => { cancelled = true; }
  );
  rendered();
  return {
    input, cellEl,
    get panel() { return panels().slice(-1)[0]; },
    result: () => ({ committed, cancelled }),
  };
}
const boxesOf = (panel) => panel.find((e) => e.tagName === "INPUT" && e.type === "checkbox");
const labelOf = (box) => box.parent.children.filter((c) => c.tagName === "SPAN")[0].textContent;
const byName = (panel) => {
  const out = {};
  boxesOf(panel).forEach((b) => { out[labelOf(b)] = b; });
  return out;
};
function tick(box, on) { box.checked = on; box.fire("change"); }
function clickButton(panel, text) {
  panel.find((e) => e.tagName === "BUTTON" && e.textContent === text)[0].fire("click");
}
function key(el, k) {
  el.fire("keydown", { key: k, preventDefault() {}, stopPropagation() {} });
}
function fireDocMouseDown(target) {
  (docListeners.mousedown || []).slice().forEach((fn) => fn({ target }));
}

// --------------------------------------------------------------------------
// Tests
// --------------------------------------------------------------------------
console.log("parseGenomeList");
eq(parseGenomeList("B73,W22"), ["B73", "W22"], "splits a pair");
eq(parseGenomeList(" B73 , W22 "), ["B73", "W22"], "trims whitespace");
// read_sample_sheet drops empty entries too, so this must not be an error.
eq(parseGenomeList("B73,,W22"), ["B73", "W22"], "drops stray comma");
eq(parseGenomeList("B73,B73"), ["B73"], "dedupes");
eq(parseGenomeList(""), [], "empty string -> []");
eq(parseGenomeList(null), [], "null -> []");

console.log("\npicker reflects the cell's current value");
{
  const e = openEditor("B73,W22", ["B73", "W22", "TAIR10"]);
  const bs = boxesOf(e.panel);
  eq(bs.map(labelOf), ["B73", "W22", "TAIR10"], "one box per defined genome");
  eq(bs.map((b) => b.checked), [true, true, false], "current selection pre-ticked");
  eq(e.input.value, "B73,W22", "cell input shows the comma list");
}
{
  // A genome present only in this cell must not vanish from the picker.
  const e = openEditor("Spombe", []);
  eq(boxesOf(e.panel).map(labelOf), ["Spombe"], "cell-only genome still listed");
  eq(boxesOf(e.panel)[0].checked, true, "and ticked");
}

console.log("\nticking preserves existing order");
{
  const e = openEditor("W22,B73", ["B73", "W22", "TAIR10"]);
  const b = byName(e.panel);
  tick(b.TAIR10, true);
  eq(e.input.value, "W22,B73,TAIR10", "new tick appended, existing order kept");
  tick(b.B73, false);
  eq(e.input.value, "W22,TAIR10", "untick removes without reshuffling");
  clickButton(e.panel, "Done");
  eq(e.result().committed, "W22,TAIR10", "Done commits");
}

console.log("\nthe text box stays live");
{
  // Keyboard and paste workflows have to keep working, so a typed value wins.
  const e = openEditor("B73", ["B73", "W22"]);
  e.input.value = "TAIR10 , ColCEN";
  clickButton(e.panel, "Done");
  eq(e.result().committed, "TAIR10,ColCEN", "typed value normalized and committed");
}
{
  const e = openEditor("B73", ["B73", "W22"]);
  tick(byName(e.panel).W22, true);
  key(e.input, "Enter");
  eq(e.result().committed, "B73,W22", "Enter commits");
}
{
  const e = openEditor("B73", ["B73", "W22"]);
  key(e.input, "Escape");
  eq(e.result(), { committed: undefined, cancelled: true }, "Escape cancels");
}

console.log("\nadding a genome from the cell");
{
  const e = openEditor("B73", ["B73"]);
  const addInput = e.panel.find(
    (x) => x.tagName === "INPUT" && x.attrs["aria-label"] === "New genome name")[0];
  addInput.value = "NewRef";
  clickButton(e.panel, "Add");
  eq(sandbox.definedGenomes.indexOf("NewRef") >= 0, true, "registered for other rows");
  eq(byName(e.panel).NewRef !== undefined, true, "appears as a box");
  eq(byName(e.panel).NewRef.checked, true, "and is ticked");
  eq(panels().length, 1, "panel rebuilt in place, not stacked");
  clickButton(e.panel, "Done");
  eq(e.result().committed, "B73,NewRef", "commits with the added genome");
}
{
  const e = openEditor("", []);
  eq(e.panel.find((x) => x.className === "genome-picker-empty").length, 1,
     "empty state hints where to add one");
  eq(boxesOf(e.panel).length, 0, "and shows no boxes");
}

console.log("\nclicking away commits rather than discarding");
{
  const e = openEditor("B73", ["B73", "W22"]);
  const w22 = byName(e.panel).W22;
  tick(w22, true);
  // Inside the panel: must be ignored, or ticking a box would end the edit.
  fireDocMouseDown(w22);
  eq(e.result(), { committed: undefined, cancelled: false }, "mousedown in panel ignored");
  fireDocMouseDown(e.cellEl);
  eq(e.result(), { committed: undefined, cancelled: false }, "mousedown on cell ignored");
  fireDocMouseDown(new El("div"));
  eq(e.result().committed, "B73,W22", "mousedown outside commits the edit");
  eq(panels().length, 0, "and removes the panel");
}

console.log("\nteardown leaves nothing behind");
{
  const base = (docListeners.mousedown || []).length;
  const e = openEditor("B73", ["B73"]);
  const during = (docListeners.mousedown || []).length - base;
  clickButton(e.panel, "Done");
  eq([during, (docListeners.mousedown || []).length - base], [1, 0],
     "document listener added then removed on commit");
}
{
  const base = (docListeners.mousedown || []).length;
  const e = openEditor("B73", ["B73"]);
  key(e.input, "Escape");
  eq([(docListeners.mousedown || []).length - base, panels().length], [0, 0],
     "cancel tears down too");
}
{
  // success()/cancel() must fire at most once: Tabulator treats a second call
  // as an edit on a cell that is no longer being edited.
  const e = openEditor("B73", ["B73"]);
  clickButton(e.panel, "Done");
  fireDocMouseDown(new El("div"));   // would re-commit if not guarded
  key(e.input, "Enter");
  eq(e.result().committed, "B73", "repeat commits are no-ops");
}

console.log("\n" + pass + " passed, " + failures.length + " failed");
if (failures.length) process.exit(1);
