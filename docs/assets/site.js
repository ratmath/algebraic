// Shared page furniture for the algebraic docs: the contents column, the
// theme switch, result rendering, expression demos (any element with
// data-expr becomes one) and a small complex-plane plot.

import * as A from "./algebraic.js?v=4213fe0658";
import { run } from "./lang.js?v=4213fe0658";
import { mountRunboxes } from "./runbox.js?v=4213fe0658";

export { A, run };

export const PAGES = [
  { href: "index.html", title: "Overview" },
  { href: "reference.html", title: "Reference" },
  { href: "explore.html", title: "Explore" },
  { group: "Numbers" },
  { href: "numbers.html", title: "Algebraic numbers", blurb: "Exact roots, minimal polynomials, equality" },
  { href: "complex.html", title: "Complex numbers", blurb: "Certified boxes, principal roots, a + bi" },
  { href: "angles.html", title: "Angles and roots of unity", blurb: "cos(π/7), cyclotomic polynomials, radicals" },
  { href: "recognize.html", title: "Recognising numbers", blurb: "From digits to a minimal polynomial (LLL)" },
  { group: "Polynomials" },
  { href: "factoring.html", title: "Factoring", blurb: "Zassenhaus: mod p, Hensel, recombine" },
  { href: "polynomials.html", title: "Polynomial toolkit", blurb: "Resultants, Sturm, partial fractions" },
  { href: "galois.html", title: "Galois groups", blurb: "The group moving the roots, its subgroups and their fields" },
  { group: "Structures" },
  { href: "number-fields.html", title: "Number fields", blurb: "Q(√2), Q(i), Q(∛2): norms and traces" },
  { href: "finite-fields.html", title: "Finite fields", blurb: "GF(pᵏ), generators, logs, mod p" },
];

export const h = (tag, attrs = {}, ...kids) => {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === "class") el.className = v;
    else if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
    else el.setAttribute(k, v === true ? "" : v);
  }
  for (const kid of kids.flat(Infinity)) if (kid !== null && kid !== undefined && kid !== false) el.append(kid instanceof Node ? kid : document.createTextNode(String(kid)));
  return el;
};

// ---- chrome -------------------------------------------------------------------

function theme() {
  let saved = null;
  try {
    saved = localStorage.getItem("ratmath-theme");
  } catch {
    /* storage may be blocked */
  }
  const apply = (t) => {
    if (t === "auto") document.documentElement.removeAttribute("data-theme");
    else document.documentElement.setAttribute("data-theme", t);
    try {
      localStorage.setItem("ratmath-theme", t);
    } catch {
      /* ignore */
    }
    for (const b of document.querySelectorAll(".theme button")) b.classList.toggle("on", b.dataset.t === t);
  };
  if (saved) apply(saved);
  return h("div", { class: "theme", role: "group", "aria-label": "Theme" },
    ["auto", "light", "dark"].map((t) => h("button", { "data-t": t, class: (saved ?? "auto") === t ? "on" : "", onclick: () => apply(t) }, t)));
}

function chrome() {
  const here = location.pathname.split("/").pop() || "index.html";
  const aside = document.querySelector("aside");
  if (!aside) return;
  const nav = h("nav", { class: "toc", "aria-label": "Contents" });
  const sel = h("select", { id: "toc-select", "aria-label": "Page", onchange: (e) => (location.href = e.target.value) });
  let group = null;
  for (const p of PAGES) {
    if (p.group) {
      nav.append(h("div", { class: "group" }, p.group));
      group = h("optgroup", { label: p.group });
      sel.append(group);
      continue;
    }
    nav.append(h("a", { href: p.href, class: p.href === here ? "on" : null, "aria-current": p.href === here ? "page" : null }, p.title));
    (group ?? sel).append(h("option", { value: p.href, selected: p.href === here }, p.title));
  }
  aside.append(
    h("header", { class: "brand" }, h("a", { href: "index.html" }, h("div", { class: "mark" }, "RatMath ", h("span", {}, "algebraic"))), h("small", {}, "Exact algebraic numbers")),
    sel, nav, theme(),
  );
}

/** Topic cards for the overview page. */
export function topicCards(el) {
  el.append(...PAGES.filter((p) => p.blurb).map((p) => h("a", { href: p.href }, h("strong", {}, p.title), h("span", {}, p.blurb))));
}

// ---- results --------------------------------------------------------------------

const fmtP = (c) => A.qpoly.formatP(c);

/** One result descriptor (from lang.run) as DOM. */
export function render(r, { detail = true } = {}) {
  if (!r) return h("span", { class: "muted" }, "—");
  switch (r.type) {
    case "number": {
      const v = r.value;
      const main = A.format(v);
      const wrap = h("div", {}, h("div", { class: "big" }, main));
      if (r.note) wrap.append(h("div", { class: "small muted" }, r.note));
      if (detail) {
        const rows = [];
        const exact = !main.endsWith("…");
        const rat = A.rationalOf(v);
        if (!rat) {
          if (exact) rows.push(["decimal", A.isReal(v) ? v.toString(30) : v.toDecimal(30)]);
          const mp = A.minimalPolynomial(v);
          rows.push(["minimal polynomial", A.qpoly.formatP(mp.map((c) => [c, 1n]))], ["degree", String(mp.length - 1)]);
          rows.push(["kind", A.isReal(v) ? "real algebraic" : "complex algebraic"]);
        }
        if (rows.length) wrap.append(h("dl", { class: "kv" }, rows.flatMap(([k, x]) => [h("dt", {}, k), h("dd", {}, x)])));
      }
      return wrap;
    }
    case "poly": return h("div", { class: "big" }, fmtP(r.c));
    case "text": return h("div", { class: "big" }, r.text);
    case "list": {
      const box = h("div", {}, r.label ? h("div", { class: "small muted" }, `${r.items.length} ${r.label}`) : null);
      const tbl = h("table", { class: "data" });
      r.items.forEach((it, k) => tbl.append(h("tr", {}, h("td", { class: "muted" }, String(k + 1)), h("td", {}, render(it, { detail: false })))));
      box.append(h("div", { class: "tablewrap" }, tbl));
      return box;
    }
    case "record": return h("div", {}, h("div", { class: "big" }, r.title), h("dl", { class: "kv" }, r.rows.flatMap(([k, x]) => [h("dt", {}, k), h("dd", {}, x)])));
    default: return h("pre", {}, JSON.stringify(r));
  }
}

/** Run `work` after the browser has painted a busy state; shows the time it took. */
export function later(out, work) {
  out.classList.add("busy");
  setTimeout(() => {
    const t0 = performance.now();
    let node;
    try {
      node = work();
      out.classList.remove("err");
    } catch (e) {
      node = h("div", { class: "err" }, e.message);
      out.classList.add("err");
    }
    const ms = performance.now() - t0;
    out.replaceChildren(node, h("div", { class: "timing" }, ms < 1 ? "< 1 ms" : `${ms.toFixed(0)} ms`));
    out.classList.remove("busy");
  }, 20);
}

/** An expression demo: input, examples, result. */
export function exprDemo(el, { expr, examples = [], title = null, rows = 1 }) {
  const env = {};
  const input = rows > 1
    ? h("textarea", { class: "text", rows, spellcheck: "false", "aria-label": "Expression" }, expr)
    : h("input", { class: "text grow", value: expr, spellcheck: "false", "aria-label": "Expression" });
  const out = h("div", { class: "out" });
  const go = () => later(out, () => {
    const lines = input.value.split(/\n|;/).map((s) => s.trim()).filter(Boolean);
    let last = null;
    for (const line of lines) last = run(line, env).result;
    return render(last);
  });
  input.addEventListener("keydown", (e) => {
    if (e.key === "Enter" && (rows === 1 || e.metaKey || e.ctrlKey)) {
      e.preventDefault();
      go();
    }
  });
  el.classList.add("demo");
  el.replaceChildren(...[
    title ? h("div", { class: "demo-title" }, title) : null,
    h("div", { class: "demo-row" }, input, h("button", { onclick: go }, "Run")),
    examples.length ? h("div", { class: "demo-row" }, examples.map((x) => h("button", { class: "chipbtn", onclick: () => { input.value = x; go(); } }, x))) : null,
    out,
  ].filter(Boolean));
  go();
}

// ---- the complex plane ----------------------------------------------------------

const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim() || "#888";

/** Plot points { re, im, label?, ring? } on a square canvas centred at 0; draws the unit circle when asked. */
export function plane(canvas, points, { unit = false, R = null } = {}) {
  const dpr = window.devicePixelRatio || 1;
  const W = canvas.clientWidth || 360;
  canvas.width = Math.round(W * dpr);
  canvas.height = Math.round(W * dpr);
  const g = canvas.getContext("2d");
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, W, W);
  const r = R ?? Math.max(1.2, ...points.map((p) => Math.max(Math.abs(p.re), Math.abs(p.im)) * 1.25));
  const sx = (x) => W / 2 + (x / r) * (W / 2 - 12);
  const sy = (y) => W / 2 - (y / r) * (W / 2 - 12);
  g.strokeStyle = css("--line");
  g.lineWidth = 1;
  g.beginPath();
  g.moveTo(0, sy(0)); g.lineTo(W, sy(0));
  g.moveTo(sx(0), 0); g.lineTo(sx(0), W);
  g.stroke();
  if (unit) {
    g.strokeStyle = css("--muted");
    g.setLineDash([4, 4]);
    g.beginPath();
    g.arc(sx(0), sy(0), sx(1) - sx(0), 0, 2 * Math.PI);
    g.stroke();
    g.setLineDash([]);
  }
  g.font = `11px ${css("--font-mono")}`;
  g.fillStyle = css("--muted");
  g.fillText(r.toPrecision(2), sx(r) - 26, sy(0) + 14);
  for (const p of points) {
    g.fillStyle = p.real ? css("--exact") : css("--accent");
    g.beginPath();
    g.arc(sx(p.re), sy(p.im), 5, 0, 2 * Math.PI);
    g.fill();
    if (p.label) {
      g.fillStyle = css("--fg");
      g.fillText(p.label, sx(p.re) + 7, sy(p.im) - 7);
    }
  }
}

/** Float coordinates of an exact number, for plotting. */
export function coords(v) {
  if (A.isReal(v)) return { re: v.toNumber(), im: 0, real: true };
  const { re, im } = v.toNumbers();
  return { re, im, real: false };
}

// ---- boot -------------------------------------------------------------------------

chrome();
for (const el of document.querySelectorAll("[data-expr]")) {
  exprDemo(el, {
    expr: el.dataset.expr.replace(/\\n/g, "\n"),
    examples: el.dataset.examples ? el.dataset.examples.split("|") : [],
    title: el.dataset.title ?? null,
    rows: el.dataset.rows ? Number(el.dataset.rows) : 1,
  });
}
mountRunboxes(document, { scope: A });
