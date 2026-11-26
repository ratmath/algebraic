// The Galois explorer: a polynomial's group acting on its roots, the subgroup
// lattice with the fixed field of each subgroup, and factor patterns mod p.

import { A, run, h, coords } from "./site.js?v=b16f7dfc21";

const SVG = "http://www.w3.org/2000/svg";
const s = (tag, attrs = {}, ...kids) => {
  const el = document.createElementNS(SVG, tag);
  for (const [k, v] of Object.entries(attrs)) if (v !== null && v !== undefined) el.setAttribute(k, v);
  for (const kid of kids) el.append(kid instanceof Node ? kid : document.createTextNode(String(kid)));
  return el;
};
const fmt = (p) => A.polys.format(p);
const subscript = (n) => String(n).split("").map((d) => "₀₁₂₃₄₅₆₇₈₉"[d]).join("");

const EXAMPLES = [
  ["x^3 - 2", "S₃"], ["x^3 - 3x + 1", "C₃"], ["x^4 - 2", "D₄"], ["x^4 + 1", "V₄"], ["x^4 - 4x^2 + 2", "C₄"], ["x^4 + 8x + 12", "A₄"], ["x^4 + x + 1", "S₄"],
  ["x^5 - 2", "F₂₀"], ["x^5 - 5x + 12", "D₅"], ["x^5 - x - 1", "S₅"], ["x^6 - 2", "D₆"], ["x^6 + x^4 - 2x^2 - 1", "A₄"], ["x^7 - 2", "F₄₂"], ["x^7 - 7x + 3", "PSL(3,2)"],
];

function parse(text) {
  const { result } = run(text);
  if (!result || result.type !== "poly") throw new Error("Type a polynomial in x with rational coefficients, such as x^4 - 2.");
  return result.c;
}

/** The roots in the plane, numbered; `perm` draws an arrow from each root to its image. */
function rootPlot(action, perm, highlight = null) {
  const pts = action.roots.map(coords);
  const R = Math.max(1.2, ...pts.map((p) => Math.hypot(p.re, p.im))) * 1.3;
  const W = 340;
  const x = (v) => W / 2 + (v / R) * (W / 2);
  const y = (v) => W / 2 - (v / R) * (W / 2);
  const svg = s("svg", { viewBox: `0 0 ${W} ${W}`, class: "gal-plot", role: "img", "aria-label": "The roots in the complex plane" });
  svg.append(
    s("defs", {}, s("marker", { id: "gal-arrow", viewBox: "0 0 10 10", refX: 9, refY: 5, markerWidth: 7, markerHeight: 7, orient: "auto-start-reverse" }, s("path", { d: "M0,0 L10,5 L0,10 z", class: "gal-arrowhead" }))),
    s("line", { x1: 0, y1: W / 2, x2: W, y2: W / 2, class: "gal-axis" }),
    s("line", { x1: W / 2, y1: 0, x2: W / 2, y2: W, class: "gal-axis" }),
  );
  if (perm) {
    perm.forEach((to, from) => {
      if (to === from) return;
      const [x1, y1, x2, y2] = [x(pts[from].re), y(pts[from].im), x(pts[to].re), y(pts[to].im)];
      const dx = x2 - x1;
      const dy = y2 - y1;
      const len = Math.hypot(dx, dy) || 1;
      // a gentle bow, so that a swap shows as two arrows; stop short of the dots
      const bow = Math.min(40, len * 0.22);
      const cx = (x1 + x2) / 2 - (dy / len) * bow;
      const cy = (y1 + y2) / 2 + (dx / len) * bow;
      const trim = (px, py, qx, qy, d) => { const l = Math.hypot(qx - px, qy - py) || 1; return [px + ((qx - px) / l) * d, py + ((qy - py) / l) * d]; };
      const [sx, sy] = trim(x1, y1, cx, cy, 11);
      const [ex, ey] = trim(x2, y2, cx, cy, 12);
      svg.append(s("path", { d: `M${sx},${sy} Q${cx},${cy} ${ex},${ey}`, class: "gal-move", "marker-end": "url(#gal-arrow)" }));
    });
  }
  pts.forEach((p, i) => {
    const on = highlight && highlight.has(i);
    svg.append(
      s("circle", { cx: x(p.re), cy: y(p.im), r: 9, class: `gal-root${p.real ? " real" : ""}${on ? " on" : ""}` }),
      s("text", { x: x(p.re), y: y(p.im) + 3.5, class: "gal-root-label" }, i + 1),
    );
  });
  return svg;
}

/** The subgroup classes as a diagram: the whole group on top, the trivial group at the bottom. */
function lattice(rows, selected, pick) {
  const orders = [...new Set(rows.map((r) => r.order))].sort((a, b) => b - a);
  const widest = Math.max(...orders.map((o) => rows.filter((r) => r.order === o).length));
  const W = Math.max(340, widest * 58);
  const H = Math.max(120, orders.length * 62);
  const pos = rows.map((r) => {
    const level = orders.indexOf(r.order);
    const same = rows.filter((q) => q.order === r.order);
    const k = same.indexOf(r);
    return { x: ((k + 1) / (same.length + 1)) * W, y: 26 + (level / Math.max(1, orders.length - 1)) * (H - 52) };
  });
  const svg = s("svg", { viewBox: `0 0 ${W} ${H}`, class: "gal-lattice", style: `max-width:${Math.round(W * 1.15)}px`, role: "img", "aria-label": "Subgroups up to conjugacy, ordered by containment" });
  rows.forEach((r, i) => {
    // draw only the covers: containments with nothing in between
    const covers = r.above.filter((j) => !r.above.some((m) => m !== j && rows[m].above.includes(j)));
    for (const j of covers) svg.append(s("line", { x1: pos[i].x, y1: pos[i].y, x2: pos[j].x, y2: pos[j].y, class: "gal-edge" }));
  });
  rows.forEach((r, i) => {
    const g = s("g", { class: `gal-node${r.normal ? " normal" : ""}${i === selected ? " on" : ""}`, tabindex: 0, role: "button", "aria-label": `Subgroup of order ${r.order}` });
    g.append(s("circle", { cx: pos[i].x, cy: pos[i].y, r: 15 }), s("text", { x: pos[i].x, y: pos[i].y + 4 }, r.order));
    g.addEventListener("click", () => pick(i));
    g.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); pick(i); } });
    svg.append(g);
  });
  return svg;
}

function mount(root) {
  const input = h("input", { class: "text grow", value: "x^4 - 2", spellcheck: "false", "aria-label": "Polynomial" });
  const go = h("button", {}, "Find the group");
  const out = h("div", { class: "gal-out" });
  const chips = h("div", { class: "demo-row" }, EXAMPLES.map(([p, name]) => h("button", { class: "chipbtn", title: name, onclick: () => { input.value = p; compute(); } }, p)));
  root.append(h("div", { class: "demo" }, h("div", { class: "demo-row" }, input, go), chips, out));
  go.addEventListener("click", () => compute());
  input.addEventListener("keydown", (e) => { if (e.key === "Enter") compute(); });

  function compute() {
    out.replaceChildren(h("p", { class: "muted" }, "Working…"));
    setTimeout(() => {
      try {
        show(parse(input.value));
      } catch (e) {
        out.replaceChildren(h("p", { class: "err" }, e.message));
      }
    }, 20);
  }

  function show(coeffs) {
    let action;
    try {
      action = A.galoisAction(coeffs);
    } catch (e) {
      const g = e instanceof RangeError ? A.galoisGroup(coeffs) : null;
      if (!g?.reducible) throw e;
      out.replaceChildren(
        h("p", {}, "This polynomial factors, so its roots fall into separate families. Each irreducible factor has its own group:"),
        h("div", { class: "tablewrap" }, h("table", { class: "data" }, h("tr", {}, h("th", {}, "factor"), h("th", {}, "group"), h("th", {}, "order")),
          g.factors.map((f) => h("tr", {}, h("td", {}, fmt(f.factor)), h("td", {}, f.group.name), h("td", {}, f.group.order))))),
        h("p", { class: "small muted" }, "Pick one factor to explore it."),
      );
      return;
    }
    const G = action.group;
    const conj = action.conjugation;
    const series = action.derivedSeries();

    // --- summary
    const summary = h("dl", { class: "kv" },
      h("dt", {}, "polynomial"), h("dd", {}, fmt(action.polynomial)),
      h("dt", {}, "group"), h("dd", {}, `${action.name}  (${action.label}), ${action.order} permutations of the ${action.degree} roots`),
      h("dt", {}, "generated by"), h("dd", {}, G.generators.map(A.cycleString).join(",  ") || "nothing: only the identity"),
      h("dt", {}, "by radicals?"), h("dd", {}, action.solvable
        ? `yes. Taking commutators again and again shrinks the group to nothing: ${series.join(" → ")}`
        : `no. The commutators stop at a group of order ${series[series.length - 1]} and never reach 1: ${series.join(" → ")}`),
      h("dt", {}, "discriminant"), h("dd", {}, `${action.discriminant}${action.discriminantIsSquare ? ", a square: every permutation is even" : ", not a square: some permutations are odd"}`),
    );

    // --- the action on the roots
    const plotBox = h("div", {});
    const caption = h("p", { class: "small" });
    const elements = h("div", { class: "gal-elements" });
    let current = null;
    const draw = (p, highlight) => {
      current = p;
      plotBox.replaceChildren(rootPlot(action, p, highlight));
      const moved = p ? p.map((to, from) => (to === from ? null : `x${subscript(from + 1)} → x${subscript(to + 1)}`)).filter(Boolean) : [];
      const isConj = p && p.join() === conj.join();
      caption.textContent = !p ? "Pick a permutation to see where it sends each root."
        : `${A.cycleString(p)}: ${moved.length ? moved.join(", ") : "every root stays put"}.${isConj && moved.length ? " This one is complex conjugation." : ""}`;
      for (const b of elements.children) b.classList.toggle("on", b.dataset.p === (p ? p.join() : ""));
    };
    const shown = G.elements.slice(0, 120);
    for (const p of shown) {
      const isConj = p.join() === conj.join() && A.cycleString(p) !== "()";
      elements.append(h("button", { class: `chipbtn${isConj ? " conj" : ""}`, "data-p": p.join(), title: isConj ? "complex conjugation" : null, onclick: () => draw(p) }, A.cycleString(p)));
    }
    const more = G.order > shown.length ? h("p", { class: "small muted" }, `Showing ${shown.length} of the ${G.order} permutations.`) : null;
    const actionPart = h("div", { class: "gal-split" },
      h("div", {}, plotBox, caption),
      h("div", {}, h("p", { class: "small muted" }, "Roots are numbered: real ones from left to right, then the others. Each button is one permutation that keeps every polynomial relation among the roots true."), elements, more),
    );
    draw(conj.some((x, i) => x !== i) ? conj : G.elements[1] ?? null);

    // --- the correspondence
    const corr = h("div", {}, h("p", { class: "muted" }, "Working out the subgroups and their fixed fields…"));
    // --- Frobenius
    const frob = h("div", {});
    const count = h("input", { type: "range", min: 50, max: 1000, step: 50, value: 200, "aria-label": "Number of primes" });
    const countLabel = h("span", { class: "small muted" });
    const drawFrob = () => {
      const { primes, patterns } = action.frobenius(Number(count.value));
      countLabel.textContent = `${primes} primes`;
      const top = Math.max(...patterns.map((x) => Math.max(x.observed, x.expected)));
      frob.replaceChildren(h("div", { class: "tablewrap" }, h("table", { class: "data gal-frob" },
        h("tr", {}, h("th", {}, "degrees of the factors mod p"), h("th", {}, "first primes"), h("th", {}, "primes"), h("th", {}, "share seen"), h("th", {}, "share of the group"), h("th", {}, "")),
        patterns.map((x) => h("tr", {},
          h("td", {}, x.type.join(" + ")),
          h("td", {}, x.first.join(", ") || "none yet"),
          h("td", {}, x.primes),
          h("td", {}, `${(100 * x.observed).toFixed(1)}%`),
          h("td", {}, `${x.elements}/${G.order} = ${(100 * x.expected).toFixed(1)}%`),
          h("td", {}, h("div", { class: "gal-bars" }, h("i", { class: "seen", style: `width:${(100 * x.observed) / top}%` }), h("i", { class: "want", style: `width:${(100 * x.expected) / top}%` }))),
        )))));
    };
    count.addEventListener("input", drawFrob);

    out.replaceChildren(
      summary,
      h("h3", {}, "The group moving the roots"), actionPart,
      h("h3", {}, "Subgroups and the fields they fix"), corr,
      h("h3", {}, "The same group, seen modulo primes"),
      h("p", { class: "small" }, "Factor the polynomial modulo a prime p. The degrees of the factors are the cycle lengths of a permutation in the group (Dedekind), and over many primes each pattern turns up as often as the group has permutations of that shape (Chebotarev). The upper bar is what was seen, the lower one what the group predicts."),
      h("div", { class: "demo-row" }, count, countLabel), frob,
    );
    drawFrob();

    setTimeout(() => {
      try {
        if (G.order > 200) {
          corr.replaceChildren(h("p", { class: "small muted" }, `This group has ${G.order} elements, too many to list every subgroup here. In code, action.fixedField(H) still gives the field fixed by any subgroup H you choose.`));
          return;
        }
        const rows = action.correspondence({ maxDegree: 24 });
        const detail = h("div", { class: "out" });
        const diagram = h("div", { class: "tablewrap" });
        const pick = (i) => {
          const r = rows[i];
          diagram.replaceChildren(lattice(rows, i, pick));
          const H = r.subgroup;
          const fixedRoots = new Set(action.roots.map((_, k) => k).filter((k) => H.elements.every((g) => g[k] === k)));
          plotBox.replaceChildren(rootPlot(action, current, fixedRoots));
          detail.replaceChildren(h("dl", { class: "kv" },
            h("dt", {}, "subgroup"), h("dd", {}, `${H.generators.length ? H.generators.map(A.cycleString).join(",  ") : "only the identity"}   (order ${r.order})`),
            h("dt", {}, "copies"), h("dd", {}, r.normal ? "1: it is a normal subgroup, so its field is itself a Galois extension of Q" : `${r.conjugates} conjugate subgroups, one for each way of relabelling the roots`),
            h("dt", {}, "fixed field"), h("dd", {}, r.index === 1 ? "Q itself: only rational numbers stay put under the whole group"
              : r.field ? `degree ${r.index} over Q, generated by a root of ${fmt(r.field.polynomial)}` : `degree ${r.index} over Q (too large to write down here)`),
            ...(r.field && r.index > 1 ? [h("dt", {}, "that root is"), h("dd", {}, `${r.field.generator}${r.field.substitution ? `, with ${r.field.substitution}` : ""}`)] : []),
            ...(fixedRoots.size ? [h("dt", {}, "roots it fixes"), h("dd", {}, [...fixedRoots].map((k) => `x${subscript(k + 1)}`).join(", "))] : []),
          ));
        };
        corr.replaceChildren(
          h("p", { class: "small" }, `Each circle is a subgroup (up to relabelling the roots), marked with its size; the whole group is on top. A subgroup of size m fixes a field of degree ${G.order}/m over Q, and bigger subgroups fix smaller fields. Filled circles are normal subgroups. Pick one.`),
          diagram, detail,
        );
        pick(rows.length - 1);
      } catch (e) {
        corr.replaceChildren(h("p", { class: "err" }, e.message));
      }
    }, 30);
  }

  compute();
}

for (const el of document.querySelectorAll("[data-galois-explorer]")) mount(el);
