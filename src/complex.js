// Certified complex roots of integer polynomials, used for complex
// eigenvalues. Floats find the roots (Durand–Kerner); each is then certified
// by the Krawczyk test applied to p(x + iy) = 0 as two real equations, which
// proves the box holds exactly one root. Boxes for distinct roots are checked
// disjoint, so every root is accounted for exactly once.

import { derivative, degree } from "./poly.js";
import { reduce, addQ, subQ, mulQ, pairToNumber, toPair, toRational } from "./q.js";
import { quadraticParts, formatComplex, certifiedDecimal } from "./radical.js";
import { iv, ivAdd, ivSub, ivMul, ivWidth, krawczyk, toCore } from "./interval.js";

/** Float roots of an integer polynomial (constant first) by Durand–Kerner. */
export function floatRoots(p) {
  const n = degree(p);
  if (n < 1) return [];
  const lc = Number(p[n]);
  const c = p.map((x) => Number(x) / lc);
  let R = 1;
  for (let i = 0; i < n; i += 1) R = Math.max(R, 1 + Math.abs(c[i]));
  let z = Array.from({ length: n }, (_, k) => ({ re: R * 0.9 * Math.cos(0.4 + (2 * Math.PI * k) / n), im: R * 0.9 * Math.sin(0.4 + (2 * Math.PI * k) / n) }));
  const evalAt = (w) => {
    let re = 1;
    let im = 0;
    for (let i = n - 1; i >= 0; i -= 1) {
      [re, im] = [re * w.re - im * w.im + c[i], re * w.im + im * w.re];
    }
    return { re, im };
  };
  for (let it = 0; it < 2000; it += 1) {
    let move = 0;
    const next = z.map((zi, i) => {
      let dr = 1;
      let di = 0;
      z.forEach((zj, j) => {
        if (i === j) return;
        const ar = zi.re - zj.re;
        const ai = zi.im - zj.im;
        [dr, di] = [dr * ar - di * ai, dr * ai + di * ar];
      });
      const v = evalAt(zi);
      const den = dr * dr + di * di || 1e-300;
      const qr = (v.re * dr + v.im * di) / den;
      const qi = (v.im * dr - v.re * di) / den;
      move = Math.max(move, Math.hypot(qr, qi) / Math.max(1, Math.hypot(zi.re, zi.im)));
      return { re: zi.re - qr, im: zi.im - qi };
    });
    z = next;
    if (move < 1e-15) break;
  }
  return z;
}

// complex values: [re, im] of pairs; complex intervals: [reIv, imIv]
function hornerPoint(p, [a, b]) {
  let re = [0n, 1n];
  let im = [0n, 1n];
  for (let i = p.length - 1; i >= 0; i -= 1) {
    [re, im] = [addQ(subQ(mulQ(re, a), mulQ(im, b)), [p[i], 1n]), addQ(mulQ(re, b), mulQ(im, a))];
  }
  return [re, im];
}

export function hornerBox(p, [X, Y]) {
  let re = [[0n, 1n], [0n, 1n]];
  let im = [[0n, 1n], [0n, 1n]];
  for (let i = p.length - 1; i >= 0; i -= 1) {
    const c = [[p[i], 1n], [p[i], 1n]];
    [re, im] = [ivAdd(ivSub(ivMul(re, X), ivMul(im, Y)), c), ivAdd(ivMul(re, Y), ivMul(im, X))];
  }
  return [re, im];
}

/**
 * A complex algebraic number: the unique root of the square-free integer
 * polynomial `poly` in the box re × im (each [lo, hi] of pairs).
 */
export class ComplexAlgebraic {
  constructor(poly, re, im) {
    this.poly = poly;
    this.re = re;
    this.im = im;
  }

  conjugate() {
    const [a, b] = this.im;
    return new ComplexAlgebraic(this.poly, this.re, [[-b[0], b[1]], [-a[0], a[1]]]);
  }

  /** Shrink the box below 2^-bits per side (Krawczyk contraction). */
  refine(bits = 53) {
    const wide = (x) => pairToNumber(...ivWidth(x)) > 2 ** -bits;
    if (!wide(this.re) && !wide(this.im)) return this;
    const { f, J } = systemFor(this.poly);
    const it = krawczyk(f, J, [this.re, this.im], { bits });
    let r;
    while (!(r = it.next()).done);
    if (r.value.status === "unique") [this.re, this.im] = r.value.box;
    return this;
  }

  /** { re: [lo, hi], im: [lo, hi] } as Core Rationals. */
  box(bits = 53) {
    this.refine(bits);
    return { re: this.re.map(toRational), im: this.im.map(toRational) };
  }

  toNumbers() {
    this.refine(60);
    const mid = ([a, b]) => (pairToNumber(...a) + pairToNumber(...b)) / 2;
    return { re: mid(this.re), im: mid(this.im) };
  }

  /**
   * Exact "a + bi" form when the number is a root of a quadratic (so every
   * Gaussian rational and every (p + q√−m)/r), else certified digits for
   * each part, the last digit off by at most one unit, ending in "…".
   */
  toString(digits = 10) {
    if (degree(this.poly) === 2) return formatComplex(quadraticParts(this.poly, this.#imSign() * (this.poly[2] > 0n ? 1 : -1)));
    return this.toDecimal(digits);
  }

  /** Certified digits for each part, "a + bi…", the last digit off by at most one unit. */
  toDecimal(digits = 10) {
    const reS = certifiedDecimal(() => this.re, (b) => this.refine(b), digits);
    const imS = certifiedDecimal(() => this.im, (b) => this.refine(b), digits);
    const neg = imS.startsWith("−");
    return `${reS} ${neg ? "−" : "+"} ${neg ? imS.slice(1) : imS}i…`;
  }

  /** Sign of the imaginary part (never 0 for a non-real number). */
  #imSign() {
    for (let bits = 8; bits <= 8192; bits *= 2) {
      if (this.im[0][0] > 0n) return 1;
      if (this.im[1][0] < 0n) return -1;
      this.refine(bits);
    }
    throw new Error("could not settle the sign of the imaginary part");
  }

  get imSign() {
    return this.#imSign();
  }
}

function systemFor(p) {
  const dp = derivative(p);
  const f = ([x, y]) => {
    const [re, im] = hornerPoint(p, [toPair(x), toPair(y)]);
    return [re, im];
  };
  const J = ([X, Y]) => {
    const [a, b] = hornerBox(dp, [iv(X), iv(Y)]);
    const negB = [[-b[1][0], b[1][1]], [-b[0][0], b[0][1]]];
    return [[a, negB], [b, a]].map((row) => row.map(toCore));
  };
  return { f, J };
}

const disjoint = (u, v) => {
  const lt = ([a, b], [c, d]) => a * d < c * b;
  const apart = (x, y) => lt(x[1], y[0]) || lt(y[1], x[0]);
  return apart(u.re, v.re) || apart(u.im, v.im);
};

/**
 * Certified non-real roots of a square-free integer polynomial, given how
 * many there are (from a Sturm count). Returns ComplexAlgebraic[] (upper half
 * plane first, each followed by its conjugate), or null if certification
 * failed (the caller then reports float approximations).
 */
export function certifyComplexRoots(p, count) {
  if (count === 0) return [];
  const approx = floatRoots(p).filter((z) => z.im > 1e-9 * Math.max(1, Math.abs(z.re))).sort((a, b) => a.re - b.re || a.im - b.im);
  if (approx.length * 2 !== count) return null;
  const { f, J } = systemFor(p);
  const out = [];
  for (const z of approx) {
    let radius = 1e-10 * Math.max(1, Math.hypot(z.re, z.im));
    let found = null;
    for (let attempt = 0; attempt < 10 && !found; attempt += 1) {
      const r = toPair(radius);
      const c = [toPair(z.re), toPair(z.im)];
      const X = c.map((m) => [subQ(m, r), addQ(m, r)]);
      const it = krawczyk(f, J, X, { bits: 40 });
      let res;
      while (!(res = it.next()).done);
      if (res.value.status === "unique") found = res.value.box;
      radius *= 10;
    }
    if (!found) return null;
    const root = new ComplexAlgebraic(p, found[0], found[1]);
    if (out.some((o) => !disjoint(o, root))) return null;
    out.push(root);
  }
  return out.flatMap((r) => [r, r.conjugate()]);
}

export { reduce };
