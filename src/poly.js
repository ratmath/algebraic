// Exact polynomials and certified real roots.
//
// A polynomial is an array of BigInt coefficients, constant term first,
// usually kept primitive (content 1, positive leading coefficient), since
// every question asked here (roots, square-free parts, Sturm counts) ignores
// constant factors. Rational points are [n, d] pairs.
//
// Real roots come back as AlgebraicReal: a square-free integer polynomial plus
// an isolating interval that can be shrunk on demand, and the exact rational
// value when the root is rational. That is the "sharpen as far as you like"
// form of an irrational number.

import { gcd, abs, reduce, toRational, pairToNumber, lcm } from "./q.js";
import { Rational } from "@ratmath/core";

export const degree = (p) => p.length - 1;

export function trim(p) {
  let n = p.length;
  while (n > 0 && p[n - 1] === 0n) n -= 1;
  return n === p.length ? p : p.slice(0, n);
}

export function primitive(p) {
  p = trim(p);
  if (!p.length) return p;
  let g = 0n;
  for (const c of p) g = gcd(g, c);
  const s = p[p.length - 1] < 0n ? -g : g;
  return s === 1n ? p : p.map((c) => c / s);
}

/** Rational coefficients ([n, d] pairs) to a primitive integer polynomial. */
export function fromPairs(coeffs) {
  let L = 1n;
  for (const [, d] of coeffs) L = lcm(L, d);
  return primitive(coeffs.map(([n, d]) => n * (L / d)));
}

export const derivative = (p) => p.slice(1).map((c, i) => c * BigInt(i + 1));

/** Sign of p(n/d): sign of Σ a_i n^i d^(deg−i), all integers. */
export function signAt(p, [n, d]) {
  let acc = 0n;
  let dp = 1n;
  // Horner on the homogenised form: acc = acc·n + a_i·d^(deg−i)
  for (let i = p.length - 1; i >= 0; i -= 1) {
    acc = acc * n + p[i] * dp;
    dp *= d;
  }
  return acc > 0n ? 1 : acc < 0n ? -1 : 0;
}

/** The numerator of p(n/d)·d^deg: p in homogeneous form. */
function homogeneous(p, [n, d]) {
  let acc = 0n;
  let dp = 1n;
  for (let i = p.length - 1; i >= 0; i -= 1) {
    acc = acc * n + p[i] * dp;
    dp *= d;
  }
  return acc;
}

/** √N for N = 2^(2^j). */
const sqrtPow2 = (N) => 1n << BigInt((N.toString(2).length - 1) >> 1);

/** p(x) as an exact pair. */
export function evalAt(p, [n, d]) {
  let acc = 0n;
  let dp = 1n;
  for (let i = p.length - 1; i >= 0; i -= 1) {
    acc = acc * n + p[i] * dp;
    dp *= d;
  }
  return reduce(acc, d ** BigInt(Math.max(0, p.length - 1)));
}

/** Pseudo-remainder of a by b, with the sign of the true remainder kept. */
function signedPrem(a, b) {
  a = [...a];
  const db = degree(b);
  const lb = b[db];
  let steps = 0;
  while (a.length - 1 >= db && a.length) {
    const da = a.length - 1;
    const la = a[da];
    for (let i = 0; i < a.length; i += 1) a[i] *= lb;
    for (let i = 0; i <= db; i += 1) a[da - db + i] -= la * b[i];
    a = trim(a);
    steps += 1;
  }
  // multiplied by lb^steps; undo the sign when that is negative
  if (lb < 0n && steps % 2 === 1) a = a.map((c) => -c);
  return a;
}

export function polyGcd(a, b) {
  a = primitive(a);
  b = primitive(b);
  while (b.length) {
    const r = signedPrem(a, b);
    a = b;
    b = primitive(r);
  }
  return primitive(a);
}

/** Exact quotient a / b up to a constant, as a primitive polynomial. */
export function divide(a, b) {
  a = [...a];
  const db = degree(b);
  const lb = b[db];
  const q = new Array(Math.max(0, a.length - db)).fill(0n);
  let scale = 1n;
  while (a.length - 1 >= db && a.length) {
    const da = a.length - 1;
    const la = a[da];
    for (let i = 0; i < a.length; i += 1) a[i] *= lb;
    for (let i = 0; i < q.length; i += 1) q[i] *= lb;
    scale *= lb;
    q[da - db] += la;
    for (let i = 0; i <= db; i += 1) a[da - db + i] -= la * b[i];
    a = trim(a);
  }
  if (a.length) throw new Error("polynomial division is not exact");
  return primitive(q);
}

/** Yun's square-free factorisation: [{ factor, multiplicity }], factors primitive, nonconstant. */
export function squareFree(p) {
  p = primitive(p);
  if (degree(p) < 1) return [];
  const out = [];
  const d = derivative(p);
  let a = polyGcd(p, d);
  let b = divide(p, a);
  // Musser: b holds each distinct factor once; peel off those of multiplicity k
  let k = 1;
  while (degree(b) > 0) {
    const y = polyGcd(b, a);
    const f = divide(b, y);
    if (degree(f) > 0) out.push({ factor: f, multiplicity: k });
    b = y;
    a = degree(y) > 0 ? divide(a, y) : a;
    k += 1;
  }
  return out;
}

/** Sturm chain of a square-free polynomial. */
function sturm(p) {
  const chain = [p, derivative(p)];
  for (;;) {
    const r = signedPrem(chain[chain.length - 2], chain[chain.length - 1]);
    if (!trim(r).length) break;
    // divide by the positive content only: the sign is what Sturm needs
    let g = 0n;
    for (const c of r) g = gcd(g, c);
    chain.push(trim(r).map((c) => -c / g));
  }
  return chain;
}

function variations(chain, x) {
  let v = 0;
  let last = 0;
  for (const q of chain) {
    const lead = q[q.length - 1] > 0n ? 1 : -1;
    const s = x === Infinity ? lead : x === -Infinity ? lead * (degree(q) % 2 ? -1 : 1) : signAt(q, x);
    if (s === 0) continue;
    if (last && s !== last) v += 1;
    last = s;
  }
  return v;
}

const half = ([a, b], [c, d]) => reduce(a * d + c * b, 2n * b * d);

/** Strict bound: every real root lies in (−B, B). */
function cauchyBound(p) {
  const lc = abs(p[p.length - 1]);
  let m = 0n;
  for (let i = 0; i < p.length - 1; i += 1) if (abs(p[i]) > m) m = abs(p[i]);
  return [m / lc + 2n, 1n];
}

/** Count of distinct real roots of a square-free polynomial. */
export function countRealRoots(p) {
  const chain = sturm(p);
  return variations(chain, -Infinity) - variations(chain, Infinity);
}

/** Isolate the real roots of a square-free polynomial, ascending. */
export function isolate(p) {
  if (degree(p) < 1) return [];
  if (degree(p) === 1) return [AlgebraicReal.rational(reduce(-p[0], p[1]))];
  const chain = sturm(p);
  const B = cauchyBound(p);
  const out = [];
  const stack = [[[-B[0], 1n], B]];
  while (stack.length) {
    const [a, b] = stack.pop();
    const n = variations(chain, a) - variations(chain, b);
    if (n === 0) continue;
    if (n === 1) {
      out.push(new AlgebraicReal(p, a, b));
      continue;
    }
    const m = half(a, b);
    if (signAt(p, m) === 0) {
      out.push(AlgebraicReal.rational(m));
      // nudge both halves off the root so endpoints stay non-roots
      const eps = [b[0] * a[1] - a[0] * b[1], b[1] * a[1] * 1024n];
      const left = reduce(m[0] * eps[1] - eps[0] * m[1], m[1] * eps[1]);
      const right = reduce(m[0] * eps[1] + eps[0] * m[1], m[1] * eps[1]);
      // shrink the gap until no other root sits inside it
      let l = left;
      let r = right;
      while (variations(chain, l) - variations(chain, r) !== 1 || signAt(p, l) === 0 || signAt(p, r) === 0) {
        l = half(l, m);
        r = half(m, r);
      }
      stack.push([a, l], [r, b]);
    } else {
      stack.push([a, m], [m, b]);
    }
  }
  return out.sort((x, y) => cmp(x.lo, y.lo));
}

const cmp = ([a, b], [c, d]) => {
  const x = a * d - c * b;
  return x > 0n ? 1 : x < 0n ? -1 : 0;
};

/** The fraction with the smallest denominator in [a, b] (a ≤ b). */
export function simplestBetween(a, b) {
  const floorQ = ([n, d]) => (n >= 0n ? n / d : -((-n + d - 1n) / d));
  const fl = floorQ(a);
  if (a[0] === fl * a[1]) return [fl, 1n];
  if (cmp([fl + 1n, 1n], b) <= 0) return [fl + 1n, 1n];
  const inv = ([n, d]) => reduce(d, n);
  const sub = ([n, d]) => reduce(n - fl * d, d);
  const [p, q] = simplestBetween(inv(sub(b)), inv(sub(a)));
  return reduce(fl * p + q, p);
}

/**
 * A real algebraic number: a root of the square-free integer polynomial
 * `poly` in the open interval (lo, hi), or an exact rational.
 */
export class AlgebraicReal {
  #rational = null;
  #checked = false;
  #N = 2n; // sub-intervals tried per refinement step

  constructor(poly, lo, hi) {
    this.poly = poly;
    this.lo = lo;
    this.hi = hi;
    this.sLo = poly ? signAt(poly, lo) : 0;
  }

  static rational(pair) {
    const r = new AlgebraicReal([-pair[0], pair[1]], pair, pair);
    r.#rational = pair;
    r.#checked = true;
    return r;
  }

  /**
   * Shrink the interval below width 2^-bits. Returns this.
   *
   * Quadratic interval refinement (Abbott): the secant through the ends
   * guesses which of N equal sub-intervals holds the root, and two exact sign
   * evaluations confirm or refute the guess. A confirmed guess squares N, so
   * the number of correct bits doubles per step once the root is close; a
   * refuted one falls back towards plain bisection (N = 2). Every accepted
   * interval is certified by its signs, exactly as in bisection.
   */
  refine(bits = 53) {
    if (this.#rational) return this;
    const target = [1n, 2n ** BigInt(bits)];
    const wide = () => cmp(reduce(this.hi[0] * this.lo[1] - this.lo[0] * this.hi[1], this.hi[1] * this.lo[1]), target) > 0;
    const exact = (m) => {
      this.#rational = m;
      this.#checked = true;
      this.lo = this.hi = m;
      return this;
    };
    let N = this.#N;
    while (wide()) {
      if (N > 2n) {
        // secant position t = f(lo)/(f(lo) − f(hi)) in [0, 1], on a common scale
        const deg = BigInt(this.poly.length - 1);
        const fl = homogeneous(this.poly, this.lo) * this.hi[1] ** deg;
        const fh = homogeneous(this.poly, this.hi) * this.lo[1] ** deg;
        if (fh === 0n) return exact(this.hi);
        if (fl !== 0n && fl !== fh) {
          let k = (N * fl) / (fl - fh);
          if (k < 0n) k = 0n;
          if (k > N - 1n) k = N - 1n;
          // grid points lo + (hi − lo)·k/N
          const at = (j) => reduce(this.lo[0] * this.hi[1] * (N - j) + this.hi[0] * this.lo[1] * j, this.lo[1] * this.hi[1] * N);
          const a = k === 0n ? this.lo : at(k);
          const c = k === N - 1n ? this.hi : at(k + 1n);
          const sa = k === 0n ? this.sLo : signAt(this.poly, a);
          if (sa === 0) return exact(a);
          const sc = k === N - 1n ? -this.sLo : signAt(this.poly, c);
          if (sc === 0) return exact(c);
          if (sa === this.sLo && sc !== this.sLo) {
            this.lo = a;
            this.hi = c;
            if (N < 1n << 4096n) N *= N;
            continue;
          }
          // the guess was wrong, but the two signs still tell which side the root is on
          if (sa !== this.sLo) this.hi = a;
          else this.lo = c;
          N = N > 4n ? sqrtPow2(N) : 2n;
          continue;
        }
      }
      const m = half(this.lo, this.hi);
      const sgn = signAt(this.poly, m);
      if (sgn === 0) return exact(m);
      if (sgn === this.sLo) this.lo = m;
      else this.hi = m;
      N = 4n;
    }
    this.#N = N;
    return this;
  }

  /** The exact value when rational, else null. Decided exactly (rational root bound). */
  get rational() {
    if (!this.#checked) {
      this.#checked = true;
      const L = abs(this.poly[this.poly.length - 1]);
      // rational roots p/q of a primitive polynomial have q | lc; two such
      // fractions differ by at least 1/lc², so once the interval is narrower
      // the simplest fraction in it is the only candidate
      this.refine(2 * L.toString(2).length + 2);
      if (!this.#rational) {
        const s = simplestBetween(this.lo, this.hi);
        if (s[1] <= L && signAt(this.poly, s) === 0) {
          this.#rational = s;
          this.lo = this.hi = s;
        }
      }
    }
    return this.#rational;
  }

  get isRational() {
    return this.rational !== null;
  }

  /** [lo, hi] as Core Rationals, width below 2^-bits. */
  interval(bits = 53) {
    this.refine(bits);
    return [toRational(this.lo), toRational(this.hi)];
  }

  /** The exact value as a Core Rational if rational, else null. */
  toRational() {
    const r = this.rational;
    return r && toRational(r);
  }

  toNumber() {
    if (this.#rational) return pairToNumber(...this.#rational);
    this.refine(60);
    return pairToNumber(...half(this.lo, this.hi));
  }

  /** Decimal digits, with the last digit certain up to one unit. */
  toString(digits = 12) {
    if (this.isRational) return toRational(this.#rational).toString();
    const bits = Math.ceil(digits * 3.33) + 4;
    const [lo] = this.interval(bits);
    return lo.toDecimal ? `${lo.toDecimal(digits)}…` : `${this.toNumber()}`;
  }
}

export { Rational };

/** Distinct real roots of the square-free p in the half-open interval (a, b] (pairs; ±Infinity allowed). */
export function countRootsBetween(p, a, b) {
  const chain = sturm(p);
  return variations(chain, a) - variations(chain, b);
}
