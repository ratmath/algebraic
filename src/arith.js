// Exact arithmetic on algebraic numbers, real (AlgebraicReal) and complex
// (ComplexAlgebraic).
//
// a ± b and a·b are roots of a resultant: the polynomial whose roots are all
// sums (or products) of a conjugate of a with a conjugate of b. That
// polynomial is factored into irreducibles, so every result carries its
// minimal polynomial and degrees stay as small as they can be. The right root
// is then picked by refining enclosures until exactly one candidate fits.

import { AlgebraicReal, isolate, primitive, degree, signAt, squareFree } from "./poly.js";
import { ComplexAlgebraic, certifyComplexRoots } from "./complex.js";
import { factor } from "./factor.js";
import { reduce, addQ, subQ, mulQ, divQ, toPair, sqrtFloor, detInt, abs } from "./q.js";
import { ivAdd, ivSub, ivMul } from "./interval.js";

export const MAX_DEGREE = 64;

const isReal = (a) => a instanceof AlgebraicReal;
const isComplex = (a) => a instanceof ComplexAlgebraic;
const cmp = ([a, b], [c, d]) => {
  const x = a * d - c * b;
  return x > 0n ? 1 : x < 0n ? -1 : 0;
};
const Q = (n, d = 1n) => reduce(BigInt(n), BigInt(d));
const ZERO = [0n, 1n];

/** A rational (pair [n, d], Core Rational, bigint, number or string) as an AlgebraicReal. */
export function fromRational(x) {
  const [n, d] = Array.isArray(x) ? x : toPair(x);
  return AlgebraicReal.rational(reduce(BigInt(n), BigInt(d)));
}

/** The exact rational value of a, or null. */
export const rationalOf = (a) => (isReal(a) ? a.rational : null);

// ---------------------------------------------------------------------------
// Roots of integer polynomials

function trimInt(p) {
  p = p.map(BigInt);
  while (p.length > 1 && p[p.length - 1] === 0n) p.pop();
  return p;
}

/** Every root of an irreducible integer polynomial: real ones ascending, then complex pairs. */
function rootsOfIrreducible(f) {
  const real = isolate(f);
  const nComplex = degree(f) - real.length;
  if (nComplex === 0) return real;
  const complex = certifyComplexRoots(f, nComplex);
  if (!complex) throw new Error(`could not certify the complex roots of a degree ${degree(f)} polynomial`);
  return [...real, ...complex];
}

/**
 * All roots of an integer or rational polynomial (coefficients constant first,
 * as bigints, pairs or Rationals), each with its multiplicity. Real roots are
 * AlgebraicReal (ascending); with { complex: true } the non-real ones follow as
 * ComplexAlgebraic. Every root carries its minimal polynomial.
 */
export function roots(coeffs, { complex = false } = {}) {
  const pairs = coeffs.map((c) => (Array.isArray(c) ? c : toPair(c)));
  let L = 1n;
  const g = (x, y) => (y ? g(y, x % y) : x);
  for (const [, d] of pairs) L = (L * d) / g(L, d);
  const ints = trimInt(pairs.map(([n, d]) => n * (L / d)));
  if (ints.length < 2) return [];
  const out = [];
  for (const { factor: f, multiplicity } of factor(ints).factors) {
    const rs = complex ? rootsOfIrreducible(f) : isolate(f);
    for (const value of rs) out.push({ value, multiplicity });
  }
  const real = out.filter((r) => isReal(r.value)).sort((x, y) => compare(x.value, y.value));
  return [...real, ...out.filter((r) => !isReal(r.value))];
}

/** Real roots only (shorthand for roots(p)). */
export const realRoots = (coeffs) => roots(coeffs);

/** The minimal polynomial of a (primitive integer, constant first). */
export function minimalPolynomial(a) {
  const r = a.rational;
  if (r) return primitive([-r[0], r[1]]);
  const { factors } = factor(a.poly);
  if (factors.length === 1) return factors[0].factor;
  return pick(a.poly, (bits) => box(a, bits)).poly;
}

/** The same number, re-expressed on its minimal polynomial. */
export function simplify(a) {
  if (a.rational) return a;
  const { factors } = factor(a.poly);
  if (factors.length === 1) return a;
  return pick(a.poly, (bits) => box(a, bits));
}

// ---------------------------------------------------------------------------
// Enclosures: boxes { re: [lo, hi], im: [lo, hi] } of pairs

function box(a, bits) {
  if (isReal(a)) {
    const r = a.rational;
    if (r) return { re: [r, r], im: [ZERO, ZERO] };
    a.refine(bits);
    return { re: [a.lo, a.hi], im: [ZERO, ZERO] };
  }
  a.refine(bits);
  return { re: a.re, im: a.im };
}
const meets = ([a, b], [c, d]) => cmp(b, c) >= 0 && cmp(d, a) >= 0;
const overlaps = (u, v) => meets(u.re, v.re) && meets(u.im, v.im);

function pick(poly, enclose) {
  const { factors } = factor(trimInt(poly));
  const candidates = [];
  for (const { factor: f } of factors) {
    if (degree(f) > MAX_DEGREE) throw new RangeError(`the result needs a polynomial of degree ${degree(f)}; the limit is ${MAX_DEGREE}`);
    candidates.push(...rootsOfIrreducible(f));
  }
  for (let bits = 8; bits <= 8192; bits *= 2) {
    const E = enclose(bits);
    const hits = candidates.filter((c) => overlaps(box(c, bits), E));
    if (hits.length === 1) return hits[0];
    if (hits.length === 0) throw new Error("algebraic: the enclosure missed every candidate root");
  }
  throw new Error("algebraic: could not separate the candidate roots");
}

// ---------------------------------------------------------------------------
// Resultants by evaluation and interpolation

function resultant(p, q) {
  const m = p.length - 1;
  const n = q.length - 1;
  const size = m + n;
  const hp = [...p].reverse();
  const hq = [...q].reverse();
  const rows = [];
  for (let i = 0; i < n; i += 1) rows.push(Array.from({ length: size }, (_, j) => (j - i >= 0 && j - i <= m ? hp[j - i] : 0n)));
  for (let i = 0; i < m; i += 1) rows.push(Array.from({ length: size }, (_, j) => (j - i >= 0 && j - i <= n ? hq[j - i] : 0n)));
  return detInt(rows);
}

const binom = (n, k) => {
  let r = 1n;
  for (let i = 0n; i < BigInt(k); i += 1n) r = (r * (BigInt(n) - i)) / (i + 1n);
  return r;
};
/** Coefficients in y of q(t − y). */
function shifted(q, t) {
  const out = Array(q.length).fill(0n);
  q.forEach((c, k) => {
    for (let j = 0; j <= k; j += 1) out[j] += c * binom(k, j) * t ** BigInt(k - j) * (j % 2 ? -1n : 1n);
  });
  return out;
}
/** Coefficients in y of yⁿ·q(t / y). */
function scaled(q, t) {
  const n = q.length - 1;
  return q.map((_, j) => q[n - j] * t ** BigInt(n - j));
}

function interpolate(D, f) {
  const xs = Array.from({ length: D + 1 }, (_, t) => BigInt(t));
  const coef = xs.map((t) => Q(f(t)));
  for (let j = 1; j <= D; j += 1) for (let i = D; i >= j; i -= 1) coef[i] = divQ(subQ(coef[i], coef[i - 1]), Q(xs[i] - xs[i - j]));
  let poly = [coef[D]];
  for (let i = D - 1; i >= 0; i -= 1) {
    const next = Array(poly.length + 1).fill(ZERO);
    poly.forEach((c, k) => {
      next[k + 1] = addQ(next[k + 1], c);
      next[k] = subQ(next[k], mulQ(c, Q(xs[i])));
    });
    next[0] = addQ(next[0], coef[i]);
    poly = next;
  }
  let L = 1n;
  const g = (a, b) => (b ? g(b, a % b) : a);
  for (const [, d] of poly) L = (L * d) / g(L, d);
  return primitive(trimInt(poly.map(([n, d]) => n * (L / d))));
}

/** The polynomial to compute with: the minimal one when cheap, without a stray root at 0. */
function usable(a) {
  const r = a.rational;
  if (r) return [-r[0], r[1]];
  let p = a.poly;
  while (p.length > 1 && p[0] === 0n) p = p.slice(1);
  return p;
}

// ---------------------------------------------------------------------------
// Operations

export function neg(a) {
  if (isReal(a)) {
    const r = a.rational;
    if (r) return fromRational([-r[0], r[1]]);
    return new AlgebraicReal(primitive(a.poly.map((c, i) => (i % 2 ? -c : c))), [-a.hi[0], a.hi[1]], [-a.lo[0], a.lo[1]]);
  }
  const p = primitive(a.poly.map((c, i) => (i % 2 ? -c : c)));
  const flip = ([lo, hi]) => [[-hi[0], hi[1]], [-lo[0], lo[1]]];
  return new ComplexAlgebraic(p, flip(a.re), flip(a.im));
}

export function add(a, b) {
  const ra = a.rational;
  const rb = b.rational;
  if (ra && rb) return fromRational(addQ(ra, rb));
  const p = usable(a);
  const q = usable(b);
  const R = interpolate((p.length - 1) * (q.length - 1), (t) => resultant(p, shifted(q, t)));
  return pick(R, (bits) => {
    const A = box(a, bits);
    const B = box(b, bits);
    return { re: ivAdd(A.re, B.re), im: ivAdd(A.im, B.im) };
  });
}

export const sub = (a, b) => add(a, neg(b));

export function mul(a, b) {
  const ra = a.rational;
  const rb = b.rational;
  if (ra && rb) return fromRational(mulQ(ra, rb));
  if ((ra && ra[0] === 0n) || (rb && rb[0] === 0n)) return fromRational([0n, 1n]);
  const p = usable(a);
  const q = usable(b);
  const R = interpolate((p.length - 1) * (q.length - 1), (t) => resultant(p, scaled(q, t)));
  return pick(R, (bits) => {
    const A = box(a, bits);
    const B = box(b, bits);
    return { re: ivSub(ivMul(A.re, B.re), ivMul(A.im, B.im)), im: ivAdd(ivMul(A.re, B.im), ivMul(A.im, B.re)) };
  });
}

/** -1, 0 or 1 for a real a, exactly. */
export function sign(a) {
  if (!isReal(a)) throw new TypeError("sign needs a real number");
  const r = a.rational;
  if (r) return r[0] > 0n ? 1 : r[0] < 0n ? -1 : 0;
  for (let bits = 8; ; bits *= 2) {
    const { re: [lo, hi] } = box(a, bits);
    if (lo[0] > 0n) return 1;
    if (hi[0] < 0n) return -1;
  }
}

export const isZero = (a) => {
  const r = a.rational;
  return r ? r[0] === 0n : false; // an irrational or non-real number is never 0
};

function ivInv([lo, hi]) {
  if (lo[0] <= 0n && hi[0] >= 0n) return null;
  return [reduce(hi[1], hi[0]), reduce(lo[1], lo[0])].sort(cmp);
}

export function inv(a) {
  if (isZero(a)) throw new RangeError("division by zero");
  const r = a.rational;
  if (r) return fromRational([r[1], r[0]]);
  const rev = primitive([...usable(a)].reverse());
  if (isReal(a)) {
    sign(a);
    const [lo, hi] = [reduce(a.hi[1], a.hi[0]), reduce(a.lo[1], a.lo[0])].sort(cmp);
    return new AlgebraicReal(rev[rev.length - 1] < 0n ? rev.map((c) => -c) : rev, lo, hi);
  }
  return pick(rev, (bits) => {
    for (let b = bits; ; b *= 2) {
      const A = box(a, b);
      const d = ivInv(ivAdd(ivMul(A.re, A.re), ivMul(A.im, A.im)));
      if (!d) continue;
      return { re: ivMul(A.re, d), im: ivMul(ivSub([ZERO, ZERO], A.im), d) };
    }
  });
}

export const div = (a, b) => mul(a, inv(b));

export function pow(a, k) {
  if (k < 0) return pow(inv(a), -k);
  let out = fromRational([1n, 1n]);
  let base = a;
  while (k > 0) {
    if (k & 1) out = mul(out, base);
    k >>= 1;
    if (k) base = mul(base, base);
  }
  return out;
}

/** -1, 0, 1 as real a compares to real b. */
export const compare = (a, b) => sign(sub(a, b));

/** Equality, exactly (works for complex numbers too). */
export function equals(a, b) {
  if (isReal(a) !== isReal(b)) return false;
  const ra = a.rational;
  const rb = b.rational;
  if (ra || rb) return !!(ra && rb && ra[0] === rb[0] && ra[1] === rb[1]);
  return isZero(sub(a, b));
}

let I_CACHE = null;
/** The imaginary unit i, the root of x² + 1 with positive imaginary part. */
export function I() {
  if (!I_CACHE) I_CACHE = rootsOfIrreducible([1n, 0n, 1n]).find((z) => isComplex(z) && z.im[0][0] > 0n);
  return I_CACHE;
}

export const conjugate = (a) => (isReal(a) ? a : a.conjugate());

/** Real part, as an exact real algebraic number. */
export function re(a) {
  if (isReal(a)) return a;
  const s = add(a, a.conjugate());
  return mul(s, fromRational([1n, 2n]));
}
/** Imaginary part, as an exact real algebraic number. */
export function im(a) {
  if (isReal(a)) return fromRational([0n, 1n]);
  const d = sub(a, a.conjugate()); // 2i·Im
  return mul(div(d, I()), fromRational([1n, 2n]));
}

// ---------------------------------------------------------------------------
// Roots

/** ⌊x^(1/n)⌋ bounds for a nonnegative rational, about `bits` bits. */
function nthBounds([num, den], n, bits) {
  const s = 2n ** BigInt(bits);
  const N = BigInt(n);
  // integer nth root of num·den^(n-1)·s^n, then divide by den·s
  const target = num * den ** (N - 1n) * s ** N;
  let lo = 0n;
  let hi = 1n;
  while (hi ** N <= target) hi *= 2n;
  while (hi - lo > 1n) {
    const mid = (lo + hi) / 2n;
    if (mid ** N <= target) lo = mid;
    else hi = mid;
  }
  return [reduce(lo, den * s), reduce(lo + 1n, den * s)];
}

/** The real n-th root (n ≥ 2) of a real a: principal for a ≥ 0; for a < 0 only odd n. */
export function nthRoot(a, n) {
  if (!Number.isInteger(n) || n < 2) throw new RangeError("nthRoot needs n ≥ 2");
  const sg = sign(a);
  if (sg === 0) return a;
  if (sg < 0) {
    if (n % 2 === 0) throw new RangeError("an even root of a negative number is not real; use sqrt for i·√");
    return neg(nthRoot(neg(a), n));
  }
  const r = a.rational;
  if (r) {
    const N = BigInt(n);
    const root = (v) => {
      let lo = 0n;
      let hi = 1n;
      while (hi ** N <= v) hi *= 2n;
      while (hi - lo > 1n) {
        const mid = (lo + hi) / 2n;
        if (mid ** N <= v) lo = mid;
        else hi = mid;
      }
      return lo;
    };
    const a1 = root(r[0]);
    const b1 = root(r[1]);
    if (a1 ** N === r[0] && b1 ** N === r[1]) return fromRational([a1, b1]);
  }
  const p = usable(a);
  const pn = [];
  p.forEach((c, i) => {
    pn[n * i] = c;
    if (i < p.length - 1) for (let k = 1; k < n; k += 1) pn[n * i + k] = 0n;
  });
  return pick(primitive(pn), (bits) => {
    const { re: [lo, hi] } = box(a, bits + 4);
    return { re: [nthBounds(lo, n, bits)[0], nthBounds(hi, n, bits)[1]], im: [ZERO, ZERO] };
  });
}

/** √a for real a: real when a ≥ 0, i·√(−a) when a < 0. */
export function sqrt(a) {
  if (!isReal(a)) throw new RangeError("square roots of non-real numbers are not supported yet");
  return sign(a) >= 0 ? nthRoot(a, 2) : mul(I(), nthRoot(neg(a), 2));
}

export { isReal, isComplex, signAt, squareFree, abs, sqrtFloor };
