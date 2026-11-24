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
import { reduce, addQ, subQ, mulQ, divQ, toPair, sqrtFloor } from "./q.js";
import { quadraticParts, realSide, formatReal, fmtQ, joinParts, nestedRadical, pureRootParts, formatPureRoot } from "./radical.js";
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
  // distinct roots, so refining their intervals always separates them (no resultants needed)
  const apart = (a, b) => {
    for (let bits = 8; ; bits *= 2) {
      const { re: A } = box(a, bits);
      const { re: B } = box(b, bits);
      if (cmp(A[1], B[0]) < 0) return -1;
      if (cmp(B[1], A[0]) < 0) return 1;
    }
  };
  const real = out.filter((r) => isReal(r.value)).sort((x, y) => apart(x.value, y.value));
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

function pick(poly, enclose, realOnly = false) {
  const { factors } = factor(trimInt(poly));
  const candidates = [];
  for (const { factor: f } of factors) {
    if (degree(f) > MAX_DEGREE) throw new RangeError(`the result needs a polynomial of degree ${degree(f)}; the limit is ${MAX_DEGREE}`);
    // a real result is a real root: the complex roots need not be certified at all
    candidates.push(...(realOnly ? isolate(f) : rootsOfIrreducible(f)));
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
// Composed sums and products by power sums
//
// The polynomial whose roots are all a_i + b_j (or a_i·b_j) is the resultant
// of p(x) and q(t − x) (or xⁿq(t/x)). Its power sums come straight from those
// of p and q:  Σ (a_i·b_j)^k = s_k(p)·s_k(q)  and
// Σ (a_i + b_j)^k = Σ_j C(k, j)·s_j(p)·s_(k−j)(q),  and Newton's identities
// turn power sums back into coefficients. That is O(D²) rational operations
// for a result of degree D, where evaluating and interpolating Sylvester
// determinants costs D of them at O(D³) each.

/** Power sums s_0 … s_K of the roots of p (pairs). */
function powerSums(p, K) {
  const n = p.length - 1;
  const lc = Q(p[n]);
  // elementary symmetric e_i = (−1)^i a_(n−i)/a_n
  const e = (i) => (i > n ? ZERO : divQ(Q(i % 2 ? -p[n - i] : p[n - i]), lc));
  const s = [Q(n)];
  for (let k = 1; k <= K; k += 1) {
    let acc = mulQ(Q(k % 2 ? k : -k), e(k));
    for (let i = 1; i < k && i <= n; i += 1) {
      const term = mulQ(e(i), s[k - i]);
      acc = i % 2 ? addQ(acc, term) : subQ(acc, term);
    }
    s.push(acc);
  }
  return s;
}

/** The primitive integer polynomial of degree D with power sums s_1 … s_D. */
function fromPowerSums(s, D) {
  const e = [Q(1)];
  for (let k = 1; k <= D; k += 1) {
    let acc = ZERO;
    for (let i = 1; i <= k; i += 1) {
      const term = mulQ(e[k - i], s[i]);
      acc = i % 2 ? addQ(acc, term) : subQ(acc, term);
    }
    e.push(divQ(acc, Q(k)));
  }
  const g = (a, b) => (b ? g(b, a % b) : a);
  let L = 1n;
  for (const [, d] of e) L = (L * d) / g(L, d);
  // x^D − e₁x^(D−1) + e₂x^(D−2) − …, constant first
  const out = [];
  for (let k = D; k >= 0; k -= 1) out.push((k % 2 ? -e[k][0] : e[k][0]) * (L / e[k][1]));
  return primitive(trimInt(out));
}

function composedSum(p, q) {
  const D = (p.length - 1) * (q.length - 1);
  const sp = powerSums(p, D);
  const sq = powerSums(q, D);
  const s = [Q(D)];
  for (let k = 1; k <= D; k += 1) {
    let acc = ZERO;
    let c = 1n; // C(k, j)
    for (let j = 0; j <= k; j += 1) {
      acc = addQ(acc, mulQ(Q(c), mulQ(sp[j], sq[k - j])));
      c = (c * BigInt(k - j)) / BigInt(j + 1);
    }
    s.push(acc);
  }
  return fromPowerSums(s, D);
}

function composedProduct(p, q) {
  const D = (p.length - 1) * (q.length - 1);
  const sp = powerSums(p, D);
  const sq = powerSums(q, D);
  return fromPowerSums(sp.map((x, k) => mulQ(x, sq[k])), D);
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
  return pick(composedSum(p, q), (bits) => {
    const A = box(a, bits);
    const B = box(b, bits);
    return { re: ivAdd(A.re, B.re), im: ivAdd(A.im, B.im) };
  }, isReal(a) && isReal(b));
}

export const sub = (a, b) => add(a, neg(b));

export function mul(a, b) {
  const ra = a.rational;
  const rb = b.rational;
  if (ra && rb) return fromRational(mulQ(ra, rb));
  if ((ra && ra[0] === 0n) || (rb && rb[0] === 0n)) return fromRational([0n, 1n]);
  const p = usable(a);
  const q = usable(b);
  return pick(composedProduct(p, q), (bits) => {
    const A = box(a, bits);
    const B = box(b, bits);
    return { re: ivSub(ivMul(A.re, B.re), ivMul(A.im, B.im)), im: ivAdd(ivMul(A.re, B.im), ivMul(A.im, B.re)) };
  }, isReal(a) && isReal(b));
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
  if (degree(a.poly) === 2) return fromRational(reduce(-a.poly[1], 2n * a.poly[2]));
  const s = add(a, a.conjugate());
  return mul(s, fromRational([1n, 2n]));
}
/** Imaginary part, as an exact real algebraic number. */
export function im(a) {
  if (isReal(a)) return fromRational([0n, 1n]);
  if (degree(a.poly) === 2) {
    // (b² − 4ac) < 0: Im = ±√(4ac − b²) / 2a, sign from the box
    const [c, b, lc] = a.poly;
    const v = nthRoot(fromRational(reduce(4n * lc * c - b * b, 4n * lc * lc)), 2);
    return a.imSign > 0 ? v : neg(v);
  }
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

/** Box product of two enclosures. */
function boxMul(A, B) {
  return { re: ivSub(ivMul(A.re, B.re), ivMul(A.im, B.im)), im: ivAdd(ivMul(A.re, B.im), ivMul(A.im, B.re)) };
}
function boxPow(A, n) {
  let out = A;
  for (let k = 1; k < n; k += 1) out = boxMul(out, A);
  return out;
}
const pairLt = (x, y) => cmp(x, y) < 0;

/**
 * The principal n-th root of a non-real z, or of a negative real when n is
 * even: the n-th root with the largest real part (ties, which happen only for
 * negative reals, go to the one with positive imaginary part).
 *
 * Candidates are the roots of m(xⁿ), m the minimal polynomial of z. Each one's
 * n-th power is a conjugate of z; interval powers decide which, and then the
 * principal one among the n roots of z itself is picked by real part.
 */
function principalRoot(z, n) {
  const m = minimalPolynomial(z);
  const conj = rootsOfIrreducible(m);
  const pn = [];
  m.forEach((c, i) => {
    pn[n * i] = c;
    if (i < m.length - 1) for (let k = 1; k < n; k += 1) pn[n * i + k] = 0n;
  });
  const candidates = [];
  for (const { factor: f } of factor(primitive(pn)).factors) {
    if (degree(f) > MAX_DEGREE) throw new RangeError(`the result needs a polynomial of degree ${degree(f)}; the limit is ${MAX_DEGREE}`);
    candidates.push(...rootsOfIrreducible(f));
  }
  // which conjugate of z is z itself, and which conjugate each candidate's power is
  const owner = (encl, bits) => {
    const hits = conj.filter((c) => overlaps(box(c, bits), encl));
    return hits.length === 1 ? hits[0] : null;
  };
  let mine = null;
  for (let bits = 8; bits <= 8192 && !mine; bits *= 2) {
    const self = owner(box(z, bits), bits);
    if (!self) continue;
    const settled = candidates.map((c) => owner(boxPow(box(c, bits), n), bits));
    if (settled.every(Boolean)) mine = candidates.filter((c, k) => settled[k] === self);
  }
  if (!mine || mine.length !== n) throw new Error("algebraic: could not sort the n-th roots");
  const negReal = isReal(z);
  for (let bits = 8; bits <= 8192; bits *= 2) {
    let pool = mine;
    if (negReal) {
      const up = pool.filter((c) => !isReal(c) && box(c, bits).im[0][0] > 0n);
      if (up.length + pool.filter((c) => isReal(c) || box(c, bits).im[1][0] < 0n).length !== pool.length) continue;
      pool = up;
    }
    const boxes = pool.map((c) => box(c, bits));
    const best = boxes.findIndex((B, k) => boxes.every((C, j) => j === k || pairLt(C.re[1], B.re[0])));
    if (best >= 0) return pool[best];
  }
  throw new Error("algebraic: could not pick the principal root");
}

/**
 * The n-th root (n ≥ 2) of a. For real a ≥ 0 the positive real root; for real
 * a < 0 the real root when n is odd and the principal complex root when n is
 * even; for non-real a the principal root (largest real part).
 */
export function nthRoot(a, n) {
  if (!Number.isInteger(n) || n < 2) throw new RangeError("nthRoot needs n ≥ 2");
  if (!isReal(a)) return principalRoot(a, n);
  const sg = sign(a);
  if (sg === 0) return a;
  if (sg < 0) {
    if (n % 2 === 0) return n === 2 ? mul(I(), nthRoot(neg(a), 2)) : principalRoot(a, n);
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

/** The principal square root: √a ≥ 0 for real a ≥ 0, i·√(−a) for a < 0, the root with positive real part otherwise. */
export const sqrt = (a) => nthRoot(a, 2);

/** The complex number x + y·i from real parts (rationals or AlgebraicReals). */
export function complex(x, y = 0) {
  const X = isReal(x) || isComplex(x) ? x : fromRational(x);
  const Y = isReal(y) || isComplex(y) ? y : fromRational(y);
  if (!isReal(X) || !isReal(Y)) throw new TypeError("complex(x, y) needs real parts");
  return isZero(Y) ? X : add(X, mul(Y, I()));
}

/** |a|², exactly: a·ā (a real algebraic number). */
export function norm(a) {
  if (isReal(a)) return mul(a, a);
  if (degree(a.poly) === 2) {
    // a·ā is the product of the two roots: c/a
    const [c, , lc] = a.poly;
    return fromRational(reduce(c, lc));
  }
  return mul(a, a.conjugate());
}

/** |a|, exactly: the absolute value of a real, the modulus of a complex number. */
export function abs(a) {
  if (isReal(a)) return sign(a) < 0 ? neg(a) : a;
  return nthRoot(norm(a), 2);
}

/**
 * A readable exact form when one is cheap: rationals as fractions, quadratic
 * irrationals as (p + q√d)/r, quadratic complex numbers as a + bi; otherwise
 * certified digits. `{ digits }` sets how many.
 */
export function format(a, { digits = 12 } = {}) {
  if (isReal(a)) {
    const r = a.rational;
    if (r) return fmtQ(r);
    const s = degree(a.poly) === 2 ? a : simplify(a);
    if (degree(s.poly) === 2) {
      const side = realSide(s.poly, () => [s.lo, s.hi], (b) => s.refine(b));
      return formatReal(quadraticParts(s.poly, side));
    }
    if (degree(s.poly) === 4 && !pureRootParts(s.poly)) {
      // a nested radical ±√(A + B√d)/C when the square is quadratic, such as cos(π/8) = √(2 + √2)/2
      const sq = mul(s, s);
      if (degree(sq.poly) === 2) {
        const side = realSide(sq.poly, () => [sq.lo, sq.hi], (b) => sq.refine(b));
        return (sign(s) < 0 ? "−" : "") + nestedRadical(quadraticParts(sq.poly, side));
      }
    }
    // p + q·ⁿ√m when the minimal polynomial is a shifted pure power, such as ∛2 or 1 − ⁴√5/2
    const pure = pureRootParts(s.poly);
    if (pure) {
      let side = 1;
      if (pure.n % 2 === 0) {
        for (let bits = 8; ; bits *= 2) {
          s.refine(bits);
          if (cmp(s.lo, pure.p) > 0) break;
          if (cmp(s.hi, pure.p) < 0) { side = -1; break; }
        }
      }
      return formatPureRoot(pure, side);
    }
    return s.toString(digits);
  }
  if (degree(a.poly) === 2) return a.toString(digits);
  if (degree(a.poly) <= 4) {
    // parts with exact forms, such as ∛i = √3/2 + (1/2)i
    const x = re(a);
    const y = im(a);
    const negIm = sign(y) < 0;
    const reS = x.rational && x.rational[0] === 0n ? "" : format(x);
    const imS = format(negIm ? neg(y) : y);
    if (!reS.endsWith("…") && !imS.endsWith("…")) return joinParts(reS, negIm, imS);
  }
  return a.toString(digits);
}

export { isReal, isComplex, signAt, squareFree, sqrtFloor };
