// Exact values of cos and sin at rational multiples of π, roots of unity and
// cyclotomic polynomials.
//
// cos(2πj/n) for j = 0 … ⌊n/2⌋ are exactly the distinct roots of Tₙ(x) − 1
// (Tₙ the Chebyshev polynomial), and they decrease with j. So the j-th largest
// real root of Tₙ − 1 is cos(2πj/n): no floating point and no π is needed to
// say which root is which.

import { isolate, primitive, degree } from "./poly.js";
import { reduce, gcd } from "./q.js";
import { fromRational, simplify as simplifyReal } from "./arith.js";
import { certifyComplexRoots } from "./complex.js";
import { divideExact } from "./factor.js";

const MAX_N = 2000;

/** Tₙ(x), the Chebyshev polynomial of the first kind (BigInt, constant first). */
export function chebyshevT(n) {
  let a = [1n];
  let b = [0n, 1n];
  if (n === 0) return a;
  for (let k = 1; k < n; k += 1) {
    const c = Array(b.length + 1).fill(0n);
    b.forEach((x, i) => { c[i + 1] += 2n * x; });
    a.forEach((x, i) => { c[i] -= x; });
    [a, b] = [b, c];
  }
  return b;
}

/** The n-th cyclotomic polynomial Φₙ (BigInt, constant first). */
export function cyclotomic(n) {
  if (!Number.isInteger(n) || n < 1) throw new RangeError("cyclotomic needs n ≥ 1");
  let p = [-1n, ...Array(n - 1).fill(0n), 1n]; // xⁿ − 1
  for (let d = 1; d < n; d += 1) if (n % d === 0) p = divideExact(p, cyclotomic(d));
  return p;
}

/** cos(2πk/n) exactly, as an AlgebraicReal on its minimal polynomial. */
function cos2pi(k, n) {
  if (n > MAX_N) throw new RangeError(`denominators up to ${MAX_N / 2} are supported`);
  k = ((k % n) + n) % n;
  const g = gcd(BigInt(k), BigInt(n));
  k = Number(BigInt(k) / (g || 1n));
  n = Number(BigInt(n) / (g || 1n));
  const j = Math.min(k, n - k);
  // the rational cases: cos 0, cos π, cos ±π/3 (= 1/2), cos ±2π/3, cos π/2
  const known = { "1": [1n, 1n], "2": [-1n, 1n], "3": [-1n, 2n], "4": [0n, 1n], "6": [1n, 2n] };
  if (j === 0) return fromRational([1n, 1n]);
  if (known[n]) return fromRational(known[n]);
  // cos(2πj/n) for j coprime to n, 0 < j < n/2, are the roots of cosMinimal(n), descending in j
  const all = isolate(cosMinimal(n)).sort((x, y) => cmpReal(y, x));
  let idx = 0;
  for (let t = 1; t < j; t += 1) if (gcd(BigInt(t), BigInt(n)) === 1n) idx += 1;
  return all[idx];
}

/**
 * The minimal polynomial of cos(2π/n), n ≥ 3: Φₙ is palindromic, so
 * x^(−h)Φₙ(x) is a polynomial in y = x + 1/x = 2cos θ, built from
 * x^k + x^(−k) = Dₖ(y) with D₀ = 2, D₁ = y, Dₖ₊₁ = y·Dₖ − Dₖ₋₁; then y = 2c.
 */
export function cosMinimal(n) {
  const phi = cyclotomic(n);
  const h = (phi.length - 1) / 2;
  const addTo = (acc, p, c) => p.forEach((x, i) => { acc[i] = (acc[i] ?? 0n) + c * x; });
  let D0 = [2n];
  let D1 = [0n, 1n];
  const psi = [phi[h]];
  for (let k = 1; k <= h; k += 1) {
    addTo(psi, D1, phi[h + k]);
    const next = [0n, ...D1];
    D0.forEach((x, i) => { next[i] -= x; });
    [D0, D1] = [D1, next];
  }
  // y = 2c
  return primitive(psi.map((x, i) => x * 2n ** BigInt(i)));
}

function cmpReal(a, b) {
  // a and b are distinct roots, so their isolating intervals separate
  for (let bits = 8; ; bits *= 2) {
    const [al, ah] = enclose(a, bits);
    const [bl, bh] = enclose(b, bits);
    if (lt(ah, bl)) return -1;
    if (lt(bh, al)) return 1;
    if (eq(al, bl) && eq(ah, bh) && a.rational && b.rational) return 0;
  }
}
const enclose = (a, bits) => (a.rational ? [a.rational, a.rational] : (a.refine(bits), [a.lo, a.hi]));
const lt = ([a, b], [c, d]) => a * d < c * b;
const eq = ([a, b], [c, d]) => a * d === c * b;

/** A rational r (pair, Rational, number or string) as [n, d]. */
function ratio(r) {
  const p = fromRational(r).rational;
  return p;
}

/** cos(πr) for rational r, exactly. */
export function cosPi(r) {
  const [p, q] = ratio(r);
  // πp/q = 2π·p/(2q)
  return cos2pi(Number(((p % (2n * q)) + 2n * q) % (2n * q)), Number(2n * q));
}

/** sin(πr) for rational r, exactly: cos(π(1/2 − r)). */
export function sinPi(r) {
  const [p, q] = ratio(r);
  const [a, b] = reduce(q - 2n * p, 2n * q);
  return cosPi([a, b]);
}

/** tan(πr) for rational r, exactly (throws at the poles). */
export function tanPi(r) {
  let [p, q] = ratio(r);
  if (q === 2n) throw new RangeError("tan has a pole there");
  if (q > BigInt(MAX_N / 2)) throw new RangeError(`denominators up to ${MAX_N / 2} are supported`);
  // (1 + i·t)^q is real when t = tan(πj/q), so the q (or q − 1) values
  // tan(πj/q), −q/2 < j < q/2, are exactly the roots of
  //   Im (1 + i·t)^q = Σ_{k odd} (−1)^((k−1)/2) C(q, k) t^k,
  // and they increase with j: the root is picked by counting, as for cos.
  const n = Number(q);
  let j = Number(((p % q) + q) % q);
  if (2 * j > n) j -= n;
  if (j === 0) return fromRational([0n, 1n]);
  const P = [];
  let c = 1n; // C(n, k)
  for (let k = 0; k <= n; k += 1) {
    P.push(k % 2 === 0 ? 0n : (k % 4 === 1 ? c : -c));
    c = (c * BigInt(n - k)) / BigInt(k + 1);
  }
  while (P[P.length - 1] === 0n) P.pop();
  P.shift(); // drop the root t = 0 (j = 0); the rest is a polynomial in t², handled as is
  const all = isolate(P).sort(cmpReal);
  const half = Math.floor((n - 1) / 2);
  const idx = j < 0 ? j + half : j + half - 1;
  return simplifyReal(all[idx]);
}

/**
 * The root of unity e^(2πik/n), exactly: an AlgebraicReal (±1) or a
 * ComplexAlgebraic on the cyclotomic polynomial Φ_(n/gcd(k,n)).
 */
export function rootOfUnity(n, k = 1) {
  if (!Number.isInteger(n) || n < 1) throw new RangeError("rootOfUnity needs n ≥ 1");
  const c = cos2pi(k, n);
  const s = cos2pi(n - 4 * k, 4 * n); // sin(2πk/n) = cos(2π(n − 4k)/(4n))
  if (s.rational && s.rational[0] === 0n) return c;
  // the primitive m-th roots of unity are the roots of Φₘ; take the one in the box cos × sin
  const g = Number(gcd(BigInt(((k % n) + n) % n), BigInt(n)));
  const phi = cyclotomic(n / g);
  const cands = certifyComplexRoots(phi, degree(phi));
  if (!cands) throw new Error(`could not certify the roots of Φ${n / g}`);
  for (let bits = 8; bits <= 8192; bits *= 2) {
    const C = enclose(c, bits);
    const S = enclose(s, bits);
    const hits = cands.filter((z) => {
      z.refine(bits);
      return meets(z.re, C) && meets(z.im, S);
    });
    if (hits.length === 1) return hits[0];
  }
  throw new Error("could not separate the roots of unity");
}
const meets = ([a, b], [c, d]) => !lt(b, c) && !lt(d, a);

/** All n-th roots of unity, e^(2πik/n) for k = 0 … n − 1. */
export function rootsOfUnity(n) {
  return Array.from({ length: n }, (_, k) => rootOfUnity(n, k));
}

export { MAX_N };
