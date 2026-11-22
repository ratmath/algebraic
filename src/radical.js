// Exact radical forms for display: a root of a quadratic as (p + q√d)/r,
// with d square-free, or as a Gaussian-style a + b·i when d < 0.
//
// These helpers only read an integer polynomial and an enclosure, so they sit
// below arith.js and can be used by AlgebraicReal and ComplexAlgebraic alike.

import { gcd, abs, reduce, sqrtFloor, pairToString, pairToNumber } from "./q.js";

/** D = s²·d with d square-free (sign kept on d). Trial division to 10⁵, then a perfect-square test of the rest. */
export function squarefreeSplit(D) {
  if (D === 0n) return [0n, 0n];
  const sg = D < 0n ? -1n : 1n;
  let n = abs(D);
  let s = 1n;
  let d = 1n;
  for (let p = 2n; p <= 100000n && p * p <= n; p += p === 2n ? 1n : 2n) {
    while (n % (p * p) === 0n) {
      n /= p * p;
      s *= p;
    }
    if (n % p === 0n) {
      n /= p;
      d *= p;
    }
  }
  const r = sqrtFloor(n);
  if (r * r === n) s *= r;
  else d *= n;
  return [s, sg * d];
}

const cmpPair = ([a, b], [c, d]) => {
  const x = a * d - c * b;
  return x > 0n ? 1 : x < 0n ? -1 : 0;
};

/**
 * The root of the quadratic poly = [c, b, a] (a ≠ 0) that `side` names
 * (+1 for the "+√" root, −1 for the "−√" one) as { p, q, d, r }: (p + q√d)/r,
 * gcd(p, q, r) = 1, r > 0, d square-free (negative for complex roots).
 */
export function quadraticParts([c, b, a], side) {
  const D = b * b - 4n * a * c;
  const [s, d] = squarefreeSplit(D);
  let p = -b;
  let q = BigInt(side) * s;
  let r = 2n * a;
  if (r < 0n) {
    p = -p;
    q = -q;
    r = -r;
  }
  const g = gcd(gcd(p, q), r) || 1n;
  return { p: p / g, q: q / g, d, r: r / g };
}

/**
 * Which root of a real quadratic the interval (lo, hi) holds: +1 if it lies
 * above the vertex −b/2a (the +√D/2a root when a > 0), else −1. `refine(bits)`
 * shrinks (lo, hi) when needed.
 */
export function realSide([, b, a], get, refine) {
  const v = reduce(-b, 2n * a);
  for (let bits = 8; bits <= 4096; bits *= 2) {
    const [lo, hi] = get();
    const up = cmpPair(lo, v) > 0;
    const down = cmpPair(hi, v) < 0;
    if (up || down) return (up ? 1 : -1) * (a > 0n ? 1 : -1);
    refine(bits);
  }
  throw new Error("could not place the root beside the vertex");
}

const SQRT = "√";
const radical = (d) => `${SQRT}${d}`;

/** "p", "q√d", "(p + q√d)/r", … for a real quadratic irrational (d > 0). */
export function formatReal({ p, q, d, r }) {
  if (q === 0n || d === 0n) return fmtQ(reduce(p, r));
  if (d === 1n) return fmtQ(reduce(p + q, r));
  const qs = abs(q) === 1n ? radical(d) : `${abs(q)}${radical(d)}`;
  if (p === 0n) return `${q < 0n ? "−" : ""}${qs}${r === 1n ? "" : `/${r}`}`;
  const body = `${fmtInt(p)} ${q < 0n ? "−" : "+"} ${qs}`;
  return r === 1n ? body : `(${body})/${r}`;
}

/** Real and imaginary parts of (p + q√d)/r with d < 0, as "a + bi" (i√m, 3i/4, …). */
export function formatComplex({ p, q, d, r }) {
  const reS = p === 0n ? "" : fmtQ(reduce(p, r));
  const imAbs = formatReal({ p: 0n, q: q < 0n ? -q : q, d: -d, r });
  return joinParts(reS, q < 0n, imAbs);
}

/** i times a positive coefficient written by formatReal: "i", "3i", "3i/4", "i√3/2", "i(1 + √5)/4". */
export function imaginary(c) {
  if (c === "1") return "i";
  if (c.startsWith("(")) return `i${c}`;
  if (c.includes(" ")) return `i(${c})`;
  if (c.startsWith("1/")) return `i${c.slice(1)}`;
  const k = c.search(/[√/]/);
  return k < 0 ? `${c}i` : `${c.slice(0, k)}i${c.slice(k)}`;
}

/** "a + bi" from a real-part string ("" for 0), the sign of the imaginary part and its absolute coefficient. */
export function joinParts(reS, negIm, imAbs) {
  const imS = imaginary(imAbs);
  return reS ? `${reS} ${negIm ? "−" : "+"} ${imS}` : `${negIm ? "−" : ""}${imS}`;
}

const fmtInt = (n) => (n < 0n ? `−${-n}` : `${n}`);
export const fmtQ = ([n, d]) => (n < 0n ? "−" : "") + pairToString([abs(n), d]);

/**
 * A decimal for a value known to lie in [lo, hi] (pairs): `digits` significant
 * digits, the last one off by at most one unit. `refine(bits)` must shrink the
 * interval when asked.
 */
export function certifiedDecimal(get, refine, digits = 12) {
  refine(24);
  let [lo, hi] = get();
  const mag = Math.abs((pairToNumber(...lo) + pairToNumber(...hi)) / 2);
  const lead = mag > 0 ? Math.floor(Math.log10(mag)) : 0;
  const frac = Math.min(Math.max(0, digits - 1 - lead), digits + 6);
  const bits = Math.ceil((frac + 2) * 3.33) + 2;
  refine(bits);
  [lo, hi] = get();
  const s = 10n ** BigInt(frac);
  // round the midpoint to `frac` fractional digits
  const num = (lo[0] * hi[1] + hi[0] * lo[1]) * s;
  const den = 2n * lo[1] * hi[1];
  let k = num >= 0n ? (2n * num + den) / (2n * den) : -((-2n * num + den) / (2n * den));
  const neg = k < 0n;
  if (neg) k = -k;
  let str = k.toString().padStart(frac + 1, "0");
  if (frac) str = `${str.slice(0, -frac)}.${str.slice(-frac)}`;
  return `${neg && /[1-9]/.test(str) ? "−" : ""}${str}`;
}
