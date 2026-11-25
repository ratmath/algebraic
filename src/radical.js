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
  const m = /^\d*√/.exec(c);
  if (m) return `${c.slice(0, m[0].length - 1)}i${c.slice(m[0].length - 1)}`;
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

/** √((p + q√d)/r) written as k√(A + B√d)/C with integers, as small as square-free splitting allows. */
export function nestedRadical({ p, q, d, r }, negative = false) {
  // smallest C with C²·p/r and C²·q/r integers
  const rr = r / (gcd(gcd(p, q), r) || 1n);
  const [s0, t0] = squarefreeSplit(rr);
  const C = s0 * t0;
  let A = (C * C * p) / r;
  let B = (C * C * q) / r;
  // pull a square out of the radicand
  const [k] = squarefreeSplit(gcd(A, B) || 1n);
  A /= k * k;
  B /= k * k;
  // denest: √(A + B√d) = √((A + c)/2) ± √((A − c)/2) when A² − B²d = c²
  const N = A * A - B * B * d;
  if (d > 0n && N >= 0n) {
    const c = sqrtFloor(N);
    if (c * c === N) {
      const [s1, t1] = squarefreeSplit(2n * (A + c));
      const [s2, t2] = squarefreeSplit(2n * (A - c));
      let terms = [[k * s1, t1], [(B < 0n ? -1n : 1n) * k * s2, t2]].filter(([a, t]) => a !== 0n && t !== 0n);
      if (terms.length === 2 && terms[0][1] === terms[1][1]) terms = [[terms[0][0] + terms[1][0], terms[0][1]]];
      if (negative) terms = terms.map(([a, t]) => [-a, t]);
      let den = 2n * C;
      const g = terms.reduce((acc, [a]) => gcd(acc, a), den) || 1n;
      terms = terms.map(([a, t]) => [a / g, t]);
      den /= g;
      const rad = ([a, t], first) => {
        const mag = a < 0n ? -a : a;
        const body = t === 1n ? `${mag}` : `${mag === 1n ? "" : mag}√${t}`;
        return first ? `${a < 0n ? "−" : ""}${body}` : ` ${a < 0n ? "−" : "+"} ${body}`;
      };
      const num = terms.map((tm, i) => rad(tm, i === 0)).join("");
      if (den === 1n) return num;
      return terms.length > 1 ? `(${num})/${den}` : `${num}/${den}`;
    }
  }
  const [kk, CC] = (() => {
    const g = gcd(k, C) || 1n;
    return [k / g, C / g];
  })();
  const inner = formatReal({ p: A, q: B, d, r: 1n });
  return `${negative ? "−" : ""}${kk === 1n ? "" : kk}√(${inner})${CC === 1n ? "" : `/${CC}`}`;
}

const ROOT = { 3: "∛", 4: "∜" };
const SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹";
const rootSign = (n) => ROOT[n] ?? `${[...String(n)].map((c) => SUP[c]).join("")}√`;

/**
 * A root of an irreducible polynomial of degree n ≥ 3 that is a shifted pure
 * power, c·(x − p)ⁿ = e: the parts { p, k, den, m, n } of p + (k/den)·ⁿ√m with
 * m free of n-th powers of small primes, or null when the polynomial is not
 * of that shape. `side` is +1 when the root lies above p, −1 below.
 */
export function pureRootParts(poly) {
  const n = poly.length - 1;
  if (n < 3) return null;
  const add = ([a, b], [c, d]) => reduce(a * d + c * b, b * d);
  const mul = ([a, b], [c, d]) => reduce(a * c, b * d);
  const p = reduce(-poly[n - 1], BigInt(n) * poly[n]);
  // g(y) = f(y + p) by repeated synthetic division
  const g = poly.map((c) => [c, 1n]);
  for (let i = 0; i < n; i += 1) for (let j = n - 1; j >= i; j -= 1) g[j] = add(g[j], mul(g[j + 1], p));
  for (let i = 1; i < n; i += 1) if (g[i][0] !== 0n) return null;
  let [N, D] = reduce(-g[0][0] * g[n][1], g[0][1] * g[n][0]); // yⁿ = N/D
  if (D < 0n) [N, D] = [-N, -D];
  if (N === 0n || (n % 2 === 0 && N < 0n)) return null;
  const negative = N < 0n;
  let m = abs(N) * D ** BigInt(n - 1); // ⁿ√(N/D) = ⁿ√(N·Dⁿ⁻¹)/D
  let k = 1n;
  for (let q = 2n; q <= 1000n && q ** BigInt(n) <= m; q += 1n) {
    const qn = q ** BigInt(n);
    while (m % qn === 0n) {
      m /= qn;
      k *= q;
    }
  }
  const [num, den] = reduce(k, D);
  return { p, k: num, den, m, n, negative };
}

/** "∛2", "−∛4/2", "1 + 2⁵√3": the text for parts from pureRootParts; `side` picks the sign for even n. */
export function formatPureRoot({ p, k, den, m, n, negative }, side) {
  const minus = n % 2 ? negative : side < 0;
  const body = `${k === 1n ? "" : k}${rootSign(n)}${m}${den === 1n ? "" : `/${den}`}`;
  if (p[0] === 0n) return `${minus ? "−" : ""}${body}`;
  return `${fmtQ(p)} ${minus ? "−" : "+"} ${body}`;
}
