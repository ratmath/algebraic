// Factoring integer polynomials into irreducibles (Zassenhaus): factor the
// square-free part modulo a small prime (Cantor–Zassenhaus), lift the factors
// to a large prime power (Hensel), then find which products of lifted factors
// divide over the integers.
//
// Polynomials are BigInt arrays, constant term first.

import { gcd as gcdInt, abs, sqrtFloor } from "./q.js";
import { squareFree, primitive, degree, derivative } from "./poly.js";

// ---------------------------------------------------------------------------
// Polynomials modulo m (m may be a prime or a prime power)

const md = (x, m) => ((x % m) + m) % m;
const trimP = (a) => {
  let n = a.length;
  while (n > 0 && a[n - 1] === 0n) n -= 1;
  return n === a.length ? a : a.slice(0, n);
};
const modP = (a, m) => trimP(a.map((c) => md(c, m)));

function invMod(a, m) {
  let [r0, r1, s0, s1] = [m, md(a, m), 0n, 1n];
  while (r1) {
    const q = r0 / r1;
    [r0, r1, s0, s1] = [r1, r0 - q * r1, s1, s0 - q * s1];
  }
  if (r0 !== 1n) throw new Error("not invertible");
  return md(s0, m);
}

function addM(a, b, m) {
  const out = Array(Math.max(a.length, b.length)).fill(0n);
  a.forEach((c, i) => { out[i] += c; });
  b.forEach((c, i) => { out[i] += c; });
  return modP(out, m);
}
function subM(a, b, m) {
  const out = Array(Math.max(a.length, b.length)).fill(0n);
  a.forEach((c, i) => { out[i] += c; });
  b.forEach((c, i) => { out[i] -= c; });
  return modP(out, m);
}
function mulM(a, b, m) {
  if (!a.length || !b.length) return [];
  const out = Array(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i += 1) {
    if (a[i] === 0n) continue;
    for (let j = 0; j < b.length; j += 1) out[i + j] += a[i] * b[j];
  }
  return modP(out, m);
}
const scaleM = (a, k, m) => modP(a.map((c) => c * k), m);

/** Quotient and remainder by b, whose leading coefficient is invertible mod m. */
function divmodM(a, b, m) {
  a = modP(a, m);
  b = modP(b, m);
  const db = b.length - 1;
  if (db < 0) throw new Error("division by zero polynomial");
  const il = invMod(b[db], m);
  const r = [...a];
  const q = Array(Math.max(0, a.length - db)).fill(0n);
  for (let k = r.length - 1 - db; k >= 0; k -= 1) {
    const c = md(r[k + db] * il, m);
    q[k] = c;
    if (c === 0n) continue;
    for (let i = 0; i <= db; i += 1) r[k + i] = md(r[k + i] - c * b[i], m);
  }
  return [trimP(q), trimP(r.slice(0, db))];
}
const remM = (a, b, m) => divmodM(a, b, m)[1];
const monicM = (a, m) => (a.length ? scaleM(a, invMod(a[a.length - 1], m), m) : a);

function gcdM(a, b, p) {
  a = modP(a, p);
  b = modP(b, p);
  while (b.length) [a, b] = [b, remM(a, b, p)];
  return monicM(a, p);
}

/** s, t with s·a + t·b = 1 mod p (a, b coprime). */
function bezoutM(a, b, p) {
  let [r0, r1] = [modP(a, p), modP(b, p)];
  let [s0, s1] = [[1n], []];
  let [t0, t1] = [[], [1n]];
  while (r1.length) {
    const [q, r] = divmodM(r0, r1, p);
    [r0, r1] = [r1, r];
    [s0, s1] = [s1, subM(s0, mulM(q, s1, p), p)];
    [t0, t1] = [t1, subM(t0, mulM(q, t1, p), p)];
  }
  const k = invMod(r0[0], p); // r0 is a nonzero constant
  return [scaleM(s0, k, p), scaleM(t0, k, p)];
}

function powModPoly(base, e, f, p) {
  let out = [1n];
  let b = remM(base, f, p);
  while (e > 0n) {
    if (e & 1n) out = remM(mulM(out, b, p), f, p);
    e >>= 1n;
    if (e) b = remM(mulM(b, b, p), f, p);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Cantor–Zassenhaus over F_p (p odd), f monic and square-free

let seed = 12345n;
const rnd = (p) => {
  seed = (seed * 6364136223846793005n + 1442695040888963407n) % 18446744073709551616n;
  return (seed >> 32n) % p; // the high bits: the low bits of a power-of-two LCG cycle quickly
};

function distinctDegree(f, p) {
  const out = [];
  let h = [0n, 1n];
  let i = 0;
  while (f.length - 1 >= 2 * (i + 1)) {
    i += 1;
    h = powModPoly(h, p, f, p);
    const g = gcdM(subM(h, [0n, 1n], p), f, p);
    if (g.length > 1) {
      out.push([g, i]);
      f = divmodM(f, g, p)[0];
      h = remM(h, f, p);
    }
  }
  if (f.length > 1) out.push([monicM(f, p), f.length - 1]);
  return out;
}

function equalDegree(f, d, p) {
  const n = f.length - 1;
  if (n === d) return [f];
  const e = (p ** BigInt(d) - 1n) / 2n;
  for (;;) {
    const a = trimP(Array.from({ length: n }, () => rnd(p)));
    if (a.length < 2) continue;
    let b;
    if (p === 2n) {
      // characteristic 2: the trace a + a² + … + a^(2^(d−1)) splits instead of a^((pᵈ−1)/2)
      let t = remM(a, f, p);
      b = t;
      for (let j = 1; j < d; j += 1) {
        t = remM(mulM(t, t, p), f, p);
        b = addM(b, t, p);
      }
    } else b = subM(powModPoly(a, e, f, p), [1n], p);
    const g = gcdM(b, f, p);
    if (g.length > 1 && g.length < f.length) {
      return [...equalDegree(g, d, p), ...equalDegree(divmodM(f, g, p)[0], d, p)];
    }
  }
}

/** Monic irreducible factors of f mod p; f must be square-free mod p (a RangeError says when it is not). */
export function factorModP(f, p) {
  p = BigInt(p);
  f = monicM(modP(f.map(BigInt), p), p);
  if (f.length < 2) return [];
  const df = modP(derivative(f), p);
  if (!df.length || gcdM(f, df, p).length > 1) throw new RangeError(`not square-free mod ${p}`);
  return distinctDegree(f, p).flatMap(([g, d]) => equalDegree(g, d, p));
}

// ---------------------------------------------------------------------------
// Hensel lifting

/** Lift F ≡ g·h (mod p), g, h monic and coprime, to modulus p^k. F monic mod p^k. */
function henselPair(F, g, h, p, k) {
  const [s, t] = bezoutM(g, h, p);
  let m = p;
  for (let j = 1; j < k; j += 1) {
    const M = m * p;
    const e = subM(F, mulM(g, h, M), M).map((c) => md(c / m, p));
    const sigma = remM(mulM(e, s, p), h, p);
    const tau = divmodM(subM(e, mulM(sigma, g, p), p), h, p)[0];
    g = addM(g, scaleM(tau, m, M), M);
    h = addM(h, scaleM(sigma, m, M), M);
    m = M;
  }
  return [g, h];
}

/** Lift monic factors of F mod p to mod p^k (F monic mod p^k). */
function henselLift(F, factors, p, k) {
  if (factors.length === 1) return [modP(F, p ** BigInt(k))];
  const mid = factors.length >> 1;
  const left = factors.slice(0, mid);
  const right = factors.slice(mid);
  const g0 = left.reduce((a, b) => mulM(a, b, p), [1n]);
  const h0 = right.reduce((a, b) => mulM(a, b, p), [1n]);
  const [g, h] = henselPair(F, g0, h0, p, k);
  return [...henselLift(g, left, p, k), ...henselLift(h, right, p, k)];
}

// ---------------------------------------------------------------------------
// Over the integers

/** Exact quotient a / b over Z, or null when b does not divide a. */
export function divideExact(a, b) {
  const db = b.length - 1;
  const r = [...a];
  const q = Array(Math.max(0, a.length - db)).fill(0n);
  for (let k = r.length - 1 - db; k >= 0; k -= 1) {
    if (r[k + db] % b[db] !== 0n) return null;
    const c = r[k + db] / b[db];
    q[k] = c;
    for (let i = 0; i <= db; i += 1) r[k + i] -= c * b[i];
  }
  return r.slice(0, db).every((c) => c === 0n) ? trimP(q) : null;
}

const symmetric = (a, m) => a.map((c) => (c > m / 2n ? c - m : c));

function* subsets(n, k, start = 0) {
  if (k === 0) {
    yield [];
    return;
  }
  for (let i = start; i <= n - k; i += 1) for (const rest of subsets(n, k - 1, i + 1)) yield [i, ...rest];
}

const SMALL_PRIMES = [3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n, 41n, 43n, 47n, 53n, 59n, 61n, 67n, 71n, 73n, 79n, 83n, 89n, 97n, 101n, 103n, 107n, 109n, 113n];

/** Irreducible factors of a primitive square-free f with positive leading coefficient. */
function factorSquareFree(f) {
  const n = f.length - 1;
  if (n <= 1) return [f];
  if (f[0] === 0n) return [[0n, 1n], ...factorSquareFree(f.slice(1))];
  // a prime keeping f square-free; among a few, the one with fewest factors
  let best = null;
  let tried = 0;
  for (const p of SMALL_PRIMES) {
    if (f[n] % p === 0n) continue;
    const fp = modP(f, p);
    if (gcdM(fp, modP(derivative(f), p), p).length > 1) continue;
    const fac = factorModP(fp, p);
    if (!best || fac.length < best.fac.length) best = { p, fac };
    if (fac.length === 1 || ++tried >= 5) break;
  }
  if (!best) throw new Error("factor: no suitable prime");
  if (best.fac.length === 1) return [f];
  const { p } = best;
  // Mignotte-style bound on the coefficients of any factor, times the leading coefficient
  const norm2 = sqrtFloor(f.reduce((s, c) => s + c * c, 0n)) + 1n;
  const B = 2n ** BigInt(n) * norm2 * abs(f[n]);
  let k = 1;
  let M = p;
  while (M <= 2n * B) {
    M *= p;
    k += 1;
  }
  const lc = f[n];
  const F = scaleM(f, invMod(lc, M), M);
  let lifted = henselLift(F, best.fac, p, k);
  // recombine: products of lifted factors that divide over Z
  const out = [];
  let g = f;
  let size = 1;
  while (2 * size <= lifted.length) {
    let found = false;
    for (const S of subsets(lifted.length, size)) {
      const lcg = g[g.length - 1];
      const prod = S.reduce((a, i) => mulM(a, lifted[i], M), [md(lcg, M)]);
      const cand = primitive(symmetric(prod, M));
      const q = divideExact(g, cand);
      if (!q) continue;
      out.push(cand);
      g = q;
      lifted = lifted.filter((_, i) => !S.includes(i));
      found = true;
      break;
    }
    if (!found) size += 1;
  }
  if (g.length > 1) out.push(primitive(g));
  return out;
}

/**
 * Factor an integer polynomial: { content, factors: [{ factor, multiplicity }] }
 * with each factor primitive, irreducible over Q, positive leading coefficient.
 */
export function factor(p) {
  p = trimP(p.map(BigInt));
  if (!p.length) return { content: 0n, factors: [] };
  let content = p.reduce((g, c) => gcdInt(g, c), 0n);
  if (p[p.length - 1] < 0n) content = -content;
  if (degree(p) < 1) return { content: p[0], factors: [] };
  const factors = [];
  for (const { factor: s, multiplicity } of squareFree(p)) {
    for (const g of factorSquareFree(primitive(s))) factors.push({ factor: g, multiplicity });
  }
  factors.sort((a, b) => a.factor.length - b.factor.length || (a.factor.join() < b.factor.join() ? -1 : 1));
  return { content, factors };
}

/** True when p (degree ≥ 1) is irreducible over Q. */
export const isIrreducible = (p) => {
  const { factors } = factor(p);
  return factors.length === 1 && factors[0].multiplicity === 1;
};

export { modP, mulM, gcdM, invMod };
