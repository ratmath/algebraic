// Number theory on BigInts, the arithmetic underneath finite fields:
// extended gcd, inverses, powers, the Chinese remainder theorem, Jacobi
// symbols, square roots mod p, integer factoring, orders, primitive roots and
// discrete logarithms. Inputs may be numbers or bigints; results are bigints.

import { isProbablePrime } from "./finite.js";

const B = (x) => BigInt(x);
const abs = (a) => (a < 0n ? -a : a);
const mod = (a, m) => ((a % m) + m) % m;

/** { g, x, y } with a·x + b·y = g = gcd(a, b) ≥ 0. */
export function egcd(a, b) {
  a = B(a); b = B(b);
  let [r0, r1, x0, x1, y0, y1] = [a, b, 1n, 0n, 0n, 1n];
  while (r1) {
    const q = r0 / r1;
    [r0, r1, x0, x1, y0, y1] = [r1, r0 - q * r1, x1, x0 - q * x1, y1, y0 - q * y1];
  }
  return r0 < 0n ? { g: -r0, x: -x0, y: -y0 } : { g: r0, x: x0, y: y0 };
}

/** a⁻¹ mod m (throws when gcd(a, m) ≠ 1). */
export function modInverse(a, m) {
  m = B(m);
  const { g, x } = egcd(mod(B(a), m), m);
  if (g !== 1n) throw new RangeError(`${a} has no inverse mod ${m}`);
  return mod(x, m);
}

/** aᵉ mod m, e ≥ 0 (negative e uses the inverse). */
export function modPow(a, e, m) {
  a = B(a); e = B(e); m = B(m);
  if (e < 0n) return modPow(modInverse(a, m), -e, m);
  let r = 1n % m;
  a = mod(a, m);
  while (e) {
    if (e & 1n) r = (r * a) % m;
    a = (a * a) % m;
    e >>= 1n;
  }
  return r;
}

/** x with x ≡ rᵢ (mod mᵢ) for all i, as { x, m } (m the lcm), or null when inconsistent. */
export function crt(residues, moduli) {
  let x = 0n;
  let m = 1n;
  for (let i = 0; i < residues.length; i += 1) {
    const r = B(residues[i]);
    const n = B(moduli[i]);
    const { g, x: u } = egcd(m, n);
    if ((r - x) % g !== 0n) return null;
    const l = (m / g) * n;
    x = mod(x + m * mod(((r - x) / g) * u, n / g), l);
    m = l;
  }
  return { x, m };
}

/** The Jacobi symbol (a/n) for odd n > 0; the Legendre symbol when n is prime. */
export function jacobi(a, n) {
  a = B(a); n = B(n);
  if (n <= 0n || n % 2n === 0n) throw new RangeError("jacobi needs an odd positive n");
  a = mod(a, n);
  let t = 1;
  while (a) {
    while (a % 2n === 0n) {
      a /= 2n;
      const r = n % 8n;
      if (r === 3n || r === 5n) t = -t;
    }
    [a, n] = [n, a];
    if (a % 4n === 3n && n % 4n === 3n) t = -t;
    a %= n;
  }
  return n === 1n ? t : 0;
}

/** A square root of a modulo an odd prime p (Tonelli–Shanks), the smaller of the two; null when a is not a square. */
export function sqrtMod(a, p) {
  a = mod(B(a), B(p)); p = B(p);
  if (a === 0n) return 0n;
  if (p === 2n) return a;
  if (jacobi(a, p) !== 1) return null;
  let q = p - 1n;
  let s = 0n;
  while (q % 2n === 0n) { q /= 2n; s += 1n; }
  let z = 2n;
  while (jacobi(z, p) !== -1) z += 1n;
  let [m, c, t, r] = [s, modPow(z, q, p), modPow(a, q, p), modPow(a, (q + 1n) / 2n, p)];
  while (t !== 1n) {
    let i = 0n;
    let t2 = t;
    while (t2 !== 1n) { t2 = (t2 * t2) % p; i += 1n; }
    const b = modPow(c, 1n << (m - i - 1n), p);
    [m, c, t, r] = [i, (b * b) % p, (t * b * b) % p, (r * b) % p];
  }
  return r < p - r ? r : p - r;
}

function rho(n) {
  if (n % 2n === 0n) return 2n;
  for (let c = 1n; ; c += 1n) {
    let [x, y, d] = [2n, 2n, 1n];
    const f = (v) => (v * v + c) % n;
    while (d === 1n) {
      x = f(x);
      y = f(f(y));
      d = egcd(abs(x - y), n).g;
    }
    if (d !== n) return d;
  }
}

/** Prime factorisation: [{ prime, exponent }] ascending (trial division, then Pollard rho). */
export function factorInteger(n) {
  n = abs(B(n));
  if (n < 2n) return [];
  const found = new Map();
  const push = (p) => found.set(p, (found.get(p) ?? 0) + 1);
  for (const p of [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n]) while (n % p === 0n) { push(p); n /= p; }
  for (let p = 41n; p * p <= n && p < 10000n; p += 2n) while (n % p === 0n) { push(p); n /= p; }
  const stack = n > 1n ? [n] : [];
  while (stack.length) {
    const m = stack.pop();
    if (m === 1n) continue;
    if (isProbablePrime(m)) { push(m); continue; }
    const d = rho(m);
    stack.push(d, m / d);
  }
  return [...found].map(([prime, exponent]) => ({ prime, exponent })).sort((a, b) => (a.prime < b.prime ? -1 : 1));
}

/** Euler's φ(n). */
export function totient(n) {
  let r = B(n);
  for (const { prime } of factorInteger(n)) r = (r / prime) * (prime - 1n);
  return r;
}

/** The order of a in (Z/nZ)ˣ: the least k > 0 with aᵏ ≡ 1. */
export function multiplicativeOrder(a, n) {
  n = B(n);
  if (egcd(B(a), n).g !== 1n) throw new RangeError(`${a} is not a unit mod ${n}`);
  return orderIn(totient(n), (e) => modPow(a, e, n) === 1n);
}

/** Least divisor k of `group` with isOne(k), given that isOne(group). */
export function orderIn(group, isOne) {
  let k = B(group);
  for (const { prime } of factorInteger(k)) while (k % prime === 0n && isOne(k / prime)) k /= prime;
  return k;
}

/** The least primitive root modulo a prime p. */
export function primitiveRoot(p) {
  p = B(p);
  if (!isProbablePrime(p)) throw new RangeError("primitiveRoot needs a prime");
  if (p === 2n) return 1n;
  const qs = factorInteger(p - 1n).map((f) => f.prime);
  for (let g = 2n; ; g += 1n) if (qs.every((q) => modPow(g, (p - 1n) / q, p) !== 1n)) return g;
}

/**
 * Discrete logarithm: the least x ≥ 0 with gˣ = h in a cyclic group of the
 * given order (baby-step giant-step). `ops` = { mul, eq, one, key, inv }; for
 * integers mod p use discreteLogMod. Returns null when there is none.
 */
export function discreteLog(g, h, order, ops) {
  const n = B(order);
  const m = BigInt(Math.ceil(Math.sqrt(Number(n)))) + 1n;
  if (m > 2n ** 22n) throw new RangeError("the group is too large for baby-step giant-step here");
  const table = new Map();
  let e = ops.one;
  for (let j = 0n; j < m; j += 1n) {
    if (!table.has(ops.key(e))) table.set(ops.key(e), j);
    e = ops.mul(e, g);
  }
  const step = ops.inv(e === undefined ? ops.one : power(g, m, ops)); // g^(−m)
  let y = h;
  for (let i = 0n; i < m; i += 1n) {
    const j = table.get(ops.key(y));
    if (j !== undefined) return mod(i * m + j, n);
    y = ops.mul(y, step);
  }
  return null;
}

function power(g, e, ops) {
  let r = ops.one;
  let b = g;
  while (e) {
    if (e & 1n) r = ops.mul(r, b);
    b = ops.mul(b, b);
    e >>= 1n;
  }
  return r;
}

/** The least x ≥ 0 with gˣ ≡ h (mod p), p prime, or null. */
export function discreteLogMod(g, h, p) {
  p = B(p);
  const ops = { one: 1n, mul: (a, b) => (a * b) % p, key: String, inv: (a) => modInverse(a, p) };
  return discreteLog(mod(B(g), p), mod(B(h), p), multiplicativeOrder(g, p), ops);
}

export { isProbablePrime };
