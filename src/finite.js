// Finite fields: integers mod a prime p, and GF(p^k) as polynomials modulo an
// irreducible polynomial over F_p. Each field is an object with
// { zero, one, add, sub, mul, div, neg, inv, sign, from } so it plugs into
// field-generic algorithms (for example @ratmath/linalg's `field` module).

import { factorModP } from "./factor.js";
import { multiplicativeOrder, primitiveRoot, discreteLog, discreteLogMod, orderIn, sqrtMod } from "./nt.js";

export function primeField(p) {
  p = BigInt(p);
  if (p < 2n || !isProbablePrime(p)) throw new RangeError(`${p} is not prime`);
  const mod = (x) => ((x % p) + p) % p;
  const inv = (a) => {
    if (a === 0n) throw new RangeError(`division by zero mod ${p}`);
    let [r0, r1, s0, s1] = [p, a, 0n, 1n];
    while (r1) {
      const q = r0 / r1;
      [r0, r1, s0, s1] = [r1, r0 - q * r1, s1, s0 - q * s1];
    }
    return mod(s0);
  };
  return {
    p,
    zero: 0n,
    one: 1n,
    add: (a, b) => mod(a + b),
    sub: (a, b) => mod(a - b),
    mul: (a, b) => mod(a * b),
    div: (a, b) => mod(a * inv(b)),
    neg: (a) => mod(-a),
    inv,
    sign: (a) => (a === 0n ? 0 : 1),
    size: p,
    pow: (a, e) => { let r = 1n; let b = mod(a); e = BigInt(e); if (e < 0n) { b = inv(b); e = -e; } while (e) { if (e & 1n) r = (r * b) % p; b = (b * b) % p; e >>= 1n; } return r; },
    equals: (a, b) => mod(a) === mod(b),
    /** Multiplicative order of a nonzero a. */
    order: (a) => multiplicativeOrder(mod(a), p),
    /** A generator of the multiplicative group (the least primitive root). */
    primitive: () => primitiveRoot(p),
    /** The least x with gˣ = a, g a primitive root by default; null if none. */
    log: (a, g = primitiveRoot(p)) => discreteLogMod(g, mod(a), p),
    /** A square root of a, or null. */
    sqrt: (a) => (p === 2n ? mod(a) : sqrtMod(mod(a), p)),
    elements: () => Array.from({ length: Number(p) }, (_, i) => BigInt(i)),
    format: (a) => `${mod(a)}`,
    from(x) {
      if (typeof x === "bigint") return mod(x);
      if (typeof x === "number") return mod(BigInt(x));
      const d = mod(x.denominator);
      if (d === 0n) throw new RangeError(`${x} has a denominator divisible by ${p}, so it has no value mod ${p}`);
      return mod(mod(x.numerator) * inv(d));
    },
  };
}

export function isProbablePrime(n) {
  if (n < 4n) return n >= 2n;
  if (n % 2n === 0n) return false;
  const powmod = (b, e, m) => {
    let r = 1n;
    b %= m;
    while (e) {
      if (e & 1n) r = (r * b) % m;
      b = (b * b) % m;
      e >>= 1n;
    }
    return r;
  };
  let d = n - 1n;
  let s = 0;
  while (d % 2n === 0n) {
    d /= 2n;
    s += 1;
  }
  for (const a of [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n]) {
    if (a % n === 0n) continue;
    let x = powmod(a, d, n);
    if (x === 1n || x === n - 1n) continue;
    let composite = true;
    for (let i = 1; i < s; i += 1) {
      x = (x * x) % n;
      if (x === n - 1n) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}


/**
 * GF(p^k): elements are arrays of k BigInt coefficients (constant first) of a
 * polynomial modulo `modulus`, a monic irreducible of degree k over F_p
 * (found automatically when not given).
 */
export function galoisField(p, k, modulus = null) {
  const F = primeField(p);
  p = F.p;
  if (!Number.isInteger(k) || k < 1) throw new RangeError("galoisField needs a degree k ≥ 1");
  if (k === 1) return { ...F, k: 1, size: p, from: (x) => F.from(x) };
  const md = (x) => ((x % p) + p) % p;
  if (!modulus) {
    // the first monic irreducible of degree k in a fixed order
    for (let n = 0n; ; n += 1n) {
      const c = [];
      let t = n;
      for (let i = 0; i < k; i += 1) {
        c.push(t % p);
        t /= p;
      }
      if (c[0] === 0n) continue;
      const cand = [...c, 1n];
      if (p === 2n ? irreducibleMod2(cand) : factorModP(cand, p).length === 1) {
        modulus = cand;
        break;
      }
    }
  }
  modulus = modulus.map((c) => md(BigInt(c)));
  const reduceP = (a) => {
    const r = a.map(md);
    for (let i = r.length - 1; i >= k; i -= 1) {
      const c = r[i];
      if (c === 0n) continue;
      for (let j = 0; j <= k; j += 1) r[i - k + j] = md(r[i - k + j] - c * modulus[j]);
    }
    const out = r.slice(0, k);
    while (out.length < k) out.push(0n);
    return out;
  };
  const mul = (a, b) => {
    const out = Array(2 * k - 1).fill(0n);
    for (let i = 0; i < k; i += 1) for (let j = 0; j < k; j += 1) out[i + j] += a[i] * b[j];
    return reduceP(out);
  };
  const isZ = (a) => a.every((c) => c === 0n);
  const pow = (a, e) => {
    let out = [1n, ...Array(k - 1).fill(0n)];
    let b = a;
    while (e > 0n) {
      if (e & 1n) out = mul(out, b);
      e >>= 1n;
      if (e) b = mul(b, b);
    }
    return out;
  };
  const size = p ** BigInt(k);
  const inv = (a) => {
    if (isZ(a)) throw new RangeError("division by zero in GF(p^k)");
    return pow(a, size - 2n);
  };
  const zero = Array(k).fill(0n);
  const one = [1n, ...Array(k - 1).fill(0n)];
  const eq = (a, b) => a.every((c, i) => c === b[i]);
  const fromIndex = (n) => {
    const c = [];
    for (let i = 0; i < k; i += 1) {
      c.push(n % p);
      n /= p;
    }
    return c;
  };
  const order = (a) => {
    if (isZ(a)) throw new RangeError("0 has no multiplicative order");
    return orderIn(size - 1n, (e) => eq(pow(a, e), one));
  };
  let prim = null;
  const primitive = () => {
    if (!prim) for (let n = 1n; n < size; n += 1n) if (order(fromIndex(n)) === size - 1n) { prim = fromIndex(n); break; }
    return prim;
  };
  return {
    p,
    k,
    size,
    modulus,
    zero,
    one,
    gen: [0n, 1n, ...Array(k - 2).fill(0n)],
    add: (a, b) => a.map((c, i) => md(c + b[i])),
    sub: (a, b) => a.map((c, i) => md(c - b[i])),
    mul,
    div: (a, b) => mul(a, inv(b)),
    neg: (a) => a.map((c) => md(-c)),
    inv,
    pow: (a, e) => (e < 0 ? pow(inv(a), BigInt(-e)) : pow(a, BigInt(e))),
    sign: (a) => (isZ(a) ? 0 : 1),
    equals: eq,
    order,
    primitive,
    log: (a, g = primitive()) => discreteLog(g, a, order(g), { one, mul, key: (x) => x.join(","), inv }),
    elements: () => {
      if (size > 100000n) throw new RangeError("too many elements to list");
      return Array.from({ length: Number(size) }, (_, n) => fromIndex(BigInt(n)));
    },
    from: (x) => [F.from(x), ...Array(k - 1).fill(0n)],
    format: (a, v = "a") => {
      const parts = [];
      for (let i = k - 1; i >= 0; i -= 1) {
        if (a[i] === 0n) continue;
        const sup = String(i).split("").map((d) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[d]).join("");
        parts.push(i === 0 ? `${a[i]}` : `${a[i] === 1n ? "" : a[i]}${v}${i > 1 ? sup : ""}`);
      }
      return parts.join(" + ") || "0";
    },
  };
}

function irreducibleMod2(f) {
  // brute force for p = 2: no factor of degree ≤ deg/2
  const n = f.length - 1;
  const asInt = (a) => a.reduce((acc, c, i) => acc | (c === 1n ? 1n << BigInt(i) : 0n), 0n);
  const F = asInt(f);
  const deg = (x) => x.toString(2).length - 1;
  const rem = (a, b) => {
    const db = deg(b);
    while (a && deg(a) >= db) a ^= b << BigInt(deg(a) - db);
    return a;
  };
  for (let g = 2n; deg(g) <= n / 2; g += 1n) if (rem(F, g) === 0n) return false;
  return true;
}
