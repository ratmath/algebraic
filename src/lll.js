// Lattice reduction (LLL, exact rational Gram–Schmidt) and its classic use:
// recognising a decimal as an algebraic number. Given x to D digits, the
// lattice spanned by the rows eᵢ ⊕ round(C·xⁱ) (C ≈ 10^D) has a short vector
// exactly when some small integer polynomial nearly vanishes at x. The
// candidate is then checked exactly: one of its roots must lie in the interval
// the digits allow.

import { reduce, addQ, subQ, mulQ, divQ, toPair, abs, bitLength } from "./q.js";
import { factor } from "./factor.js";
import { isolate, degree } from "./poly.js";

const Z = [0n, 1n];
const dotQ = (u, v) => u.reduce((s, x, i) => addQ(s, mulQ(x, v[i])), Z);
const cmpQ = ([a, b], [c, d]) => {
  const x = a * d - c * b;
  return x > 0n ? 1 : x < 0n ? -1 : 0;
};
const roundQ = ([n, d]) => {
  const t = 2n * n + d;
  const q = t / (2n * d);
  return t % (2n * d) !== 0n && t < 0n ? q - 1n : q;
};

/** LLL-reduce a basis of integer row vectors (BigInt arrays), δ = 3/4. Returns new rows. */
export function lll(rows, delta = [3n, 4n]) {
  const b = rows.map((r) => r.map(BigInt));
  const n = b.length;
  let bs;
  let mu;
  let Bn;
  const gso = () => {
    bs = [];
    mu = Array.from({ length: n }, () => Array(n).fill(Z));
    Bn = [];
    for (let i = 0; i < n; i += 1) {
      let v = b[i].map((x) => [x, 1n]);
      for (let j = 0; j < i; j += 1) {
        mu[i][j] = divQ(dotQ(b[i].map((x) => [x, 1n]), bs[j]), Bn[j]);
        v = v.map((x, k) => subQ(x, mulQ(mu[i][j], bs[j][k])));
      }
      bs.push(v);
      Bn.push(dotQ(v, v));
    }
  };
  gso();
  let k = 1;
  while (k < n) {
    for (let j = k - 1; j >= 0; j -= 1) {
      const q = roundQ(mu[k][j]);
      if (q !== 0n) {
        b[k] = b[k].map((x, i) => x - q * b[j][i]);
        gso();
      }
    }
    // Lovász condition
    const lhs = Bn[k];
    const rhs = mulQ(subQ(delta, mulQ(mu[k][k - 1], mu[k][k - 1])), Bn[k - 1]);
    if (cmpQ(lhs, rhs) >= 0) k += 1;
    else {
      [b[k], b[k - 1]] = [b[k - 1], b[k]];
      gso();
      k = Math.max(k - 1, 1);
    }
  }
  return b;
}

/** Read x as { mid, rad } pairs: a decimal string's last digit gives the radius (a trailing "…" widens it to one unit). */
function readApprox(x) {
  if (typeof x === "string") {
    const s = x.trim().replace(/[…]|\.\.\.$/g, "");
    const loose = s !== x.trim();
    const m = /^(-?)(\d*)\.?(\d*)$/.exec(s);
    if (!m) throw new TypeError(`cannot read ${x} as a decimal`);
    const frac = m[3].length;
    const mid = toPair(s);
    return { mid, rad: reduce(1n, 10n ** BigInt(frac)), digits: (m[2] + m[3]).replace(/^0+/, "").length };
  }
  if (Array.isArray(x) && x.length === 2 && !(typeof x[0] === "bigint")) {
    const lo = toPair(x[0]);
    const hi = toPair(x[1]);
    const mid = reduce(lo[0] * hi[1] + hi[0] * lo[1], 2n * lo[1] * hi[1]);
    const rad = subQ(hi, mid);
    const digits = Math.max(1, Math.floor(-Math.log10(Number(rad[0]) / Number(rad[1]) || 1e-300)));
    return { mid, rad, digits };
  }
  if (typeof x === "number") return readApprox(String(x));
  throw new TypeError("recognize needs a decimal string or an interval [lo, hi]");
}

/**
 * Guess an algebraic number from its digits. Returns
 * { value, poly, degree, height, digitsUsed, digitsNeeded } for the simplest
 * polynomial found (degree ≤ maxDegree, an irreducible factor with a root in
 * the interval the digits allow), or null. `digitsNeeded` is about how many
 * digits a polynomial of that size takes to pin down; a fit is only accepted
 * when the digits given exceed it by `margin` (3), so that it is a real
 * recognition and not a coincidence. It is still a guess: more digits can
 * refute it, never prove it.
 */
export function recognize(x, { maxDegree = 4, margin = 3 } = {}) {
  const { mid, rad, digits } = readApprox(x);
  const lo = subQ(mid, rad);
  const hi = addQ(mid, rad);
  const C = 10n ** BigInt(Math.max(4, digits));
  for (let d = 1; d <= maxDegree; d += 1) {
    const rows = [];
    let pow = [1n, 1n];
    for (let i = 0; i <= d; i += 1) {
      const row = Array(d + 1).fill(0n);
      row[i] = 1n;
      rows.push([...row, roundQ(mulQ(pow, [C, 1n]))]);
      pow = mulQ(pow, mid);
    }
    const red = lll(rows);
    for (const v of red.slice(0, 2)) {
      const coeffs = v.slice(0, d + 1);
      if (coeffs.every((c) => c === 0n) || coeffs[d] === 0n && coeffs.slice(1).every((c) => c === 0n)) continue;
      const found = rootIn(coeffs, lo, hi);
      if (!found) continue;
      const height = found.poly.reduce((m, c) => (abs(c) > m ? abs(c) : m), 0n);
      const deg = degree(found.poly);
      const digitsNeeded = Math.ceil((deg + 1) * Math.log10(Number(height) + 1)) + 1;
      // a fit that spends nearly all the digits is a coincidence, not a recognition
      if (digits - digitsNeeded < margin) continue;
      return { value: found.value, poly: found.poly, degree: deg, height, digitsUsed: digits, digitsNeeded };
    }
  }
  return null;
}

function rootIn(coeffs, lo, hi) {
  let p = coeffs;
  while (p.length > 1 && p[p.length - 1] === 0n) p = p.slice(0, -1);
  if (p.length < 2) return null;
  for (const { factor: f } of factor(p).factors) {
    for (const r of isolate(f)) {
      const v = r.rational;
      if (v) {
        if (cmpQ(lo, v) <= 0 && cmpQ(v, hi) <= 0) return { value: r, poly: f };
        continue;
      }
      for (let bits = 8; bits <= 4 * bitLength(hi[1]) + 64; bits *= 2) {
        r.refine(bits);
        if (cmpQ(r.hi, lo) < 0 || cmpQ(hi, r.lo) < 0) break;
        if (cmpQ(lo, r.lo) <= 0 && cmpQ(r.hi, hi) <= 0) return { value: r, poly: f };
      }
    }
  }
  return null;
}

