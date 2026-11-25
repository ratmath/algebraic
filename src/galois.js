// Galois groups of irreducible polynomials of degree ≤ 5 over Q.
//
// Degree 3 and 4 use the discriminant and the resolvent cubic (with the
// Kappe–Warren test for C4 against D4). Degree 5 uses resolvents computed
// from certified root enclosures: a resolvent's coefficients are integers
// (symmetric functions of the roots of a monic integer polynomial), so
// interval arithmetic on the roots, refined until each coefficient's interval
// holds exactly one integer, gives the resolvent exactly.
//
//   F20 resolvent: θ = (Σ_pentagon xᵢxⱼ − Σ_pentagram xᵢxⱼ)², six values.
//   C5 inside F20: φ = Σ xᵥ₍ᵢ₎²·xᵥ₍ᵢ₊₁₎ along the oriented pentagon, four values.
// G ⊆ F20 (solvable) iff the sextic has a rational root; then G ⊆ C5 iff the
// quartic for that pentagon does (both resolvents square-free, else a
// Tschirnhaus transformation y = x² + kx is applied first).

import { factor } from "./factor.js";
import { primitive, degree, squareFree } from "./poly.js";
import { reduce, subQ, mulQ, sqrtFloor, toPair } from "./q.js";
import { ivAdd, ivSub, ivMul } from "./interval.js";
import { roots } from "./arith.js";
import { discriminant } from "./polys.js";

const GROUPS = {
  1: { "1": [1, "1T1"] },
  2: { C2: [2, "2T1"] },
  3: { C3: [3, "3T1"], S3: [6, "3T2"] },
  4: { C4: [4, "4T1"], V4: [4, "4T2"], D4: [8, "4T3"], A4: [12, "4T4"], S4: [24, "4T5"] },
  5: { C5: [5, "5T1"], D5: [10, "5T2"], F20: [20, "5T3"], A5: [60, "5T4"], S5: [120, "5T5"] },
};
const SOLVABLE = new Set(["1", "C2", "C3", "S3", "C4", "V4", "D4", "A4", "S4", "C5", "D5", "F20"]);

const isSquareInt = (n) => n >= 0n && sqrtFloor(n) ** 2n === n;
const isSquareQ = ([n, d]) => isSquareInt(n) && isSquareInt(d);

/** Monic integer polynomial with the same splitting field: a^(n−1)·f(x/a). */
function monic(f) {
  const n = f.length - 1;
  const a = f[n];
  return f.map((c, i) => c * a ** BigInt(Math.max(0, n - 1 - i))).map((c, i) => (i === n ? 1n : c));
}

function result(n, name, extra = {}) {
  const [order, label] = GROUPS[n][name];
  return { degree: n, name, order, label, solvable: SOLVABLE.has(name), ...extra };
}

/**
 * The Galois group of an integer or rational polynomial of degree ≤ 5:
 * { degree, name ("S4", "D5", …), order, label (transitive-group "nTk"),
 * solvable, discriminant (Core Rational), discriminantIsSquare }. A reducible
 * input gives { reducible: true, factors: [group of each irreducible factor] }.
 */
export function galoisGroup(coeffs) {
  const pairs = coeffs.map((c) => (Array.isArray(c) ? reduce(BigInt(c[0]), BigInt(c[1])) : toPair(c)));
  let L = 1n;
  const g = (a, b) => (b ? g(b, a % b) : a);
  for (const [, d] of pairs) L = (L * d) / g(L, d);
  let f = primitive(pairs.map(([n, d]) => n * (L / d)));
  const { factors } = factor(f);
  if (factors.length !== 1 || factors[0].multiplicity !== 1) {
    return { reducible: true, factors: factors.map(({ factor: h, multiplicity }) => ({ factor: h, multiplicity, group: galoisGroup(h) })) };
  }
  const n = degree(f);
  if (n > 7) throw new RangeError("Galois groups are computed up to degree 7");
  if (n > 5) {
    const a = new GaloisAction(f);
    return { degree: n, name: a.name, order: a.order, label: a.label, solvable: a.solvable, discriminant: a.discriminant, discriminantIsSquare: a.discriminantIsSquare };
  }
  if (n === 1) return result(1, "1");
  const disc = discriminant(f.map((c) => [c, 1n]));
  const D = [disc.numerator, disc.denominator];
  const sq = isSquareQ(D);
  const info = { discriminant: disc, discriminantIsSquare: sq };
  if (n === 2) return result(2, "C2", info);
  if (n === 3) return result(3, sq ? "C3" : "S3", info);
  f = monic(f);
  if (n === 4) return quartic(f, D, sq, info);
  return quintic(f, sq, info);
}

function quartic(f, D, sq, info) {
  const [e, d, c, b] = f; // x⁴ + b x³ + c x² + d x + e
  const R = [-(b * b * e - 4n * c * e + d * d), b * d - 4n * e, -c, 1n];
  const { factors } = factor(R);
  const linear = factors.filter((x) => degree(x.factor) === 1);
  const resolvent = { resolventCubic: R };
  if (!linear.length) return result(4, sq ? "A4" : "S4", { ...info, ...resolvent });
  const nRoots = linear.reduce((s, x) => s + x.multiplicity, 0);
  if (nRoots === 3) return result(4, "V4", { ...info, ...resolvent });
  const [r0, r1] = linear[0].factor; // r1·r + r0 = 0
  const r = reduce(-r0, r1);
  // Kappe–Warren: C4 iff x² − r x + e and x² + b x + (c − r) both split over Q(√D)
  const splits = (u, v) => {
    const disc = subQ(mulQ(u, u), mulQ([4n, 1n], v));
    return disc[0] === 0n || isSquareQ(disc) || isSquareQ(mulQ(disc, D));
  };
  const c4 = splits([-r[0], r[1]], [e, 1n]) && splits([b, 1n], subQ([c, 1n], r));
  return result(4, c4 ? "C4" : "D4", { ...info, ...resolvent });
}

// ---------------------------------------------------------------------------
// Complex interval arithmetic on boxes { re: [lo, hi], im: [lo, hi] }

const Z = [0n, 1n];
const ZI = [Z, Z];
const cAdd = (a, b) => ({ re: ivAdd(a.re, b.re), im: ivAdd(a.im, b.im) });
const cSub = (a, b) => ({ re: ivSub(a.re, b.re), im: ivSub(a.im, b.im) });
const cMul = (a, b) => ({ re: ivSub(ivMul(a.re, b.re), ivMul(a.im, b.im)), im: ivAdd(ivMul(a.re, b.im), ivMul(a.im, b.re)) });
const cConst = (k) => ({ re: [[k, 1n], [k, 1n]], im: ZI });

function enclosures(rs, bits) {
  return rs.map((x) => {
    if (x.rational) return { re: [x.rational, x.rational], im: ZI };
    x.refine(bits);
    if (x.im) return { re: x.re, im: x.im };
    return { re: [x.lo, x.hi], im: ZI };
  });
}

/** Coefficients of Π (y − vₖ) as boxes, constant first. */
function expand(values) {
  let poly = [cConst(1n)];
  for (const v of values) {
    const next = Array(poly.length + 1).fill(null).map(() => cConst(0n));
    poly.forEach((c, i) => {
      next[i + 1] = cAdd(next[i + 1], c);
      next[i] = cSub(next[i], cMul(c, v));
    });
    poly = next;
  }
  return poly;
}

/** The unique integer in a box, undefined while the box is wider than 1/2, null if it holds none. */
function integerIn(b) {
  const [lo, hi] = b.re;
  const width = subQ(hi, lo);
  if (width[0] * 2n >= width[1]) return undefined;
  const [ilo, ihi] = b.im;
  if (ilo[0] > 0n || ihi[0] < 0n) throw new Error("galois: a resolvent coefficient is not real");
  const k = hi[0] >= 0n ? hi[0] / hi[1] : -((-hi[0] + hi[1] - 1n) / hi[1]);
  return k * lo[1] >= lo[0] ? k : null;
}

/**
 * An integer polynomial from the values the function `vals(boxes)` returns
 * for the root boxes, refined until every coefficient is an integer.
 * Returns { poly, boxes } (boxes of the values at the final precision).
 */
function integerResolvent(rs, vals) {
  for (let bits = 32; bits <= 8192; bits *= 2) {
    const boxes = vals(enclosures(rs, bits));
    const coeffs = expand(boxes).map(integerIn);
    if (coeffs.some((c) => c === null)) throw new Error("galois: a resolvent coefficient is not an integer");
    if (coeffs.every((c) => c !== undefined)) return { poly: coeffs, boxes, bits };
  }
  throw new Error("galois: could not pin down the resolvent");
}

const squarefreeInt = (p) => squareFree(primitive(p)).every((s) => s.multiplicity === 1) && degree(primitive(p)) === p.length - 1;

/** The 5-cycles of K5 up to reversal, paired with their complements: six pentagon/pentagram pairs. */
function pentagons() {
  const perms = [];
  const rec = (rest, acc) => {
    if (!rest.length) perms.push([0, ...acc]);
    rest.forEach((x, i) => rec([...rest.slice(0, i), ...rest.slice(i + 1)], [...acc, x]));
  };
  rec([1, 2, 3, 4], []);
  const edges = (cyc) => cyc.map((v, i) => [v, cyc[(i + 1) % 5]].sort().join("")).sort().join(",");
  const seen = new Set();
  const out = [];
  for (const cyc of perms) {
    const E = edges(cyc);
    // the pentagram is the cycle 0, c2, c4, c1, c3
    const star = [cyc[0], cyc[2], cyc[4], cyc[1], cyc[3]];
    const key = [E, edges(star)].sort().join("|");
    if (seen.has(key)) continue;
    seen.add(key);
    out.push(cyc);
  }
  return out; // six cycles
}

function theta(x, cyc) {
  let s = cConst(0n);
  for (let i = 0; i < 5; i += 1) {
    s = cAdd(s, cMul(x[cyc[i]], x[cyc[(i + 1) % 5]])); // pentagon edges
    s = cSub(s, cMul(x[cyc[i]], x[cyc[(i + 2) % 5]])); // pentagram edges
  }
  return cMul(s, s);
}

function phi(x, cyc, a) {
  let s = cConst(0n);
  for (let i = 0; i < 5; i += 1) {
    const u = x[cyc[(a * i) % 5]];
    const v = x[cyc[(a * (i + 1)) % 5]];
    s = cAdd(s, cMul(cMul(u, u), v));
  }
  return s;
}

const hasIntegerRoot = (p) => factor(primitive(p)).factors.some((x) => degree(x.factor) === 1);
function integerRootIndex(p, boxes) {
  const lin = factor(primitive(p)).factors.find((x) => degree(x.factor) === 1);
  if (!lin) return -1;
  const [c0, c1] = lin.factor;
  const r = reduce(-c0, c1);
  const hits = boxes.map((b, k) => (cmpQ(b.re[0], r) <= 0 && cmpQ(r, b.re[1]) <= 0 ? k : -1)).filter((k) => k >= 0);
  return hits.length === 1 ? hits[0] : -2;
}
const cmpQ = ([a, b], [c, d]) => {
  const x = a * d - c * b;
  return x > 0n ? 1 : x < 0n ? -1 : 0;
};

function quintic(f, sq, info) {
  const base = roots(f, { complex: true }).map((r) => r.value);
  const cycles = pentagons();
  // Tschirnhaus y = x² + kx until both resolvents we need are square-free
  for (let k = 0; k < 12; k += 1) {
    const tr = (xs) => (k === 0 ? xs : xs.map((x) => cAdd(cMul(x, x), cMul(cConst(BigInt(k)), x))));
    let six = integerResolvent(base, (xs) => cycles.map((c) => theta(tr(xs), c)));
    if (!squarefreeInt(six.poly)) continue;
    if (!hasIntegerRoot(six.poly)) return result(5, sq ? "A5" : "S5", { ...info, resolvent: six.poly });
    if (!sq) return result(5, "F20", { ...info, resolvent: six.poly });
    // which pentagon: refine until one θ box holds the root
    let idx = integerRootIndex(six.poly, six.boxes);
    for (let bits = six.bits * 2; idx < 0 && bits <= 16384; bits *= 2) {
      const boxes = cycles.map((c) => theta(tr(enclosures(base, bits)), c));
      idx = integerRootIndex(six.poly, boxes);
    }
    if (idx < 0) continue;
    const four = integerResolvent(base, (xs) => [1, 2, 3, 4].map((a) => phi(tr(xs), cycles[idx], a)));
    if (!squarefreeInt(four.poly)) continue;
    return result(5, hasIntegerRoot(four.poly) ? "C5" : "D5", { ...info, resolvent: six.poly });
  }
  throw new Error("galois: no Tschirnhaus transformation made the resolvents square-free");
}


// ---------------------------------------------------------------------------
// The group as explicit permutations of the roots (any degree up to 7)
//
// 1. Orbits on pairs and triples of roots. The polynomial Π (y − (xᵢ + xⱼ))
//    over all pairs has integer coefficients; its irreducible factors are the
//    orbits of the group on pairs (Soicher–McKay). The same for triples.
// 2. H₀ = the permutations that keep every such orbit (and are even when the
//    discriminant is a square). The group G lies in H₀.
// 3. G is conjugate to one transitive group from the table and has exactly
//    those orbits. Usually one subgroup of H₀ qualifies and G is found.
// 4. Otherwise each candidate K is tested, smallest first, with the invariant
//    F_K = Σ_{h ∈ K} Π xₕ₍ᵢ₎ⁱ, which exactly K fixes: G ⊆ K exactly when F_K is a
//    simple integer root of Π (y − F_{cK}) over the cosets cK in H₀.
// Everything is computed from certified root boxes in fixed-point interval
// arithmetic, refined until each integer coefficient is pinned down and the
// values that must differ lie in disjoint boxes.

import { PermGroup, compose, inverse, identity, cycleString, isEven, subsets, allPermutations } from "./permgroup.js";
import { factorModP } from "./factor.js";

// generators (images of 0…n−1) of the transitive groups that are not Aₙ or Sₙ
const TABLE = {
  1: [["1", "1T1", "S"]],
  2: [["C2", "2T1", "S"]],
  3: [["C3", "3T1", "A"], ["S3", "3T2", "S"]],
  4: [
    ["C4", "4T1", [[1, 2, 3, 0]]],
    ["V4", "4T2", [[1, 0, 3, 2], [2, 3, 0, 1]]],
    ["D4", "4T3", [[0, 1, 3, 2], [2, 3, 0, 1]]],
    ["A4", "4T4", "A"], ["S4", "4T5", "S"],
  ],
  5: [
    ["C5", "5T1", [[1, 2, 3, 4, 0]]],
    ["D5", "5T2", [[0, 2, 1, 4, 3], [1, 0, 3, 2, 4]]],
    ["F20", "5T3", [[0, 2, 1, 4, 3], [1, 2, 3, 0, 4]]],
    ["A5", "5T4", "A"], ["S5", "5T5", "S"],
  ],
  6: [
    ["C6", "6T1", [[1, 2, 3, 4, 5, 0]]],
    ["S3", "6T2", [[1, 0, 3, 2, 5, 4], [2, 4, 0, 5, 1, 3]]],
    ["D6", "6T3", [[0, 1, 3, 2, 5, 4], [2, 4, 0, 5, 1, 3]]],
    ["A4", "6T4", [[0, 1, 3, 2, 5, 4], [2, 3, 4, 5, 0, 1]]],
    ["C3×S3", "6T5", [[0, 1, 2, 4, 5, 3], [3, 4, 5, 0, 1, 2]]],
    ["C2×A4", "6T6", [[0, 1, 2, 3, 5, 4], [1, 4, 3, 5, 0, 2]]],
    ["S4", "6T7", [[0, 1, 3, 2, 5, 4], [2, 4, 1, 5, 0, 3]]],
    ["S4", "6T8", [[0, 1, 3, 4, 5, 2], [2, 4, 0, 5, 1, 3]]],
    ["S3×S3", "6T9", [[0, 1, 3, 2, 5, 4], [1, 2, 4, 5, 3, 0]]],
    ["F36", "6T10", [[0, 1, 2, 4, 5, 3], [3, 4, 5, 0, 2, 1]]],
    ["C2×S4", "6T11", [[0, 1, 3, 2, 5, 4], [2, 4, 1, 3, 0, 5]]],
    ["A5", "6T12", [[0, 1, 3, 2, 5, 4], [1, 2, 0, 4, 5, 3]]],
    ["F72", "6T13", [[0, 1, 2, 3, 5, 4], [1, 2, 4, 5, 0, 3]]],
    ["S5", "6T14", [[0, 1, 3, 2, 5, 4], [1, 2, 3, 4, 5, 0]]],
    ["A6", "6T15", "A"], ["S6", "6T16", "S"],
  ],
  7: [
    ["C7", "7T1", [[1, 2, 3, 4, 5, 6, 0]]],
    ["D7", "7T2", [[0, 2, 1, 4, 3, 6, 5], [1, 0, 3, 2, 5, 4, 6]]],
    ["F21", "7T3", [[0, 2, 3, 1, 5, 6, 4], [1, 2, 4, 6, 3, 0, 5]]],
    ["F42", "7T4", [[0, 2, 1, 4, 3, 6, 5], [1, 2, 3, 5, 0, 4, 6]]],
    ["PSL(3,2)", "7T5", [[0, 1, 2, 4, 3, 6, 5], [1, 3, 4, 0, 5, 2, 6]]],
    ["A7", "7T6", "A"], ["S7", "7T7", "S"],
  ],
};
const NOT_SOLVABLE = new Set(["5T4", "5T5", "6T12", "6T14", "6T15", "6T16", "7T5", "7T6", "7T7"]);

const tableGroups = new Map();
/** The transitive groups of degree n: [{ name, label, group, full, signature }]. */
function transitiveGroups(n) {
  if (!tableGroups.has(n)) {
    tableGroups.set(n, TABLE[n].map(([name, label, gens]) => {
      const full = typeof gens === "string" ? gens : null;
      const group = full === "S" ? PermGroup.symmetric(n) : full === "A" ? PermGroup.alternating(n) : new PermGroup(gens, n);
      return { name, label, group, full, even: group.isEven, signature: signature(group) };
    }));
  }
  return tableGroups.get(n);
}
const setSizes = (n) => (n >= 6 ? [2, 3] : n >= 4 ? [2] : []);
const signature = (G) => setSizes(G.degree).map((k) => G.orbitsOnSets(k).map((o) => o.length).sort((a, b) => a - b).join("+")).join("|");

// ---- fixed-point complex boxes [reLo, reHi, imLo, imHi] at scale 2^P

function fixed(P) {
  const one = 1n << P;
  const up = (a) => -((-a) >> P);
  const imul = (a0, a1, b0, b1) => {
    const p = [a0 * b0, a0 * b1, a1 * b0, a1 * b1];
    let lo = p[0];
    let hi = p[0];
    for (let i = 1; i < 4; i += 1) { if (p[i] < lo) lo = p[i]; if (p[i] > hi) hi = p[i]; }
    return [lo >> P, up(hi)];
  };
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]];
  const sub = (a, b) => [a[0] - b[1], a[1] - b[0], a[2] - b[3], a[3] - b[2]];
  const mul = (a, b) => {
    const rr = imul(a[0], a[1], b[0], b[1]);
    const ii = imul(a[2], a[3], b[2], b[3]);
    const ri = imul(a[0], a[1], b[2], b[3]);
    const ir = imul(a[2], a[3], b[0], b[1]);
    return [rr[0] - ii[1], rr[1] - ii[0], ri[0] + ir[0], ri[1] + ir[1]];
  };
  const int = (k) => { const v = BigInt(k) * one; return [v, v, 0n, 0n]; };
  const fromPair = ([lo, hi], a) => [(lo[0] * a * one) / lo[1] - 1n, (hi[0] * a * one) / hi[1] + 1n];
  /** The integer in a box: undefined while the box is too wide, null when it holds none. */
  const integer = (b) => {
    if (b[1] - b[0] >= one >> 1n) return undefined;
    if (b[2] > 0n || b[3] < 0n) return null;
    const k = (b[0] + (one >> 1n)) >> P;
    for (const c of [k, k + 1n]) if (c * one >= b[0] && c * one <= b[1]) return c;
    return null;
  };
  const disjoint = (a, b) => a[1] < b[0] || b[1] < a[0] || a[3] < b[2] || b[3] < a[2];
  const holds = (b, k) => { const v = k * one; return b[0] <= v && v <= b[1] && b[2] <= 0n && b[3] >= 0n; };
  const hasZero = (b) => b[0] <= 0n && b[1] >= 0n && b[2] <= 0n && b[3] >= 0n;
  /** Coefficients of Π (y − vₖ), constant first. */
  const expand = (values) => {
    let poly = [int(1)];
    for (const v of values) {
      const next = Array(poly.length + 1);
      next[poly.length] = poly[poly.length - 1];
      for (let i = poly.length - 1; i >= 0; i -= 1) {
        const t = mul(poly[i], v);
        next[i] = i > 0 ? sub(poly[i - 1], t) : sub(int(0), t);
      }
      poly = next;
    }
    return poly;
  };
  const horner = (p, x) => {
    let acc = int(p[p.length - 1]);
    for (let i = p.length - 2; i >= 0; i -= 1) acc = add(mul(acc, x), int(p[i]));
    return acc;
  };
  return { P, add, sub, mul, int, fromPair, integer, disjoint, holds, hasZero, expand, horner };
}

const MAX_BITS = 1 << 16;
// Tschirnhaus transformations y = c₁x + c₂x² + …, tried in turn; the first is the identity
const SHIFTS = [[1], [1, 1], [2, 1], [1, 1, 1], [1, 0, 1], [2, 1, 1], [1, 2, 1], [1, 1, 1, 1], [1, 0, 1, 1], [1, 2, 3, 1], [1, 1, 1, 1, 1], [2, 1, 3, 1, 1]];
const TRANSFORMS = SHIFTS.length;

/** Certified boxes around the roots of the monic form of f, optionally moved by the k-th Tschirnhaus transformation. */
class RootBoxes {
  constructor(f) {
    this.f = f;
    this.n = f.length - 1;
    this.lead = f[this.n];
    this.values = roots(f, { complex: true }).map((r) => r.value);
    this.cache = new Map();
  }

  at(bits, k = 0) {
    const id = `${bits}:${k}`;
    if (this.cache.has(id)) return this.cache.get(id);
    const F = fixed(BigInt(bits));
    let xs;
    if (k === 0) {
      xs = this.values.map((x) => {
        if (x.rational) { const r = F.fromPair([x.rational, x.rational], this.lead); return [r[0], r[1], 0n, 0n]; }
        x.refine(bits + 8);
        const re = F.fromPair(x.im ? x.re : [x.lo, x.hi], this.lead);
        const im = x.im ? F.fromPair(x.im, this.lead) : [0n, 0n];
        return [re[0], re[1], im[0], im[1]];
      });
    } else {
      // Tschirnhaus: y = c₁x + c₂x² + …, which keeps the action on the numbered roots
      const c = SHIFTS[k];
      xs = this.at(bits, 0).map((x) => {
        let acc = F.int(c[c.length - 1]);
        for (let i = c.length - 2; i >= 0; i -= 1) acc = F.add(F.mul(acc, x), F.int(c[i]));
        return F.mul(acc, x);
      });
    }
    if (this.cache.size > 24) this.cache.clear();
    this.cache.set(id, xs);
    return xs;
  }
}

/**
 * The integer polynomial Π (y − vₖ) for the values `vals(F, boxes)` of the
 * root boxes, once every coefficient is pinned to one integer and the values
 * lie in disjoint boxes. Returns null when the values could not be separated
 * (they are probably equal), so the caller can change the transformation.
 */
function valuePolynomial(R, k, vals, { separate = true } = {}) {
  let found = null;
  for (let bits = 64; bits <= MAX_BITS; bits *= 2) {
    const F = fixed(BigInt(bits));
    const boxes = vals(F, R.at(bits, k));
    if (!found) {
      const coeffs = F.expand(boxes).map(F.integer);
      if (coeffs.some((c) => c === null)) throw new Error("galois: a resolvent coefficient is not an integer");
      if (coeffs.some((c) => c === undefined)) continue;
      found = { poly: coeffs, at: bits };
    }
    let apart = true;
    if (separate) {
      outer: for (let i = 0; i < boxes.length; i += 1) {
        for (let j = i + 1; j < boxes.length; j += 1) {
          if (!F.disjoint(boxes[i], boxes[j])) { apart = false; break outer; }
        }
      }
    }
    if (apart) return { poly: found.poly, boxes, bits, F };
    if (bits >= found.at * 4) return null;
  }
  throw new Error("galois: could not pin down a resolvent");
}


/** The orbits of the Galois group on k-element sets of roots, as lists of index sets. */
function orbitsOnRootSets(R, size) {
  const S = subsets(R.n, size);
  for (let k = 0; k < TRANSFORMS; k += 1) {
    const vals = (F, xs) => S.map((s) => s.slice(1).reduce((acc, i) => F.add(acc, xs[i]), xs[s[0]]));
    const res = valuePolynomial(R, k, vals);
    if (!res) continue;
    const facs = factor(primitive(res.poly)).factors.map((x) => x.factor);
    if (facs.length === 1) return [S];
    for (let bits = res.bits; bits <= MAX_BITS; bits *= 2) {
      const F = bits === res.bits ? res.F : fixed(BigInt(bits));
      const boxes = bits === res.bits ? res.boxes : vals(F, R.at(bits, k));
      const blocks = facs.map(() => []);
      let ok = true;
      for (let s = 0; s < S.length && ok; s += 1) {
        const hit = [];
        facs.forEach((p, j) => { if (F.hasZero(F.horner(p, boxes[s]))) hit.push(j); });
        if (hit.length !== 1) ok = false;
        else blocks[hit[0]].push(S[s]);
      }
      if (ok && blocks.every((b, j) => b.length === degree(facs[j]))) return blocks;
    }
    throw new Error("galois: could not match the roots of a resolvent to its factors");
  }
  throw new Error("galois: no Tschirnhaus transformation separated the resolvent's roots");
}

/** Σ over the permutations in `perms` of Π x_{p(i)}^i, as a box. */
function monomialSum(F, xs, perms) {
  const n = xs.length;
  const pw = xs.map((x) => { const row = [F.int(1), x]; for (let e = 2; e < n; e += 1) row.push(F.mul(row[e - 1], x)); return row; });
  let total = F.int(0);
  for (const p of perms) {
    let m = pw[p[1]][1];
    for (let i = 2; i < n; i += 1) m = F.mul(m, pw[p[i]][i]);
    total = F.add(total, m);
  }
  return total;
}

/** Is the Galois group inside K (a subgroup of H, which is known to contain it)? */
function liesIn(R, H, K) {
  const reps = H.cosets(K);
  if (reps.length === 1) return true;
  const cosets = reps.map((c) => K.elements.map((h) => compose(c, h)));
  for (let k = 0; k < TRANSFORMS; k += 1) {
    const res = valuePolynomial(R, k, (F, xs) => cosets.map((c) => monomialSum(F, xs, c)));
    if (!res) continue;
    // the identity's coset is first (the identity has order 1 and elements are sorted by order)
    const r = res.F.integer(res.boxes[0]);
    if (r === undefined || r === null) return false;
    let v = 0n;
    for (let i = res.poly.length - 1; i >= 0; i -= 1) v = v * r + res.poly[i];
    return v === 0n;
  }
  throw new Error("galois: no Tschirnhaus transformation separated the coset values");
}

function findGroup(R, discSquare) {
  const n = R.n;
  const table = transitiveGroups(n);
  const sizes = setSizes(n);
  const partitions = sizes.map((k) => orbitsOnRootSets(R, k));
  const sig = partitions.map((p) => p.map((o) => o.length).sort((a, b) => a - b).join("+")).join("|");
  const types = table.filter((t) => t.even === discSquare && t.signature === sig);
  if (!types.length) throw new Error("galois: no transitive group matches the orbit lengths");
  if (types.length === 1 && types[0].full) return { type: types[0], group: types[0].group, partitions };

  // H0: keeps every orbit of sets, and the parity
  const blockOf = partitions.map((p) => {
    const m = new Map();
    p.forEach((orbit, b) => orbit.forEach((s) => m.set(s.join(","), b)));
    return m;
  });
  const keeps = (g) => blockOf.every((m) => {
    for (const [s, b] of m) {
      if (m.get(s.split(",").map((x) => g[Number(x)]).sort((x, y) => x - y).join(",")) !== b) return false;
    }
    return true;
  });
  const all = allPermutations(n);
  const h0 = all.filter((g) => (!discSquare || isEven(g)) && keeps(g));
  const H0 = new PermGroup(h0, n).withSmallGenerators();
  const inH0 = new Set(h0.map((g) => g.join(",")));

  const candidates = [];
  const seen = new Set();
  for (const type of types) {
    if (type.full) { candidates.push({ type, group: type.group }); continue; }
    const tried = new Set();
    for (const t of all) {
      const ti = inverse(t);
      const gens = type.group.generators.map((g) => compose(t, compose(g, ti)));
      if (!gens.every((g) => inH0.has(g.join(",")))) continue;
      const gk = gens.map((g) => g.join(",")).join("|");
      if (tried.has(gk)) continue;
      tried.add(gk);
      const K = new PermGroup(gens, n);
      const id = K.elements.map((g) => g.join(",")).sort().join("|");
      if (seen.has(id)) continue;
      seen.add(id);
      // K must have exactly the observed orbits (it is inside H0, so its orbits can only be finer)
      if (sizes.every((k, i) => K.orbitsOnSets(k).length === partitions[i].length)) candidates.push({ type, group: K });
    }
  }
  if (!candidates.length) throw new Error("galois: no candidate group fits the orbits");
  candidates.sort((a, b) => a.group.order - b.group.order);
  for (let i = 0; i < candidates.length - 1; i += 1) {
    if (liesIn(R, H0, candidates[i].group)) return { ...candidates[i], partitions };
  }
  return { ...candidates[candidates.length - 1], partitions };
}

function integerPrimitive(coeffs) {
  const pairs = coeffs.map((c) => (Array.isArray(c) ? reduce(BigInt(c[0]), BigInt(c[1])) : toPair(c)));
  let L = 1n;
  const g = (a, b) => (b ? g(b, a % b) : a);
  for (const [, d] of pairs) L = (L * d) / g(L, d);
  let f = primitive(pairs.map(([n, d]) => n * (L / d)));
  if (f[f.length - 1] < 0n) f = f.map((c) => -c);
  return f;
}

const sub = "₀₁₂₃₄₅₆₇₈₉";
const xName = (i, v = "x") => `${v}${String(i + 1).split("").map((d) => sub[d]).join("")}`;
const sup = (e) => (e === 1 ? "" : String(e).split("").map((d) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[d]).join(""));
const monomialName = (e, v) => e.map((k, i) => (k ? xName(i, v) + sup(k) : "")).join("") || "1";

/** Exponent vectors to try as invariants, simplest first. */
function* exponentVectors(n) {
  const rec = function* (i, left, acc) {
    if (i === n) { if (left === 0) yield acc; return; }
    for (let e = Math.min(left, 3); e >= 0; e -= 1) yield* rec(i + 1, left - e, [...acc, e]);
  };
  for (let total = 1; total <= 4; total += 1) yield* rec(0, total, []);
  yield identity(n);
}

/**
 * The Galois group of an irreducible polynomial of degree ≤ 7 acting on its
 * numbered roots, with the tools of the Galois correspondence.
 */
export class GaloisAction {
  #R;
  #fields = new Map();

  constructor(coeffs) {
    const f = integerPrimitive(coeffs);
    const n = degree(f);
    if (n < 1) throw new RangeError("galoisAction needs a polynomial of degree at least 1");
    if (n > 7) throw new RangeError("Galois groups are computed up to degree 7");
    const { factors } = factor(f);
    if (factors.length !== 1 || factors[0].multiplicity !== 1) throw new RangeError("galoisAction needs an irreducible polynomial; galoisGroup reports the factors of a reducible one");
    const disc = discriminant(f.map((c) => [c, 1n]));
    const sq = isSquareQ([disc.numerator, disc.denominator]);
    this.#R = new RootBoxes(f);
    const found = findGroup(this.#R, sq);
    /** Primitive integer coefficients, constant first. */
    this.polynomial = f;
    this.degree = n;
    /** The roots in the order the permutations use (real roots ascending, then the others). */
    this.roots = this.#R.values;
    /** The group as permutations of the root indices 0 … n−1. */
    this.group = found.group.withSmallGenerators();
    this.name = found.type.name;
    this.label = found.type.label;
    this.order = this.group.order;
    this.solvable = !NOT_SOLVABLE.has(this.label);
    this.discriminant = disc;
    this.discriminantIsSquare = sq;
    /** Orbits of the group on pairs and on triples of roots (triples from degree 6). */
    this.orbitsOnPairs = found.partitions[0] ?? null;
    this.orbitsOnTriples = found.partitions[1] ?? null;
  }

  /** Every element in cycle notation, counted from 1. */
  permutations() {
    return this.group.elements.map(cycleString);
  }

  /** Complex conjugation as an element of the group (the identity when all roots are real). */
  get conjugation() {
    for (let bits = 32; bits <= MAX_BITS; bits *= 2) {
      const F = fixed(BigInt(bits));
      const xs = this.#R.at(bits, 0);
      const p = xs.map((x) => {
        const c = [x[0], x[1], -x[3], -x[2]];
        const hits = xs.map((y, j) => (F.disjoint(c, y) ? -1 : j)).filter((j) => j >= 0);
        return hits.length === 1 ? hits[0] : -1;
      });
      if (p.every((j) => j >= 0)) return p;
    }
    throw new Error("galois: could not pair the roots with their conjugates");
  }

  /**
   * The field fixed by a subgroup H (a PermGroup inside the group, or
   * generators in cycle notation): { degree, polynomial, generator, normal }.
   * `polynomial` is the minimal polynomial (integers, constant first) of the
   * element `generator`, a sum of products of roots that H leaves in place and
   * that generates the fixed field. Its degree is the index of H.
   */
  fixedField(H) {
    const G = this.group;
    if (!(H instanceof PermGroup)) H = G.subgroup(Array.isArray(H) ? H : [H]);
    else if (!G.contains(H)) throw new RangeError("fixedField needs a subgroup of the Galois group");
    const id = H.elements.map((g) => g.join(",")).sort().join("|");
    if (this.#fields.has(id)) return this.#fields.get(id);
    const n = this.degree;
    const reps = G.cosets(H);
    const normal = G.isNormal(H);
    let out = null;
    if (reps.length === 1) out = { degree: 1, polynomial: [0n, 1n], generator: "1", normal };
    for (let k = 0; k < TRANSFORMS && !out; k += 1) {
      const seen = new Set();
      for (const e of exponentVectors(n)) {
        // the H-orbit of the monomial with exponents e
        const orbit = new Map();
        for (const h of H.elements) {
          const v = Array(n).fill(0);
          e.forEach((x, i) => { v[h[i]] = x; });
          orbit.set(v.join(","), v);
        }
        const canon = [...orbit.keys()].sort()[0];
        if (seen.has(canon)) continue;
        seen.add(canon);
        const monos = [...orbit.values()];
        const maxE = Math.max(...e);
        const vals = (F, xs) => {
          const pw = xs.map((x) => { const row = [F.int(1), x]; for (let j = 2; j <= maxE; j += 1) row.push(F.mul(row[j - 1], x)); return row; });
          return reps.map((c) => {
            let total = F.int(0);
            for (const v of monos) {
              let m = null;
              for (let i = 0; i < n; i += 1) {
                if (!v[i]) continue;
                const t = pw[c[i]][v[i]];
                m = m ? F.mul(m, t) : t;
              }
              total = F.add(total, m);
            }
            return total;
          });
        };
        const res = valuePolynomial(this.#R, k, vals);
        if (!res) continue;
        const plain = k === 0 && this.#R.lead === 1n;
        const names = monos.map((v) => monomialName(v, plain ? "x" : "y")).sort();
        out = { degree: reps.length, polynomial: res.poly, generator: names.join(" + "), normal };
        if (!plain) {
          const base = this.#R.lead === 1n ? "x" : `${this.#R.lead}x`;
          const wrap = this.#R.lead === 1n ? "x" : `(${base})`;
          out.substitution = `y = ${SHIFTS[k].map((c, i) => (c ? `${c === 1 ? "" : c}${i === 0 ? base : wrap + sup(i + 1)}` : "")).filter(Boolean).join(" + ")}`;
        }
        break;
      }
    }
    if (!out) throw new Error("galois: no invariant separated the cosets");
    this.#fields.set(id, out);
    return out;
  }

  /**
   * The Galois correspondence: one entry per conjugacy class of subgroups,
   * { subgroup, order, index, conjugates, normal, field, above }, from the
   * trivial subgroup (whose fixed field is the splitting field) to the whole
   * group (fixed field Q). `field` is the fixedField of the subgroup, left
   * null when its degree exceeds `maxDegree` (default 60). `above` lists the
   * entries whose subgroup contains a conjugate of this one.
   */
  correspondence({ maxDegree = 60, maxOrder = 400 } = {}) {
    const classes = this.group.subgroupClasses({ maxOrder });
    const G = this.group;
    const out = classes.map((c) => ({
      subgroup: c.group.withSmallGenerators(),
      order: c.order,
      index: c.index,
      conjugates: c.conjugates,
      normal: c.normal,
      field: c.index <= maxDegree ? this.fixedField(c.group) : null,
      above: [],
    }));
    out.forEach((a, i) => {
      out.forEach((b, j) => {
        if (i === j || b.order <= a.order || b.order % a.order) return;
        if (G.elements.some((g) => b.subgroup.contains(G.conjugate(a.subgroup, g)))) a.above.push(j);
      });
    });
    return out;
  }

  /**
   * Dedekind's theorem in numbers: factor the polynomial modulo the first
   * `count` usable primes and tally the degrees of the factors. Each pattern
   * is the cycle type of an element of the group, and in the long run the
   * patterns appear as often as the group has such elements (Chebotarev).
   * Returns { primes, patterns: [{ type, primes, observed, expected }] }.
   */
  frobenius(count = 200) {
    const f = this.polynomial;
    const tally = new Map();
    let used = 0;
    for (let p = 2; used < count; p += 1) {
      if (!isPrime(p) || f[f.length - 1] % BigInt(p) === 0n) continue;
      let parts;
      try { parts = factorModP(f, p); } catch { continue; }
      const t = parts.map((g) => g.length - 1).sort((a, b) => b - a).join(",");
      if (!tally.has(t)) tally.set(t, []);
      tally.get(t).push(p);
      used += 1;
    }
    const expected = this.group.cycleTypes();
    const patterns = expected.map(({ type, count: c, share }) => {
      const ps = tally.get(type.join(",")) || [];
      tally.delete(type.join(","));
      return { type, elements: c, expected: share, primes: ps.length, observed: ps.length / used, first: ps.slice(0, 5) };
    });
    if (tally.size) throw new Error("galois: a factorization pattern is not a cycle type of the group");
    return { primes: used, patterns };
  }

  /** G ⊇ G′ ⊇ G″ ⊇ … as orders; it reaches 1 exactly when the roots can be written with radicals. */
  derivedSeries() {
    return this.group.derivedSeries().map((H) => H.order);
  }

  toString() {
    return `${this.name} (${this.label}), order ${this.order}, ${this.solvable ? "solvable" : "not solvable"}`;
  }
}

function isPrime(p) {
  if (p < 2) return false;
  for (let d = 2; d * d <= p; d += 1) if (p % d === 0) return false;
  return true;
}

/** The Galois group of an irreducible polynomial of degree ≤ 7 as permutations of its roots. */
export const galoisAction = (coeffs) => new GaloisAction(coeffs);
