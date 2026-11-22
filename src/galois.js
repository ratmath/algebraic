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
  if (n > 5) throw new RangeError("Galois groups are computed up to degree 5");
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

