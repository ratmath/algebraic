// A friendly polynomial toolkit over Q. Polynomials are coefficient arrays,
// constant term first; entries may be numbers, bigints, strings ("3/4"),
// Core Rationals or [n, d] pairs. Results come back as arrays of Core
// Rationals (or BigInt arrays where the text says integer polynomial).

import { Rational } from "@ratmath/core";
import { toPair, reduce, addQ, subQ, mulQ, divQ, detInt, lcm, toRational } from "./q.js";
import { addP, subP, mulP, divmodP, gcdP, trimQ, monicP, evalP, formatP, toIntPoly } from "./qpoly.js";
import { squareFree, countRootsBetween, derivative as dInt, primitive } from "./poly.js";
import { factor } from "./factor.js";

const Z = [0n, 1n];
const ONE = [1n, 1n];
const pairs = (p) => trimQ(p.map((c) => (Array.isArray(c) && typeof c[0] === "bigint" ? reduce(c[0], c[1]) : toPair(c))));
const out = (p) => p.map(toRational);

/** p + q, p − q, p·q. */
export const add = (p, q) => out(addP(pairs(p), pairs(q)));
export const sub = (p, q) => out(subP(pairs(p), pairs(q)));
export const mul = (p, q) => out(mulP(pairs(p), pairs(q)));

/** Quotient and remainder: { quotient, remainder }. */
export function divmod(p, q) {
  const [Q, R] = divmodP(pairs(p), pairs(q));
  return { quotient: out(Q), remainder: out(R) };
}

/** Monic greatest common divisor. */
export const gcd = (p, q) => out(gcdP(pairs(p), pairs(q)));

/** Extended gcd: { gcd, s, t } with s·p + t·q = gcd (monic). */
export function gcdExt(p, q) {
  let [r0, r1] = [pairs(p), pairs(q)];
  let [s0, s1] = [[ONE], []];
  let [t0, t1] = [[], [ONE]];
  while (r1.length) {
    const [Q, R] = divmodP(r0, r1);
    [r0, r1] = [r1, R];
    [s0, s1] = [s1, subP(s0, mulP(Q, s1))];
    [t0, t1] = [t1, subP(t0, mulP(Q, t1))];
  }
  const lc = r0[r0.length - 1] ?? ONE;
  const k = (a) => trimQ(a.map((c) => divQ(c, lc)));
  return { gcd: out(k(r0)), s: out(k(s0)), t: out(k(t0)) };
}

/** p(x) exactly, at a rational x. */
export const evaluate = (p, x) => toRational(evalP(pairs(p), toPair(x)));

/** p(q(x)). */
export function compose(p, q) {
  const P = pairs(p);
  const Qp = pairs(q);
  let acc = [];
  for (let i = P.length - 1; i >= 0; i -= 1) acc = addP(mulP(acc, Qp), [P[i]]);
  return out(acc);
}

export const derivative = (p) => out(pairs(p).slice(1).map((c, i) => mulQ(c, [BigInt(i + 1), 1n])));

/** The monic polynomial with the given roots (rationals), repeated as listed. */
export function fromRoots(rs) {
  return out(rs.reduce((acc, r) => mulP(acc, [subQ(Z, toPair(r)), ONE]), [ONE]));
}

/** Integer primitive polynomial with the same roots (BigInt, constant first). */
export const toInteger = (p) => toIntPoly(pairs(p));

/** Human-readable form, highest degree first. */
export const format = (p, v = "x") => formatP(pairs(p), v);

// ---------------------------------------------------------------------------

function sylvesterInt(p, q) {
  const m = p.length - 1;
  const n = q.length - 1;
  if (m < 0 || n < 0) return 0n;
  if (m === 0 && n === 0) return 1n;
  const size = m + n;
  const hp = [...p].reverse();
  const hq = [...q].reverse();
  const rows = [];
  for (let i = 0; i < n; i += 1) rows.push(Array.from({ length: size }, (_, j) => (j - i >= 0 && j - i <= m ? hp[j - i] : 0n)));
  for (let i = 0; i < m; i += 1) rows.push(Array.from({ length: size }, (_, j) => (j - i >= 0 && j - i <= n ? hq[j - i] : 0n)));
  return detInt(rows);
}

/** Integer coefficients and the scale c with p = ints / c. */
function scaled(p) {
  let L = 1n;
  for (const [, d] of p) L = lcm(L, d);
  return { ints: p.map(([n, d]) => n * (L / d)), L };
}

/** The resultant Res(p, q), exactly: zero iff p and q share a root. */
export function resultant(p, q) {
  const P = pairs(p);
  const Qp = pairs(q);
  const a = scaled(P);
  const b = scaled(Qp);
  const r = sylvesterInt(a.ints, b.ints);
  // Res(p/L, q/M) = Res(p, q) / (L^deg q · M^deg p)
  const den = a.L ** BigInt(Math.max(0, Qp.length - 1)) * b.L ** BigInt(Math.max(0, P.length - 1));
  return new Rational(...reduce(r, den));
}

/** The discriminant: (−1)^(n(n−1)/2)·Res(p, p′)/lc. Zero iff p has a repeated root. */
export function discriminant(p) {
  const P = pairs(p);
  const n = P.length - 1;
  if (n < 1) throw new RangeError("discriminant needs degree ≥ 1");
  const r = resultant(P, derivative(P).map((c) => [c.numerator, c.denominator]));
  const lc = P[n];
  const sgn = (n * (n - 1)) / 2 % 2 ? -1n : 1n;
  return new Rational(...divQ([sgn * r.numerator, r.denominator], lc));
}

/** Number of distinct real roots in (a, b] (rationals or ±Infinity), by a Sturm sequence. */
export function countRoots(p, a = -Infinity, b = Infinity) {
  const ints = primitive(toIntPoly(pairs(p)));
  const end = (x) => (x === Infinity || x === -Infinity ? x : toPair(x));
  let total = 0;
  for (const { factor: f } of squareFree(ints)) total += countRootsBetween(f, end(a), end(b));
  return total;
}

/** The polynomial of least degree through the points [[x, y], …] (Newton's divided differences). */
export function interpolate(points) {
  const xs = points.map(([x]) => toPair(x));
  const c = points.map(([, y]) => toPair(y));
  const n = xs.length;
  for (let j = 1; j < n; j += 1) for (let i = n - 1; i >= j; i -= 1) c[i] = divQ(subQ(c[i], c[i - 1]), subQ(xs[i], xs[i - j]));
  let poly = [c[n - 1]];
  for (let i = n - 2; i >= 0; i -= 1) poly = addP(mulP(poly, [subQ(Z, xs[i]), ONE]), [c[i]]);
  return out(trimQ(poly));
}

/**
 * Partial fractions over Q of num/den:
 * { polynomial, terms: [{ factor, power, numerator }] } with
 * num/den = polynomial + Σ numerator / factorᵖᵒʷᵉʳ, each factor irreducible
 * over Q and deg numerator < deg factor.
 */
export function partialFractions(num, den) {
  const N = pairs(num);
  const D = pairs(den);
  if (!D.length) throw new RangeError("division by the zero polynomial");
  const [poly, rem] = divmodP(N, D);
  const { factors } = factor(toIntPoly(D));
  // den = c · Π fᵢ^eᵢ with fᵢ integer primitive; work with monic-free rational forms
  const blocks = factors.map(({ factor: f, multiplicity: e }) => ({ f: f.map((x) => [x, 1n]), e }));
  const blockPoly = (b) => Array.from({ length: b.e }, () => b.f).reduce((a, f) => mulP(a, f), [ONE]);
  const prodAll = blocks.reduce((a, b) => mulP(a, blockPoly(b)), [ONE]);
  const cst = divQ(D[D.length - 1], prodAll[prodAll.length - 1]);
  let R = trimQ(rem.map((x) => divQ(x, cst))); // rem/den = R / Π blocks
  const terms = [];
  let rest = prodAll;
  for (let k = 0; k < blocks.length; k += 1) {
    const B = blockPoly(blocks[k]);
    const others = divmodP(rest, B)[0];
    // R/(B·others) = A/B + C/others with A·others + C·B = R
    let A;
    if (k === blocks.length - 1) A = R;
    else {
      const { s } = gcdExt(others.map((x) => toRational(x)), B.map((x) => toRational(x)));
      const sP = pairs(s);
      A = divmodP(mulP(R, sP), B)[1];
      R = divmodP(subP(R, mulP(A, others)), B)[0];
    }
    rest = others;
    // A / f^e = Σ A_j / f^j: f-adic expansion of A
    const { f, e } = blocks[k];
    let cur = A;
    const digits = [];
    for (let j = 0; j < e; j += 1) {
      const [Q, Rm] = divmodP(cur, f);
      digits.push(Rm);
      cur = Q;
    }
    // A = Σ digits[j]·f^j  →  A/f^e = Σ digits[j] / f^(e−j)
    digits.forEach((d, j) => {
      if (d.length) terms.push({ factor: f.map(([x]) => x), power: e - j, numerator: out(d) });
    });
  }
  terms.sort((a, b) => a.factor.length - b.factor.length || a.power - b.power);
  return { polynomial: out(poly), terms };
}

export { dInt as derivativeInt };
