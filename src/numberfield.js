// Number fields Q(α) = Q[x] / (m(x)), m irreducible: exact, closed arithmetic
// where every element is a rational polynomial in α of degree < deg m.
// Quadratic fields Q(√d) and the Gaussian rationals Q(i) are the common cases.
//
// A NumberField is also a field object { zero, one, add, sub, mul, div, neg,
// sign } for field-generic algorithms, and each element can be turned into an
// exact AlgebraicReal / ComplexAlgebraic through the chosen embedding of α.

import { Rational } from "@ratmath/core";
import { toPair, reduce, addQ, subQ, mulQ, divQ, pairToString } from "./q.js";
import { addP, subP, mulP, divmodP, monicP, trimQ, degQ, formatP } from "./qpoly.js";
import { factor } from "./factor.js";
import { roots, add as aAdd, mul as aMul, fromRational } from "./arith.js";

const Z = [0n, 1n];
const ONE = [1n, 1n];
const isZq = ([n]) => n === 0n;
const asPair = (x) => (Array.isArray(x) ? reduce(BigInt(x[0]), BigInt(x[1])) : toPair(x));

export class NumberField {
  /**
   * @param minpoly  coefficients, constant first (integers, pairs, Rationals or strings); must be irreducible
   * @param options  { name: how α prints, root: index into roots (real ascending, then complex) or "positive" }
   */
  constructor(minpoly, { name = "α", root = 0, check = true } = {}) {
    const m = trimQ(minpoly.map(asPair));
    if (degQ(m) < 1) throw new RangeError("a number field needs a polynomial of degree ≥ 1");
    if (check) {
      let L = 1n;
      const g = (a, b) => (b ? g(b, a % b) : a);
      for (const [, d] of m) L = (L * d) / g(L, d);
      const ints = m.map(([n, d]) => n * (L / d));
      const { factors } = factor(ints);
      if (factors.length !== 1 || factors[0].multiplicity !== 1) throw new RangeError(`${formatP(m)} is not irreducible over Q`);
    }
    this.m = monicP(m);
    this.n = degQ(this.m);
    this.name = name;
    this.rootChoice = root;
    this.zero = new NumberFieldElement(this, []);
    this.one = new NumberFieldElement(this, [ONE]);
    this.gen = new NumberFieldElement(this, this.n === 1 ? [subQ(Z, this.m[0])] : [Z, ONE]);
    // bound field-interface methods
    this.add = (a, b) => this.of(a).add(b);
    this.sub = (a, b) => this.of(a).sub(b);
    this.mul = (a, b) => this.of(a).mul(b);
    this.div = (a, b) => this.of(a).div(b);
    this.neg = (a) => this.of(a).neg();
    this.inv = (a) => this.of(a).inv();
    this.sign = (a) => (this.of(a).isZero() ? 0 : 1);
  }

  /** Q(√d) for an integer d that is not a square; Q(i) when d = −1. */
  static quadratic(d) {
    d = BigInt(d);
    const name = d === -1n ? "i" : `√${d < 0n ? `(${d})` : d}`;
    return new NumberField([-d, 0n, 1n], { name, root: d > 0n ? 1 : 0 });
  }

  /** The Gaussian rationals Q(i). */
  static gaussian() {
    return NumberField.quadratic(-1);
  }

  /** An element from a rational, a list of coefficients in α (constant first), or an element. */
  of(x) {
    if (x instanceof NumberFieldElement) {
      if (x.field !== this) throw new TypeError("that element belongs to another field");
      return x;
    }
    if (Array.isArray(x) && x.length && typeof x[0] !== "bigint") return this.reduce(x.map(asPair));
    return new NumberFieldElement(this, isZq(asPair(x)) ? [] : [asPair(x)]);
  }

  reduce(coeffs) {
    return new NumberFieldElement(this, trimQ(divmodP(trimQ(coeffs), this.m)[1]));
  }

  /** The embedding of α: an exact AlgebraicReal or ComplexAlgebraic. */
  get root() {
    if (!this._root) {
      if (typeof this.rootChoice === "function") return (this._root = this.rootChoice());
      if (typeof this.rootChoice === "object" && this.rootChoice !== null) return (this._root = this.rootChoice);
      const all = roots(this.m, { complex: true }).map((r) => r.value);
      const upper = (z) => z.im && z.im[0][0] > 0n;
      this._root = this.rootChoice === "positive" ? all.find((r) => r.lo && r.lo[0] > 0n) : this.n === 2 && all.every((z) => z.im) ? all.find(upper) : all[this.rootChoice];
      if (!this._root) throw new RangeError("no such root of the defining polynomial");
    }
    return this._root;
  }

  toString() {
    return `Q(${this.name}), ${this.name} a root of ${formatP(this.m)}`;
  }
}

export class NumberFieldElement {
  constructor(field, coeffs) {
    this.field = field;
    this.c = coeffs; // pairs, constant first, degree < n, trimmed
  }

  get F() {
    return this.field;
  }

  add(b) {
    return new NumberFieldElement(this.field, trimQ(addP(this.c, this.field.of(b).c)));
  }
  sub(b) {
    return new NumberFieldElement(this.field, trimQ(subP(this.c, this.field.of(b).c)));
  }
  neg() {
    return new NumberFieldElement(this.field, this.c.map(([n, d]) => [-n, d]));
  }
  mul(b) {
    return this.field.reduce(mulP(this.c, this.field.of(b).c));
  }
  isZero() {
    return this.c.length === 0;
  }
  equals(b) {
    return this.sub(b).isZero();
  }
  /** Inverse by the extended Euclidean algorithm in Q[x]. */
  inv() {
    if (this.isZero()) throw new RangeError("division by zero");
    let [r0, r1] = [this.field.m, this.c];
    let [t0, t1] = [[], [ONE]];
    while (degQ(r1) > 0) {
      const [q, r] = divmodP(r0, r1);
      [r0, r1] = [r1, r];
      [t0, t1] = [t1, trimQ(subP(t0, mulP(q, t1)))];
    }
    const k = trimQ(r1)[0];
    return this.field.reduce(t1.map((c) => divQ(c, k)));
  }
  div(b) {
    return this.mul(this.field.of(b).inv());
  }
  pow(k) {
    if (k < 0) return this.inv().pow(-k);
    let out = this.field.one;
    let base = this;
    while (k > 0) {
      if (k & 1) out = out.mul(base);
      k >>= 1;
      if (k) base = base.mul(base);
    }
    return out;
  }

  /** The matrix of multiplication by this element in the basis 1, α, …, α^(n−1). */
  #matrix() {
    const n = this.field.n;
    const cols = [];
    let e = this;
    for (let j = 0; j < n; j += 1) {
      cols.push(Array.from({ length: n }, (_, i) => e.c[i] ?? Z));
      e = e.mul(this.field.gen);
    }
    return Array.from({ length: n }, (_, i) => cols.map((c) => c[i]));
  }
  /** The field norm (product of all conjugates), a rational. */
  norm() {
    const A = this.#matrix();
    const n = A.length;
    let det = ONE;
    for (let c = 0; c < n; c += 1) {
      const p = A.findIndex((r, i) => i >= c && !isZq(r[c]));
      if (p < 0) return new Rational(0n, 1n);
      if (p !== c) {
        [A[c], A[p]] = [A[p], A[c]];
        det = [-det[0], det[1]];
      }
      det = mulQ(det, A[c][c]);
      for (let i = c + 1; i < n; i += 1) {
        const f = divQ(A[i][c], A[c][c]);
        A[i] = A[i].map((x, j) => subQ(x, mulQ(f, A[c][j])));
      }
    }
    return new Rational(...det);
  }
  /** The field trace (sum of all conjugates), a rational. */
  trace() {
    const A = this.#matrix();
    return new Rational(...A.reduce((s, r, i) => addQ(s, r[i]), Z));
  }

  /** Coefficients in the basis 1, α, … as Core Rationals (length n). */
  coefficients() {
    return Array.from({ length: this.field.n }, (_, i) => new Rational(...(this.c[i] ?? Z)));
  }
  /** The rational value when this element lies in Q, else null. */
  toRational() {
    return this.c.length <= 1 ? new Rational(...(this.c[0] ?? Z)) : null;
  }

  /** This element as an exact AlgebraicReal or ComplexAlgebraic, through the field's embedding. */
  toAlgebraic() {
    if (this.c.length <= 1) return fromRational(this.c[0] ?? Z);
    const a = this.field.root;
    let out = fromRational(this.c[this.c.length - 1]);
    for (let i = this.c.length - 2; i >= 0; i -= 1) out = aAdd(aMul(out, a), fromRational(this.c[i]));
    return out;
  }

  toString() {
    const v = this.field.name;
    if (!this.c.length) return "0";
    const parts = [];
    for (let i = 0; i < this.c.length; i += 1) {
      const c = this.c[i];
      if (isZq(c)) continue;
      const neg = c[0] < 0n;
      const a = neg ? [-c[0], c[1]] : c;
      const one = a[0] === 1n && a[1] === 1n;
      const coef = pairToString(a);
      const sup = String(i).split("").map((d) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[d]).join("");
      const term = i === 0 ? coef : `${one ? "" : a[1] === 1n ? coef : `(${coef})`}${v}${i > 1 ? sup : ""}`;
      parts.push((parts.length ? (neg ? " − " : " + ") : neg ? "−" : "") + term);
    }
    return parts.join("");
  }
}
