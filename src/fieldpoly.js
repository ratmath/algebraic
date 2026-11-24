// Polynomials over a number field K = Q(α): gcd, factoring (Trager's method),
// roots in K, extensions K(β), composite fields and splitting fields.
//
// Trager: for a square-free g in K[x], the norm N(x) = Π_σ g^σ(x − s·σ(α)) is a
// rational polynomial. For all but finitely many integers s it is square-free,
// and then its irreducible factors N_i over Q match the irreducible factors of
// g over K one to one:  g_i(x) = gcd(g(x), N_i(x + s·α)).
// When g is irreducible, N is the minimal polynomial of β + s·α for a root β
// of g: a primitive element of K(β), which gives the extension as a plain
// NumberField again.

import { toPair, reduce } from "./q.js";
import { toIntPoly, trimQ } from "./qpoly.js";
import { factor } from "./factor.js";
import { interpolate } from "./polys.js";
import { NumberField, NumberFieldElement } from "./numberfield.js";
import { add as aAdd, mul as aMul, fromRational, minimalPolynomial, equals as aEquals, roots as aRoots } from "./arith.js";

const MAX_ABSOLUTE_DEGREE = 48;

// ---- polynomials with coefficients in K (arrays of elements, constant first)

const trim = (p) => {
  let n = p.length;
  while (n > 0 && p[n - 1].isZero()) n -= 1;
  return n === p.length ? p : p.slice(0, n);
};
const deg = (p) => trim(p).length - 1;
const lead = (p) => p[p.length - 1];
const monic = (p) => {
  p = trim(p);
  if (!p.length) return p;
  const inv = lead(p).inv();
  return p.map((c) => c.mul(inv));
};
function divmod(K, a, b) {
  a = trim(a).slice();
  b = trim(b);
  if (!b.length) throw new RangeError("division by the zero polynomial");
  const q = Array(Math.max(0, a.length - b.length + 1)).fill(K.zero);
  const inv = lead(b).inv();
  for (let i = a.length - b.length; i >= 0; i -= 1) {
    const f = a[i + b.length - 1].mul(inv);
    q[i] = f;
    if (f.isZero()) continue;
    for (let j = 0; j < b.length; j += 1) a[i + j] = a[i + j].sub(f.mul(b[j]));
  }
  return [trim(q), trim(a.slice(0, b.length - 1))];
}
function gcd(K, a, b) {
  a = trim(a);
  b = trim(b);
  while (b.length) [a, b] = [b, divmod(K, a, b)[1]];
  return monic(a);
}
const derivative = (K, p) => trim(p.slice(1).map((c, i) => c.mul(K.of(BigInt(i + 1)))));
const evalAt = (K, p, x) => p.reduceRight((acc, c) => acc.mul(x).add(c), K.zero);
/** p(x + t) for an element t. */
function shift(K, p, t) {
  const g = p.slice();
  const n = g.length - 1;
  for (let i = 0; i < n; i += 1) for (let j = n - 1; j >= i; j -= 1) g[j] = g[j].add(g[j + 1].mul(t));
  return g;
}

function lift(K, coeffs) {
  return trim(coeffs.map((c) => K.of(c)));
}

const intPoly = (rationals) => toIntPoly(trimQ(rationals.map((r) => toPair(r))));

/** Norm of the square-free monic g(x − s·α): a primitive integer polynomial of degree n·deg g. */
function normPoly(K, g, s) {
  const D = K.n * deg(g);
  const sa = K.gen.mul(K.of(BigInt(s)));
  const pts = [];
  for (let i = 0; i <= D; i += 1) {
    const x0 = BigInt(i - (D >> 1));
    pts.push([x0, evalAt(K, g, K.of(x0).sub(sa)).norm()]);
  }
  return intPoly(interpolate(pts));
}

/** Irreducible factors of the square-free monic g over K, with the s and norm factors used. */
function factorSquareFree(K, g) {
  if (deg(g) === 1) return { parts: [g], s: 0, norms: null };
  for (let s = 0; s < 200; s += 1) {
    const N = normPoly(K, g, s);
    const { factors } = factor(N);
    if (factors.some((f) => f.multiplicity > 1)) continue; // not square-free: try the next shift
    if (factors.length === 1) return { parts: [g], s, norms: [factors[0].factor] };
    const sa = K.gen.mul(K.of(BigInt(s)));
    const parts = [];
    const norms = [];
    for (const { factor: f } of factors) {
      const h = gcd(K, g, shift(K, lift(K, f), sa));
      parts.push(h);
      norms.push(f);
    }
    return { parts, s, norms };
  }
  throw new Error("number field factoring: no square-free norm found");
}

/**
 * Factor a polynomial with coefficients in K into irreducibles over K.
 * Returns { unit, factors: [{ factor, multiplicity }] } with monic factors
 * (arrays of elements, constant first), lowest degree first.
 */
export function factorOver(K, coeffs) {
  const p = lift(K, coeffs);
  if (!p.length) throw new RangeError("cannot factor the zero polynomial");
  const unit = lead(p);
  let f = monic(p);
  const out = [];
  if (deg(f) >= 1) {
    const sqf = divmod(K, f, gcd(K, f, derivative(K, f)))[0];
    for (const h of factorSquareFree(K, monic(sqf)).parts) {
      let multiplicity = 0;
      for (;;) {
        const [q, r] = divmod(K, f, h);
        if (r.length) break;
        f = q;
        multiplicity += 1;
      }
      out.push({ factor: h, multiplicity });
    }
  }
  out.sort((a, b) => a.factor.length - b.factor.length);
  return { unit, factors: out };
}

/** The roots of a polynomial that lie in K: [{ value, multiplicity }]. */
export function rootsIn(K, coeffs) {
  return factorOver(K, coeffs).factors.filter((f) => f.factor.length === 2).map((f) => ({ value: f.factor[0].neg(), multiplicity: f.multiplicity }));
}

/**
 * K(β) for a root β of g, irreducible over K, as a NumberField L = Q(γ) with
 * γ = β + s·α. Returns { field, embed, root }: `embed` maps elements of K into
 * L and `root` is β in L. `pick` chooses the embedding of γ when it matters
 * (an exact algebraic number); by default the first root of γ's polynomial.
 */
export function extend(K, coeffs, { name = "γ", pick = null } = {}) {
  const g = monic(lift(K, coeffs));
  if (deg(g) < 1) throw new RangeError("extend needs a polynomial of degree ≥ 1");
  if (deg(g) * K.n > MAX_ABSOLUTE_DEGREE) throw new RangeError(`the extension has degree ${deg(g) * K.n} over Q; the limit is ${MAX_ABSOLUTE_DEGREE}`);
  const { parts, s, norms } = factorSquareFree(K, g);
  if (parts.length !== 1) throw new RangeError("extend needs a polynomial that is irreducible over the field (use factor first)");
  if (deg(g) === 1) return { field: K, embed: (x) => K.of(x), root: g[0].neg() };
  const m = norms[0];
  const L = new NumberField(m, { name, root: pick ?? 0, check: false });
  const autoPick = !pick && K.n > 1;
  // α in L: the common root of m_K(y) and g(γ − s·y) (coefficients of g read in y)
  let alpha;
  if (K.n === 1) alpha = L.of(K.gen.c[0] ?? [0n, 1n]);
  else {
    const mK = K.m.map((c) => L.of(c));
    // G(y) = Σ_j g_j(y)·(γ − s·y)^j in L[y]
    const lin = [L.gen, L.of(BigInt(-s))];
    const mulL = (a, b) => {
      const o = Array(a.length + b.length - 1).fill(L.zero);
      a.forEach((x, i) => b.forEach((y, j) => { o[i + j] = o[i + j].add(x.mul(y)); }));
      return o;
    };
    let G = [L.zero];
    let pw = [L.one];
    for (let j = 0; j < g.length; j += 1) {
      const cj = g[j].c.map((c) => L.of(c)); // polynomial in y
      if (cj.length) {
        const term = mulL(cj, pw);
        const o = Array(Math.max(G.length, term.length)).fill(L.zero);
        G.forEach((x, i) => { o[i] = o[i].add(x); });
        term.forEach((x, i) => { o[i] = o[i].add(x); });
        G = o;
      }
      pw = mulL(pw, lin);
    }
    const h = gcd(L, mK, G);
    if (h.length !== 2) throw new Error("number field extension: could not express the old generator in the new field");
    alpha = h[0].neg();
  }
  const embed = (x) => {
    const e = K.of(x);
    let acc = L.zero;
    for (let i = e.c.length - 1; i >= 0; i -= 1) acc = acc.mul(alpha).add(L.of(e.c[i]));
    return acc;
  };
  if (autoPick) {
    // embed γ so that α keeps the value it has in K
    L.rootChoice = () => {
      for (const { value: r } of aRoots(m, { complex: true })) {
        let v = fromRational([0n, 1n]);
        for (let i = alpha.c.length - 1; i >= 0; i -= 1) v = aAdd(aMul(v, r), fromRational(alpha.c[i]));
        if (aEquals(v, K.root)) return r;
      }
      throw new Error("number field extension: no embedding agrees with the base field");
    };
  }
  const root = L.gen.sub(alpha.mul(L.of(BigInt(s))));
  return { field: L, embed, root };
}

/** The value of a NumberFieldElement of L under a chosen exact root. */
const valueOf = (x) => x.toAlgebraic();

/**
 * One field holding several exact algebraic numbers: { field, elements } with
 * elements[i] the i-th number as an element of `field`, whose generator is
 * embedded so that every element's toAlgebraic() equals the number given.
 */
export function fieldContaining(numbers, { name = "γ" } = {}) {
  let K = new NumberField([0n, 1n], { name, check: false }); // Q
  let elements = [];
  for (const a of numbers) {
    const r = a.rational;
    if (r) {
      elements.push(K.of(reduce(r[0], r[1])));
      continue;
    }
    const mp = minimalPolynomial(a);
    // the factor of a's minimal polynomial over K that a satisfies
    const { factors } = factorOver(K, mp);
    let done = false;
    for (const { factor: g } of factors) {
      if (g.length === 2) {
        if (aEquals(valueOf(g[0].neg()), a)) {
          elements.push(g[0].neg());
          done = true;
          break;
        }
        continue;
      }
      // is a a root of g? evaluate with exact numbers
      let v = fromRational([0n, 1n]);
      for (let i = g.length - 1; i >= 0; i -= 1) v = aAdd(aMul(v, a), valueOf(g[i]));
      if (!(v.rational && v.rational[0] === 0n)) continue;
      // adjoin it, with γ = a + s·α embedded accordingly
      const { parts, s } = factorSquareFree(K, g);
      if (parts.length !== 1) continue;
      const gamma = K.n === 1 ? a : aAdd(a, aMul(fromRational([BigInt(s), 1n]), K.root));
      const ext = extend(K, g, { name, pick: gamma });
      elements = elements.map((e) => ext.embed(e));
      elements.push(ext.root);
      K = ext.field;
      done = true;
      break;
    }
    if (!done) throw new Error("fieldContaining: could not place a number in the field");
  }
  return { field: K, elements };
}

/**
 * The splitting field of a rational polynomial: { field, roots, degree } with
 * every root as an element of `field` (repeated roots listed once).
 */
export function splittingField(coeffs, { name = "γ" } = {}) {
  let K = new NumberField([0n, 1n], { name, check: false });
  let poly = lift(K, coeffs.map((c) => toPair(c)));
  let found = [];
  for (;;) {
    const { factors } = factorOver(K, poly);
    found = factors.filter((f) => f.factor.length === 2).map((f) => f.factor[0].neg());
    const big = factors.find((f) => f.factor.length > 2);
    if (!big) return { field: K, roots: found, degree: K.n };
    const ext = extend(K, big.factor, { name });
    poly = trim(poly).map((c) => ext.embed(c));
    K = ext.field;
  }
}

NumberField.prototype.factor = function (coeffs) {
  return factorOver(this, coeffs);
};
NumberField.prototype.roots = function (coeffs) {
  return rootsIn(this, coeffs);
};
NumberField.prototype.extend = function (coeffs, opts) {
  return extend(this, coeffs, opts);
};
NumberField.containing = fieldContaining;
NumberField.splitting = splittingField;

