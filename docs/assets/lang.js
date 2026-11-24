// A small exact expression language for the docs: numbers are exact
// algebraic numbers, x makes polynomials, and the library's functions are
// available by short names. Evaluation returns plain result descriptors
// ({ type: "number" | "poly" | "list" | "text" | "record", … }) that the
// page renders; nothing here touches the DOM.
//
//   sqrt(2) + sqrt(3)        cospi(1/5)        roots(x^5 - x - 1)
//   (1 + i)^(1/3)            galois(x^4 - 2)   cf(2^(1/3), 15)

import * as A from "./algebraic.js?v=1ae5f5eda2";

const { qpoly } = A;

// ---------------------------------------------------------------------------
// Tokens and parsing

const TOKEN = /\s*(?:(\d+(?:\.\d+)?|\.\d+)|("[^"]*")|([A-Za-z_][A-Za-z_0-9]*)|(√|\*\*|[-+*/^(),=;·×]))/y;

function tokenize(src) {
  const out = [];
  TOKEN.lastIndex = 0;
  let pos = 0;
  while (pos < src.length) {
    if (/^\s*$/.test(src.slice(pos))) break;
    TOKEN.lastIndex = pos;
    const m = TOKEN.exec(src);
    if (!m) throw new SyntaxError(`I can't read "${src.slice(pos).trim().slice(0, 12)}"`);
    pos = TOKEN.lastIndex;
    if (m[1]) out.push({ t: "num", v: m[1] });
    else if (m[2]) out.push({ t: "str", v: m[2].slice(1, -1) });
    else if (m[3]) out.push({ t: "id", v: m[3] });
    else out.push({ t: "op", v: m[4] === "**" ? "^" : m[4] === "·" || m[4] === "×" ? "*" : m[4] });
  }
  out.push({ t: "end" });
  return out;
}

function parse(src) {
  const toks = tokenize(src);
  let k = 0;
  const peek = () => toks[k];
  const take = () => toks[k++];
  const isOp = (v) => peek().t === "op" && peek().v === v;
  const expect = (v) => {
    if (!isOp(v)) throw new SyntaxError(`expected "${v}"`);
    take();
  };
  const startsPrimary = () => {
    const t = peek();
    return t.t === "num" || t.t === "id" || t.t === "str" || (t.t === "op" && (t.v === "(" || t.v === "√"));
  };
  function primary() {
    const t = take();
    if (t.t === "num") return { k: "num", v: t.v };
    if (t.t === "str") return { k: "str", v: t.v };
    if (t.t === "op" && t.v === "(") {
      const e = expr();
      expect(")");
      return e;
    }
    if (t.t === "op" && t.v === "√") return { k: "call", name: "sqrt", args: [power()] };
    if (t.t === "id") {
      if (isOp("(")) {
        take();
        const args = [];
        if (!isOp(")")) {
          args.push(expr());
          while (isOp(",")) {
            take();
            args.push(expr());
          }
        }
        expect(")");
        return { k: "call", name: t.v, args };
      }
      return { k: "var", name: t.v };
    }
    throw new SyntaxError(t.t === "end" ? "the expression ends too soon" : `unexpected "${t.v}"`);
  }
  function power() {
    const base = primary();
    if (isOp("^")) {
      take();
      return { k: "^", a: base, b: unary() };
    }
    return base;
  }
  function unary() {
    if (isOp("-")) {
      take();
      return { k: "neg", a: unary() };
    }
    if (isOp("+")) {
      take();
      return unary();
    }
    return power();
  }
  function term() {
    let a = unary();
    for (;;) {
      if (isOp("*") || isOp("/")) {
        const op = take().v;
        a = { k: op, a, b: unary() };
      } else if (startsPrimary()) {
        a = { k: "*", a, b: power() }; // implicit: 2x, 3i, 2sqrt(2)
      } else return a;
    }
  }
  function expr() {
    let a = term();
    while (isOp("+") || isOp("-")) {
      const op = take().v;
      a = { k: op, a, b: term() };
    }
    return a;
  }
  function statement() {
    if (peek().t === "id" && toks[k + 1].t === "op" && toks[k + 1].v === "=") {
      const name = take().v;
      take();
      return { k: "let", name, e: expr() };
    }
    return expr();
  }
  const s = statement();
  if (peek().t !== "end") throw new SyntaxError(`unexpected "${peek().v}"`);
  return s;
}

// ---------------------------------------------------------------------------
// Values

const isNum = (v) => v instanceof A.AlgebraicReal || v instanceof A.ComplexAlgebraic;
const isPoly = (v) => v && v.poly && Array.isArray(v.c);
const P = (c) => ({ poly: true, c: qpoly.trimQ(c) });
const ONE = [1n, 1n];

function rationalOf(v) {
  if (isNum(v)) return A.rationalOf(v);
  if (isPoly(v) && v.c.length <= 1) return v.c[0] ?? [0n, 1n];
  return null;
}
function asPoly(v) {
  if (isPoly(v)) return v;
  const r = rationalOf(v);
  if (!r) throw new TypeError("polynomial coefficients must be rational here");
  return P([r]);
}
function asNum(v) {
  if (isNum(v)) return v;
  if (isPoly(v) && v.c.length <= 1) return A.fromRational(v.c[0] ?? [0n, 1n]);
  throw new TypeError("expected a number, got a polynomial");
}
function asInt(v, what = "an integer") {
  const r = rationalOf(v);
  if (!r || r[1] !== 1n) throw new TypeError(`expected ${what}`);
  return r[0];
}
function asRat(v) {
  const r = rationalOf(v);
  if (!r) throw new TypeError("expected a rational number");
  return r;
}
/** A polynomial argument as integer-or-rational coefficient pairs. */
const coeffs = (v) => asPoly(v).c;

function binary(op, a, b) {
  if (isPoly(a) || isPoly(b)) {
    if (op === "/" && !isPoly(b)) {
      const r = asRat(b);
      if (r[0] === 0n) throw new RangeError("division by zero");
      return P(qpoly.scaleP(asPoly(a).c, [r[1], r[0]]));
    }
    const pa = asPoly(a).c;
    const pb = asPoly(b).c;
    if (op === "+") return P(qpoly.addP(pa, pb));
    if (op === "-") return P(qpoly.subP(pa, pb));
    if (op === "*") return P(qpoly.mulP(pa, pb));
    if (op === "/") {
      const [q, r] = qpoly.divmodP(pa, pb);
      if (r.length) throw new RangeError("the division leaves a remainder; use pf(p, q) for partial fractions");
      return P(q);
    }
  }
  const x = asNum(a);
  const y = asNum(b);
  if (op === "+") return A.add(x, y);
  if (op === "-") return A.sub(x, y);
  if (op === "*") return A.mul(x, y);
  return A.div(x, y);
}

function power(a, b) {
  const e = asRat(b);
  if (isPoly(a)) {
    if (e[1] !== 1n || e[0] < 0n) throw new TypeError("polynomials take whole-number powers");
    let out = [ONE];
    for (let k = 0n; k < e[0]; k += 1n) out = qpoly.mulP(out, a.c);
    return P(out);
  }
  if (e[0] > 2000n || e[0] < -2000n) throw new RangeError("that power is too large");
  const base = A.pow(a, Number(e[0]));
  return e[1] === 1n ? base : A.nthRoot(base, Number(e[1]));
}

// ---------------------------------------------------------------------------
// Results

export const num = (value, note) => ({ type: "number", value, note });
const poly = (c) => ({ type: "poly", c });
const text = (s) => ({ type: "text", text: s });
const list = (items, label) => ({ type: "list", items, label });
const record = (title, rows) => ({ type: "record", title, rows });

const fmtP = (c) => qpoly.formatP(c);
const fmtInts = (p) => qpoly.formatP(p.map((c) => [c, 1n]));
const toResult = (v) => (isPoly(v) ? (v.c.length <= 1 ? num(A.fromRational(v.c[0] ?? [0n, 1n])) : poly(v.c)) : isNum(v) ? num(v) : v);

// ---------------------------------------------------------------------------
// Functions

const F = {
  // numbers
  sqrt: { args: 1, doc: "principal square root", f: ([a]) => A.sqrt(asNum(a)) },
  cbrt: { args: 1, doc: "real cube root of a real, principal for complex", f: ([a]) => A.nthRoot(asNum(a), 3) },
  root: { args: 2, doc: "root(a, n): n-th root", f: ([a, n]) => A.nthRoot(asNum(a), Number(asInt(n))) },
  abs: { args: 1, doc: "absolute value or modulus", f: ([a]) => A.abs(asNum(a)) },
  norm: { args: 1, doc: "|a|²", f: ([a]) => A.norm(asNum(a)) },
  re: { args: 1, doc: "real part", f: ([a]) => A.re(asNum(a)) },
  im: { args: 1, doc: "imaginary part", f: ([a]) => A.im(asNum(a)) },
  conj: { args: 1, doc: "complex conjugate", f: ([a]) => A.conjugate(asNum(a)) },
  minpoly: { args: 1, doc: "minimal polynomial", f: ([a]) => P(A.minimalPolynomial(asNum(a)).map((c) => [c, 1n])) },
  degree: { args: 1, doc: "degree of a number or polynomial", f: ([a]) => A.fromRational(BigInt(isPoly(a) ? a.c.length - 1 : A.minimalPolynomial(asNum(a)).length - 1)) },
  digits: {
    args: [1, 2], doc: "digits(a, n): n certified digits",
    f: ([a, n]) => text(asNum(a).toString(n ? Number(asInt(n)) : 40)),
  },
  // angles and unity
  cospi: { args: 1, doc: "cospi(r) = cos(πr), exactly", f: ([r]) => A.cosPi(asRat(r)) },
  sinpi: { args: 1, doc: "sinpi(r) = sin(πr)", f: ([r]) => A.sinPi(asRat(r)) },
  tanpi: { args: 1, doc: "tanpi(r) = tan(πr)", f: ([r]) => A.tanPi(asRat(r)) },
  unity: { args: [1, 2], doc: "unity(n, k): e^(2πik/n)", f: ([n, k]) => A.rootOfUnity(Number(asInt(n)), k ? Number(asInt(k)) : 1) },
  cyclotomic: { args: 1, doc: "the n-th cyclotomic polynomial", f: ([n]) => P(A.cyclotomic(Number(asInt(n))).map((c) => [c, 1n])) },
  chebyshev: { args: 1, doc: "Chebyshev polynomial Tₙ", f: ([n]) => P(A.chebyshevT(Number(asInt(n))).map((c) => [c, 1n])) },
  // polynomials
  roots: {
    args: 1, doc: "real roots of a polynomial",
    f: ([p]) => list(A.roots(coeffs(p)).map((r) => num(r.value, r.multiplicity > 1 ? `multiplicity ${r.multiplicity}` : null)), "real roots"),
  },
  croots: {
    args: 1, doc: "all roots, complex too",
    f: ([p]) => list(A.roots(coeffs(p), { complex: true }).map((r) => num(r.value, r.multiplicity > 1 ? `multiplicity ${r.multiplicity}` : null)), "roots"),
  },
  factor: {
    args: 1, doc: "factor over the integers (or an integer into primes)",
    f: ([p]) => {
      if (isNum(p) || (isPoly(p) && p.c.length <= 1)) {
        const n = asInt(p, "an integer or a polynomial");
        const fs = A.factorInteger(n);
        return text((n < 0n ? "−" : "") + (fs.map(({ prime, exponent }) => (exponent > 1 ? `${prime}^${exponent}` : `${prime}`)).join(" · ") || String(n)));
      }
      const ints = A.polys.toInteger(coeffs(p));
      const { content, factors } = A.factor(ints);
      // keep the overall constant: p = (lc ratio) · content · Π
      const lcIn = coeffs(p)[coeffs(p).length - 1];
      const c = qpoly.scaleP(qpoly.scaleP([[content, 1n]], lcIn), [1n, ints[ints.length - 1]])[0];
      const k = c[0] === 1n && c[1] === 1n ? "" : `${c[1] === 1n ? c[0] : `(${c[0]}/${c[1]})`} · `;
      return text(k + factors.map(({ factor: f, multiplicity }) => `(${fmtInts(f)})${multiplicity > 1 ? `^${multiplicity}` : ""}`).join(" · "));
    },
  },
  irreducible: { args: 1, doc: "is it irreducible over Q?", f: ([p]) => text(A.isIrreducible(A.polys.toInteger(coeffs(p))) ? "irreducible over Q" : "reducible over Q") },
  disc: { args: 1, doc: "discriminant", f: ([p]) => A.fromRational(A.discriminant(coeffs(p))) },
  res: { args: 2, doc: "resultant", f: ([p, q]) => A.fromRational(A.resultant(coeffs(p), coeffs(q))) },
  gcd: {
    args: 2, doc: "gcd of polynomials or integers",
    f: ([p, q]) => (isPoly(p) && p.c.length > 1) || (isPoly(q) && q.c.length > 1)
      ? P(A.polys.gcd(coeffs(p), coeffs(q)).map((r) => [r.numerator, r.denominator]))
      : A.fromRational(A.nt.egcd(asInt(p), asInt(q)).g),
  },
  deriv: { args: 1, doc: "derivative", f: ([p]) => P(A.polys.derivative(coeffs(p)).map((r) => [r.numerator, r.denominator])) },
  at: {
    args: 2, doc: "at(p, a): p evaluated at a number, exactly",
    f: ([p, a]) => {
      const c = coeffs(p);
      const x = asNum(a);
      let acc = A.fromRational([0n, 1n]);
      for (let i = c.length - 1; i >= 0; i -= 1) acc = A.add(A.mul(acc, x), A.fromRational(c[i]));
      return acc;
    },
  },
  countroots: {
    args: [1, 3], doc: "countroots(p, a, b): distinct real roots in (a, b]",
    f: ([p, a, b]) => A.fromRational(BigInt(A.countRoots(coeffs(p), a ? asRat(a) : -Infinity, b ? asRat(b) : Infinity))),
  },
  pf: {
    args: 2, doc: "pf(p, q): partial fractions of p/q",
    f: ([p, q]) => {
      const r = A.partialFractions(coeffs(p), coeffs(q));
      const pieces = [];
      const pol = r.polynomial.map((x) => [x.numerator, x.denominator]);
      if (qpoly.trimQ(pol).length) pieces.push(fmtP(pol));
      for (const t of r.terms) {
        const n = fmtP(t.numerator.map((x) => [x.numerator, x.denominator]));
        const fs = fmtInts(t.factor);
        const d = `${fs.includes(" ") ? `(${fs})` : fs}${t.power > 1 ? `^${t.power}` : ""}`;
        pieces.push(`${/[ /]/.test(n) ? `(${n})` : n}/${d}`);
      }
      return text(pieces.join(" + ").replace(/\+ −/g, "− ").replace(/\+ \(−/g, "− ("));
    },
  },
  galois: {
    args: 1, doc: "Galois group (degree ≤ 5)",
    f: ([p]) => {
      const g = A.galoisGroup(coeffs(p));
      if (g.reducible) return list(g.factors.map((x) => record(`(${fmtInts(x.factor)})`, [["group", x.group.name], ["order", String(x.group.order)]])), "reducible: one group per factor");
      return record(`Galois group ${g.name}`, [
        ["order", String(g.order)],
        ["transitive group", g.label],
        ["solvable by radicals", g.solvable ? "yes" : "no"],
        ...(g.discriminant ? [["discriminant", String(g.discriminant)], ["square?", g.discriminantIsSquare ? "yes" : "no"]] : []),
      ]);
    },
  },
  // continued fractions and quadratic fields
  cf: {
    args: [1, 2], doc: "cf(a, n): continued fraction terms",
    f: ([a, n]) => text(`[${A.continuedFraction(asNum(a), n ? Number(asInt(n)) : 20).map(String).join(", ").replace(", ", "; ")}${A.rationalOf(asNum(a)) ? "" : ", …"}]`),
  },
  period: {
    args: 1, doc: "periodic continued fraction of a quadratic irrational",
    f: ([a]) => {
      const { preperiod, period } = A.periodicContinuedFraction(asNum(a));
      return text(`[${preperiod.map(String).join(", ")}${preperiod.length ? "; " : ""}(${period.map(String).join(", ")}) repeating]`);
    },
  },
  convergents: {
    args: [1, 2], doc: "convergents(a, n)",
    f: ([a, n]) => text(A.convergents(A.continuedFraction(asNum(a), n ? Number(asInt(n)) : 8)).map(String).join(", ")),
  },
  pell: {
    args: 1, doc: "smallest solution of x² − D y² = 1",
    f: ([d]) => {
      const r = A.pell(asInt(d));
      return record(`x² − ${asInt(d)}·y² = 1`, [["x", String(r.x)], ["y", String(r.y)], ["x² − Dy² = −1", r.negative ? `x = ${r.negative.x}, y = ${r.negative.y}` : "no solution"]]);
    },
  },
  unit: {
    args: 1, doc: "fundamental unit of Q(√d)",
    f: ([d]) => {
      const r = A.fundamentalUnit(asInt(d));
      return record(`fundamental unit of Q(√${asInt(d)})`, [["ε", String(r.unit)], ["norm", String(r.norm)], ["≈", r.unit.toAlgebraic().toString(15)]]);
    },
  },
  // integers
  phi: { args: 1, doc: "Euler φ", f: ([n]) => A.fromRational(A.totient(asInt(n))) },
  isprime: { args: 1, doc: "primality (Miller–Rabin, certain below 3·10²⁴)", f: ([n]) => text(A.isProbablePrime(asInt(n)) ? "prime" : "not prime") },
  primroot: { args: 1, doc: "least primitive root mod a prime", f: ([p]) => A.fromRational(A.primitiveRoot(asInt(p))) },
  order: { args: 2, doc: "order(a, n): multiplicative order", f: ([a, n]) => A.fromRational(A.multiplicativeOrder(asInt(a), asInt(n))) },
  sqrtmod: {
    args: 2, doc: "sqrtmod(a, p)",
    f: ([a, p]) => {
      const r = A.sqrtMod(asInt(a), asInt(p));
      return r === null ? text(`${asInt(a)} is not a square mod ${asInt(p)}`) : A.fromRational(r);
    },
  },
  // guessing
  recognize: {
    args: [1, 2], doc: "recognize(\"1.41421356\"): guess an algebraic number from digits",
    raw: true,
    f: ([a, d], raw) => {
      const src = raw[0].k === "str" || raw[0].k === "num" ? raw[0].v : null;
      if (!src) throw new TypeError("recognize needs the digits as written, e.g. recognize(\"1.6180339887\")");
      const r = A.recognize(src, { maxDegree: d ? Number(asInt(d)) : 6 });
      if (!r) return text("no convincing algebraic number of small degree fits those digits");
      return num(r.value, `degree ${r.degree}, minimal polynomial ${fmtInts(r.poly)}; ${r.digitsUsed} digits given, about ${r.digitsNeeded} needed`);
    },
  },
};

export const FUNCTIONS = F;
export const CONSTANTS = { i: "the imaginary unit", x: "the polynomial variable" };

// ---------------------------------------------------------------------------
// Evaluation

function evaluate(node, env) {
  switch (node.k) {
    case "num": return A.fromRational(node.v);
    case "str": return { str: node.v };
    case "var": {
      if (node.name === "i") return A.I();
      if (node.name === "x") return P([[0n, 1n], ONE]);
      if (node.name in env) return env[node.name];
      throw new ReferenceError(`unknown name "${node.name}"`);
    }
    case "neg": {
      const a = evaluate(node.a, env);
      return isPoly(a) ? P(a.c.map(([n, d]) => [-n, d])) : A.neg(asNum(a));
    }
    case "+": case "-": case "*": case "/":
      return binary(node.k, evaluate(node.a, env), evaluate(node.b, env));
    case "^": return power(evaluate(node.a, env), evaluate(node.b, env));
    case "call": {
      const fn = F[node.name];
      if (!fn) throw new ReferenceError(`unknown function "${node.name}"`);
      const [lo, hi] = Array.isArray(fn.args) ? fn.args : [fn.args, fn.args];
      if (node.args.length < lo || node.args.length > hi) throw new TypeError(`${node.name} takes ${lo === hi ? lo : `${lo} to ${hi}`} argument${hi > 1 ? "s" : ""}`);
      const args = fn.raw ? node.args.map((a) => (a.k === "str" ? null : evaluate(a, env))) : node.args.map((a) => evaluate(a, env));
      return fn.f(args, node.args);
    }
    default: throw new Error("cannot evaluate that");
  }
}

/** Run one line: { result, name? } or throws. `env` keeps assignments. */
export function run(src, env = {}) {
  const ast = parse(src);
  if (ast.k === "let") {
    if (ast.name === "i" || ast.name === "x" || F[ast.name]) throw new SyntaxError(`"${ast.name}" is taken`);
    const v = evaluate(ast.e, env);
    env[ast.name] = v;
    return { name: ast.name, result: toResult(v) };
  }
  return { result: toResult(evaluate(ast, env)) };
}

export { parse };
