import { describe, expect, test } from "bun:test";
import {
  fromRational, roots, minimalPolynomial, simplify, add, sub, mul, div, inv, pow, sign, compare, equals, isZero,
  I, conjugate, re, im, nthRoot, sqrt, isReal, factor, isIrreducible, factorModP,
  NumberField, primeField, galoisField, AlgebraicReal,
} from "../index.js";

const q = (x) => fromRational(x);
const big = (a) => a.map(BigInt);
const polyMul = (a, b) => {
  const o = Array(a.length + b.length - 1).fill(0n);
  a.forEach((x, i) => b.forEach((y, j) => { o[i + j] += x * y; }));
  return o;
};

describe("factoring over the integers", () => {
  const show = (r) => r.factors.map((f) => `${f.factor.join(",")}^${f.multiplicity}`).sort();
  test("small classics", () => {
    expect(show(factor(big([-1, 0, 1])))).toEqual(["-1,1^1", "1,1^1"].sort());
    expect(show(factor(big([1, 0, 0, 0, 1])))).toEqual(["1,0,0,0,1^1"]);
    expect(show(factor(big([-1, 0, 0, 0, 0, 0, 1])))).toEqual(["-1,1^1", "1,1^1", "1,-1,1^1", "1,1,1^1"].sort());
    expect(show(factor(big([1, 5, 6])))).toEqual(["1,2^1", "1,3^1"].sort());
    expect(factor(big([0, -4, 0, 2])).content).toBe(2n);
  });
  test("multiplicities and a Swinnerton-Dyer polynomial", () => {
    const p = polyMul(polyMul(polyMul(big([-1, 1]), big([-1, 1])), big([-1, 1])), big([2, 1]));
    expect(show(factor(p))).toEqual(["-1,1^3", "2,1^1"].sort());
    expect(isIrreducible(big([576, 0, -960, 0, 352, 0, -40, 0, 1]))).toBe(true);
  });
  test("random products come apart into their factors", () => {
    let s = 11n;
    const r = () => (s = (s * 48271n) % 2147483647n) % 19n - 9n;
    for (let trial = 0; trial < 20; trial += 1) {
      const a = [...Array.from({ length: 4 }, r), 1n + (s % 3n)];
      const b = [...Array.from({ length: 3 }, r), 2n];
      if (a[0] === 0n || b[0] === 0n) continue;
      const p = polyMul(a, b);
      const { content, factors } = factor(p);
      const back = factors.reduce((acc, f) => {
        let x = acc;
        for (let k = 0; k < f.multiplicity; k += 1) x = polyMul(x, f.factor);
        return x;
      }, [content]);
      expect(back).toEqual(p);
      for (const f of factors) expect(factor(f.factor).factors.length).toBe(1);
    }
  });
  test("mod p", () => {
    expect(factorModP(big([1, 0, 1]), 5n).length).toBe(2); // x² + 1 = (x − 2)(x − 3) mod 5
    expect(factorModP(big([1, 0, 1]), 7n).length).toBe(1);
    expect(factorModP(big([0, 1, ...Array(14).fill(0), 1]), 2n).map((f) => f.length - 1)).toEqual([1, 1, 2, 4, 4, 4]); // x¹⁶ + x over F₂
    expect(() => factorModP(big([1, 0, 0, 0, 1]), 2n)).toThrow(); // (x + 1)⁴ mod 2
  });
});

describe("real and complex algebraic arithmetic", () => {
  const s2 = sqrt(q(2));
  const s3 = sqrt(q(3));
  test("results carry minimal polynomials", () => {
    expect(add(s2, s3).poly).toEqual(big([1, 0, -10, 0, 1]));
    expect(pow(add(s2, s3), 2).poly).toEqual(big([1, -10, 1]));
    expect(mul(s2, s3).poly).toEqual(big([-6, 0, 1]));
    expect(mul(s2, s2).toRational().toString()).toBe("2");
    expect(isZero(sub(s2, s2))).toBe(true);
    expect(minimalPolynomial(nthRoot(q(2), 3))).toEqual(big([-2, 0, 0, 1]));
  });
  test("comparison and equality are exact", () => {
    expect(compare(s2, q("7/5"))).toBe(1);
    expect(compare(add(s2, s3), q("3.1462643"))).toBe(1);
    expect(equals(mul(s2, s3), sqrt(q(6)))).toBe(true);
    expect(sign(sub(s3, s2))).toBe(1);
    expect(div(q(1), s2).toNumber()).toBeCloseTo(Math.SQRT1_2, 14);
    expect(inv(add(s2, s3)).toNumber()).toBeCloseTo(1 / (Math.SQRT2 + Math.sqrt(3)), 14);
  });
  test("complex numbers", () => {
    const i = I();
    expect(mul(i, i).toRational().toString()).toBe("-1");
    const z = add(q(1), i);
    expect(pow(z, 2).poly).toEqual(big([4, 0, 1]));
    const w = inv(z);
    expect(re(w).toRational().toString()).toBe("1/2");
    expect(im(w).toRational().toString()).toBe("-1/2");
    expect(equals(conjugate(z), sub(q(1), i))).toBe(true);
    expect(sqrt(q(-4)).toNumbers().im).toBeCloseTo(2, 12);
    expect(mul(sqrt(q(-2)), sqrt(q(-2))).toRational().toString()).toBe("-2");
  });
  test("roots of a polynomial, with complex ones on request", () => {
    expect(roots([-2, 0, 0, 1]).length).toBe(1);
    const all = roots([-2, 0, 0, 1], { complex: true });
    expect(all.length).toBe(3);
    expect(isReal(all[0].value)).toBe(true);
    const r = roots([2n, -2n, -1n, 1n]); // (x² − 2)(x − 1)
    expect(r.map((x) => x.value.poly.length - 1)).toEqual([2, 1, 2]);
    expect(simplify(new AlgebraicReal(big([2, -2, -1, 1]), [5n, 4n], [3n, 2n])).poly).toEqual(big([-2, 0, 1]));
  });
});

describe("number fields", () => {
  test("Q(√2)", () => {
    const K = NumberField.quadratic(2);
    const a = K.of([1, 1]);
    expect(String(a)).toBe("1 + √2");
    expect(String(a.inv())).toBe("−1 + √2");
    expect(a.norm().toString()).toBe("-1");
    expect(a.trace().toString()).toBe("2");
    expect(a.toAlgebraic().toNumber()).toBeCloseTo(1 + Math.SQRT2, 14);
  });
  test("Gaussian rationals and a cubic field", () => {
    const G = NumberField.gaussian();
    const z = G.of([1, 2]);
    expect(String(z.mul(z))).toBe("−3 + 4i");
    expect(z.norm().toString()).toBe("5");
    expect(z.mul(z.inv()).equals(G.one)).toBe(true);
    const C = new NumberField([-2, 0, 0, 1], { name: "∛2" });
    expect(C.gen.pow(3).toRational().toString()).toBe("2");
    expect(() => new NumberField([-4, 0, 1])).toThrow();
  });
  test("works as a field object", () => {
    const K = NumberField.quadratic(5);
    expect(K.sign(K.sub(K.gen, K.gen))).toBe(0);
    expect(K.mul(K.gen, K.gen).toRational().toString()).toBe("5");
  });
});

describe("finite fields", () => {
  test("F_p", () => {
    const F = primeField(7);
    expect(F.mul(3n, F.inv(3n))).toBe(1n);
    expect(F.from({ numerator: 1n, denominator: 3n })).toBe(5n);
    expect(() => primeField(9)).toThrow();
  });
  test("GF(p^k)", () => {
    const F8 = galoisField(2, 3);
    expect(F8.mul(F8.gen, F8.inv(F8.gen))).toEqual(F8.one);
    const F9 = galoisField(3, 2);
    expect(F9.pow(F9.gen, 8)).toEqual(F9.one);
    let x = F9.one;
    const seen = new Set();
    for (let k = 0; k < 8; k += 1) {
      seen.add(x.join(","));
      x = F9.mul(x, F9.add(F9.gen, F9.one));
    }
    expect(seen.size).toBeGreaterThan(1);
  });
});
