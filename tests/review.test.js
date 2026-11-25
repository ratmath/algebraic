// Fixes and additions from the 2026-11 review.
import { describe, test, expect } from "bun:test";
import {
  roots, isProbablePrime, primeField, tanPi, sinPi, cosPi, div, equals, minimalPolynomial, isolate, nthRoot, sqrt,
  fromRational as q, add, mul, neg, format, interpolate, NumberField, splittingField, fieldContaining, factorOver,
} from "../index.js";

describe("isProbablePrime (Baillie–PSW)", () => {
  test("agrees with a sieve below 20000", () => {
    const N = 20000;
    const composite = new Uint8Array(N + 1);
    for (let i = 2; i * i <= N; i += 1) if (!composite[i]) for (let j = i * i; j <= N; j += i) composite[j] = 1;
    for (let i = 0; i <= N; i += 1) expect(isProbablePrime(BigInt(i))).toBe(i >= 2 && !composite[i]);
  });
  test("rejects strong pseudoprimes to the first 12 prime bases", () => {
    // 318665857834031151167461 = 399165290221 × 798330580441 passes Miller–Rabin for bases 2 … 37
    expect(isProbablePrime(318665857834031151167461n)).toBe(false);
    expect(isProbablePrime(3825123056546413051n)).toBe(false);
    expect(isProbablePrime(3317044064679887385961981n)).toBe(false);
    expect(() => primeField(318665857834031151167461n)).toThrow();
  });
  test("rejects strong Lucas pseudoprimes and Carmichael numbers", () => {
    for (const n of [5459n, 5777n, 10877n, 561n, 1729n, 294409n]) expect(isProbablePrime(n)).toBe(false);
  });
  test("accepts large primes, and plain Numbers", () => {
    for (const p of [2n ** 61n - 1n, 2n ** 127n - 1n, 2n ** 521n - 1n, 18446744073709551557n]) expect(isProbablePrime(p)).toBe(true);
    expect(isProbablePrime(97)).toBe(true);
    expect(isProbablePrime(91)).toBe(false);
    expect(isProbablePrime(1)).toBe(false);
    expect(() => isProbablePrime(1.5)).toThrow(RangeError);
  });
});

describe("tanPi, directly from Im (1 + it)^q", () => {
  test("matches Math.tan for every p/q with q ≤ 16", () => {
    for (let d = 1; d <= 16; d += 1) {
      for (let p = -d; p <= d; p += 1) {
        if (d % 2 === 0 && (((p % d) + d) % d) === d / 2) continue;
        const v = tanPi([BigInt(p), BigInt(d)]).toNumber();
        const w = Math.tan((Math.PI * p) / d);
        expect(Math.abs(v - w)).toBeLessThan(1e-9 * Math.max(1, Math.abs(w)));
      }
    }
  });
  test("is the same number as sin/cos, on its minimal polynomial, and fast", () => {
    expect(equals(tanPi("1/7"), div(sinPi("1/7"), cosPi("1/7")))).toBe(true);
    expect(minimalPolynomial(tanPi("1/12"))).toEqual([1n, -4n, 1n]);
    expect(format(tanPi("1/3"))).toBe("√3");
    const t = Date.now();
    expect(tanPi("1/17").toNumber()).toBeCloseTo(Math.tan(Math.PI / 17), 12);
    expect(Date.now() - t).toBeLessThan(1000);
    expect(() => tanPi("1/2")).toThrow(/pole/);
  });
});

describe("refinement and arithmetic speed", () => {
  test("5000 digits of ∛2 and the root of x⁵ − x − 1 are right at the far end", () => {
    // digits 4991–5002 of ∛2 and 1991–2002 of √2 + √3, from an independent 80-digit-plus reference (mpmath)
    expect(nthRoot(q(2), 3).toString(5000).slice(4990, 5002)).toBe("692116759868");
    expect(add(sqrt(q(2)), sqrt(q(3))).toString(2000).slice(1990, 2002)).toBe("132971488269");
    const r = isolate([-1n, -1n, 0n, 0n, 0n, 1n])[0];
    expect(r.toString(40)).toBe("1.1673039782614186842560458998548421807205…");
  });
  test("sums and products of degree 40 and 49 stay real and exact", () => {
    const s = (n) => sqrt(q(n));
    const a = mul(nthRoot(q(2), 5), add(add(s(2), s(3)), s(5)));
    expect(minimalPolynomial(a).length - 1).toBe(40);
    expect(a.toNumber()).toBeCloseTo(2 ** 0.2 * (Math.SQRT2 + Math.sqrt(3) + Math.sqrt(5)), 10);
    const b = add(nthRoot(q(2), 7), nthRoot(q(3), 7));
    expect(minimalPolynomial(b).length - 1).toBe(49);
    expect(minimalPolynomial(add(s(2), s(3)))).toEqual([1n, 0n, -10n, 0n, 1n]);
    expect(minimalPolynomial(mul(s(2), s(3)))).toEqual([-6n, 0n, 1n]);
  });
});

describe("radical output for pure roots", () => {
  const c = (n, k) => nthRoot(q(n), k);
  test.each([
    [() => c(2, 3), "∛2"],
    [() => neg(c(2, 3)), "−∛2"],
    [() => c("1/2", 3), "∛4/2"],
    [() => c(16, 3), "2∛2"],
    [() => c(-5, 3), "−∛5"],
    [() => c(7, 5), "⁵√7"],
    [() => c(5, 4), "∜5"],
    [() => neg(c(5, 4)), "−∜5"],
    [() => add(q(1), c(2, 3)), "1 + ∛2"],
    [() => add(q("-1/2"), mul(q(3), c(2, 3))), "−1/2 + 3∛2"],
    [() => cosPi("1/8"), "√(2 + √2)/2"],
  ])("format %#", (make, text) => {
    expect(format(make())).toBe(text);
  });
});

describe("interpolate", () => {
  test("names a repeated x", () => {
    expect(() => interpolate([[1, 2], [1, 3]])).toThrow(/distinct x values; 1 appears twice/);
  });
});

describe("polynomials over number fields", () => {
  const show = (F) => F.factors.map((f) => `(${f.factor.map(String).join(", ")})^${f.multiplicity}`).join(" ");
  const K = NumberField.quadratic(2);
  test("factoring by Trager's method", () => {
    expect(show(K.factor([-4, 0, 0, 0, 1]))).toBe("(−√2, 1)^1 (√2, 1)^1 (2, 0, 1)^1");
    expect(show(K.factor([1, 0, 0, 0, 1]))).toBe("(1, √2, 1)^1 (1, −√2, 1)^1");
    expect(show(NumberField.gaussian().factor([1, 0, 0, 0, 1]))).toBe("(i, 0, 1)^1 (−i, 0, 1)^1");
    const a = K.gen;
    expect(show(factorOver(K, [a.mul(a), a.mul(-2), K.of(3), a.mul(-2), K.one]))).toBe("(−√2, 1)^2 (1, 0, 1)^1");
    expect(show(new NumberField([-2, 0, 0, 1]).factor([-2, 0, 0, 1]))).toBe("(−α, 1)^1 (α², α, 1)^1");
  });
  test("the factors multiply back", () => {
    const F = K.factor([3, 1, -7, 0, 2, 1]);
    let prod = [F.unit];
    for (const { factor: f, multiplicity } of F.factors) {
      for (let k = 0; k < multiplicity; k += 1) {
        const o = Array(prod.length + f.length - 1).fill(K.zero);
        prod.forEach((x, i) => f.forEach((y, j) => { o[i + j] = o[i + j].add(x.mul(y)); }));
        prod = o;
      }
    }
    expect(prod.map(String)).toEqual(["3", "1", "−7", "0", "2", "1"]);
  });
  test("roots in the field", () => {
    expect(K.roots([-8, 0, 1]).map((r) => String(r.value))).toEqual(["2√2", "−2√2"]);
    expect(K.roots([-3, 0, 1])).toEqual([]);
  });
  test("extend keeps the old generator's value", () => {
    const E = K.extend([-3, 0, 1]);
    expect(E.field.n).toBe(4);
    expect(String(E.root.mul(E.root))).toBe("3");
    expect(equals(E.embed(K.gen).toAlgebraic(), sqrt(q(2)))).toBe(true);
    expect(() => K.extend([-2, 0, 1])).toThrow(/irreducible/);
  });
  test("one field for several numbers", () => {
    const xs = [sqrt(q(2)), sqrt(q(3)), q("1/2"), sqrt(q(8))];
    const C = fieldContaining(xs);
    expect(C.field.n).toBe(4);
    C.elements.forEach((e, i) => expect(equals(e.toAlgebraic(), xs[i])).toBe(true));
    expect(format(C.elements[0].mul(C.elements[1]).toAlgebraic())).toBe("√6");
    expect(String(C.elements[3].sub(C.elements[0].mul(2)))).toBe("0");
  });
  test("splitting fields", () => {
    const S3 = splittingField([-2, 0, 0, 1]);
    expect(S3.degree).toBe(6);
    expect(S3.roots.length).toBe(3);
    for (const r of S3.roots) expect(String(r.pow(3))).toBe("2");
    expect(splittingField([-2, 0, 0, 0, 1]).degree).toBe(8);
    expect(NumberField.splitting([1, 0, -10, 0, 1]).degree).toBe(4);
  });
});

describe("format of negated denested radicals", () => {
  test("the four roots of x⁴ − 10x² + 1", () => {
    expect(roots([1, 0, -10, 0, 1]).map((r) => format(r.value))).toEqual(["−√3 − √2", "−√3 + √2", "√3 − √2", "√3 + √2"]);
    expect(roots([2, 0, -4, 0, 1]).map((r) => format(r.value))).toEqual(["−√(2 + √2)", "−√(2 − √2)", "√(2 − √2)", "√(2 + √2)"]);
  });
});

