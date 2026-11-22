import { describe, expect, test } from "bun:test";
import { factorInteger, crt, modPow, modInverse, jacobi, sqrtMod, primitiveRoot, multiplicativeOrder, discreteLogMod, totient, primeField, galoisField } from "../index.js";

describe("number theory", () => {
  test("factoring and φ", () => {
    expect(factorInteger(2n ** 67n - 1n)).toEqual([{ prime: 193707721n, exponent: 1 }, { prime: 761838257287n, exponent: 1 }]);
    expect(factorInteger(360)).toEqual([{ prime: 2n, exponent: 3 }, { prime: 3n, exponent: 2 }, { prime: 5n, exponent: 1 }]);
    expect(totient(36)).toBe(12n);
  });
  test("modular arithmetic", () => {
    expect(crt([2, 3, 2], [3, 5, 7])).toEqual({ x: 23n, m: 105n });
    expect(crt([1, 2], [4, 6])).toBe(null);
    expect(modPow(3, 200, 1000)).toBe(1n);
    expect(modInverse(3, 7)).toBe(5n);
    expect(jacobi(2, 7)).toBe(1);
    expect(sqrtMod(10, 13)).toBe(6n);
    expect(sqrtMod(5, 13)).toBe(null);
    expect(primitiveRoot(7)).toBe(3n);
    expect(multiplicativeOrder(2, 7)).toBe(3n);
    expect(discreteLogMod(3, 13, 17)).toBe(4n);
  });
  test("orders, generators and logs in finite fields", () => {
    const F = primeField(101);
    expect(F.order(2n)).toBe(100n);
    expect(F.pow(F.primitive(), F.log(5n))).toBe(5n);
    const G = galoisField(2, 8);
    const g = G.primitive();
    expect(G.order(g)).toBe(255n);
    expect(G.log(G.mul(g, g))).toBe(2n);
    expect(galoisField(3, 2).elements().length).toBe(9);
  });
});
