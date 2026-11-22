import { describe, expect, test } from "bun:test";
import { cosPi, sinPi, tanPi, rootOfUnity, rootsOfUnity, cyclotomic, format, pow, equals, fromRational, add } from "../index.js";

describe("exact trigonometric values", () => {
  test("cos and sin at rational multiples of π", () => {
    expect(format(cosPi("1/3"))).toBe("1/2");
    expect(format(cosPi("1/4"))).toBe("√2/2");
    expect(format(cosPi("1/5"))).toBe("(1 + √5)/4");
    expect(format(cosPi("7/6"))).toBe("−√3/2");
    expect(format(cosPi("1/8"))).toBe("√(2 + √2)/2");
    expect(format(cosPi("1/12"))).toBe("(√6 + √2)/4");
    expect(format(sinPi("1/12"))).toBe("(√6 − √2)/4");
    expect(format(sinPi("1/5"))).toBe("√(10 − 2√5)/4");
    expect(format(tanPi("1/12"))).toBe("2 − √3");
    expect(cosPi("2/7").poly).toEqual([-1n, -4n, 4n, 8n]); // 8x³ + 4x² − 4x − 1
    expect(cosPi("1/100").toString(20)).toBe("0.99950656036573155700…");
    expect(() => tanPi("1/2")).toThrow();
  });
  test("Pythagoras holds exactly", () => {
    for (const r of ["1/7", "2/9", "3/10"]) {
      expect(equals(add(pow(cosPi(r), 2), pow(sinPi(r), 2)), fromRational(1))).toBe(true);
    }
  });
});

describe("roots of unity", () => {
  test("cyclotomic polynomials", () => {
    expect(cyclotomic(12)).toEqual([1n, 0n, -1n, 0n, 1n]);
    expect(cyclotomic(105).some((c) => c === -2n)).toBe(true);
  });
  test("e^(2πik/n)", () => {
    expect(rootsOfUnity(8).map((z) => format(z))).toEqual(["1", "√2/2 + i√2/2", "i", "−√2/2 + i√2/2", "−1", "−√2/2 − i√2/2", "−i", "√2/2 − i√2/2"]);
    expect(format(rootOfUnity(5))).toBe("(−1 + √5)/4 + i√(10 + 2√5)/4");
    const z = rootOfUnity(17, 3);
    expect(pow(z, 17).toRational().toString()).toBe("1");
  });
});
