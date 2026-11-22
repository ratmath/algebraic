import { describe, expect, test } from "bun:test";
import { recognize, lll, format } from "../index.js";

describe("recognising decimals", () => {
  test("LLL finds the short vector", () => {
    const red = lll([[1n, 0n, 0n, 10000n], [0n, 1n, 0n, 14142n], [0n, 0n, 1n, 20000n]]);
    expect(red[0].slice(0, 3).map(String)).toEqual(["-2", "0", "1"].map((x, i, a) => (red[0][2] < 0n ? String(-BigInt(x)) : x)));
  });
  test("small algebraic numbers come back exactly", () => {
    expect(format(recognize("1.41421356237").value)).toBe("√2");
    expect(format(recognize("1.6180339887").value)).toBe("(1 + √5)/2");
    expect(recognize("1.2599210498948732").poly).toEqual([-2n, 0n, 0n, 1n]);
    expect(format(recognize("3.14626436994197234232").value)).toBe("√3 + √2");
    expect(format(recognize("0.9238795325112867").value)).toBe("√(2 + √2)/2");
  });
  test("coincidences are refused", () => {
    expect(recognize("3.14159265358979")).toBe(null);
  });
});
