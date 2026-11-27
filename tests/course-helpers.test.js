// Helpers added for the Course Test roadmap, batch A6 (ratmath/shelf course-test/ROADMAP.md).
import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import { sqrt, nthRoot, add, fromRational, AlgebraicReal } from "../index.js";

// independent check: floor(√(n·10^(2k))) by integer square root
const isqrt = (n) => {
  if (n < 2n) return n;
  let x = n;
  let y = (x + 1n) / 2n;
  while (y < x) [x, y] = [y, (y + n / y) / 2n];
  return x;
};

describe("rounding algebraic numbers to k places", () => {
  test("√2 to 40 places matches the integer square root", () => {
    const k = 40;
    const s = isqrt(2n * 10n ** BigInt(2 * k + 4)); // two extra digits
    const lastTwo = s % 100n;
    const truncated = s / 100n;
    const rounded = lastTwo >= 50n ? truncated + 1n : truncated; // √2 is irrational, so no tie
    const digits = rounded.toString();
    expect(sqrt(2n).toFixed(k)).toBe(`${digits[0]}.${digits.slice(1)}`);
    expect(sqrt(2n).toFixed(k, "floor")).toBe(`1.${(isqrt(2n * 10n ** BigInt(2 * k)) % 10n ** BigInt(k)).toString().padStart(k, "0")}`);
  });
  test("modes on negative and positive irrationals", () => {
    const neg = add(fromRational(0), sqrt(3n)); // √3 = 1.7320508…
    expect(neg.toFixed(3)).toBe("1.732");
    expect(neg.toFixed(3, "ceil")).toBe("1.733");
    expect(neg.toFixed(0)).toBe("2");
    const minus = add(fromRational(-1), nthRoot(fromRational(2), 3)); // ∛2 − 1 = 0.2599…
    expect(minus.toFixed(2)).toBe("0.26");
    expect(minus.toFixed(2, "toward-zero")).toBe("0.25");
    const below = add(fromRational(1), sqrt(2n)).roundTo(4, "floor");
    expect(below.equals(new Rational(24142n, 10000n))).toBe(true);
  });
  test("rational values: ties follow the mode exactly", () => {
    const half = AlgebraicReal.rational([5n, 2n]);
    expect(half.toFixed(0)).toBe("2"); // half-even
    expect(half.toFixed(0, "half-up")).toBe("3");
    const neg = AlgebraicReal.rational([-1n, 8n]);
    expect(neg.toFixed(2)).toBe("-0.12");
    expect(neg.toFixed(2, "half-up")).toBe("-0.13");
    expect(neg.toFixed(2, "floor")).toBe("-0.13");
    expect(AlgebraicReal.rational([1n, 3n]).toFixed(5)).toBe("0.33333");
    // a number that looks irrational but is the rational 3/2 after squaring back
    expect(sqrt(fromRational("9/4")).toFixed(0, "half-up")).toBe("2");
  });
  test("bad arguments are refused", () => {
    expect(() => sqrt(2n).toFixed(-1)).toThrow(RangeError);
    expect(() => sqrt(2n).toFixed(2, "up")).toThrow(/rounding mode/);
  });
});
