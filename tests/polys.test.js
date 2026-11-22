import { describe, expect, test } from "bun:test";
import { polys as P } from "../index.js";

const s = (p) => p.map(String);
const show = (pf) => pf.terms.map((t) => `${s(t.numerator).join(",")}/(${t.factor.join(",")})^${t.power}`);

describe("polynomial toolkit", () => {
  test("resultants and discriminants", () => {
    expect(String(P.resultant([-2, 0, 1], [-3, 0, 1]))).toBe("1");
    expect(String(P.resultant([-1, 1], [-1, 0, 1]))).toBe("0");
    expect(String(P.discriminant([1, 1, 1]))).toBe("-3");
    expect(String(P.discriminant([-1, -1, 0, 1]))).toBe("-23");
    expect(String(P.discriminant(["1/2", 0, 2]))).toBe("-4");
  });
  test("root counting, interpolation, gcd", () => {
    expect(P.countRoots([-2, 0, 1], 0, 2)).toBe(1);
    expect(P.countRoots([0, -1, 0, 1])).toBe(3);
    expect(P.countRoots([1, -2, 1])).toBe(1);
    expect(s(P.interpolate([[0, 1], [1, 3], [2, 7]]))).toEqual(["1", "1", "1"]);
    const g = P.gcdExt([-1, 0, 1], [1, 1]);
    expect(s(g.gcd)).toEqual(["1", "1"]);
    expect(P.format(P.compose([0, 0, 1], [1, 1]))).toBe("x² + 2x + 1");
    expect(s(P.fromRoots([1, 2]))).toEqual(["2", "-3", "1"]);
  });
  test("partial fractions", () => {
    expect(show(P.partialFractions([1], [0, -1, 0, 1]))).toEqual(["1/2/(-1,1)^1", "-1/(0,1)^1", "1/2/(1,1)^1"]);
    const pf = P.partialFractions([1, 0, 0, 0, 1], [0, 0, 1, 0, 1]);
    expect(s(pf.polynomial)).toEqual(["1"]);
    expect(show(pf)).toEqual(["1/(0,1)^2", "-2/(1,0,1)^1"]);
    expect(show(P.partialFractions([3, 2], [-1, 3, -3, 1]))).toEqual(["2/(-1,1)^2", "5/(-1,1)^3"]);
  });
});
