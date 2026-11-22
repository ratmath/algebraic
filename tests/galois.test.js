import { describe, expect, test } from "bun:test";
import { galoisGroup } from "../index.js";

const name = (p) => galoisGroup(p).name;

describe("Galois groups", () => {
  test("degrees 2 to 4", () => {
    expect(name([-2, 0, 1])).toBe("C2");
    expect(name([-2, 0, 0, 1])).toBe("S3");
    expect(name([1, -3, 0, 1])).toBe("C3");
    expect(name([-2, 0, 0, 0, 1])).toBe("D4");
    expect(name([1, 0, 0, 0, 1])).toBe("V4");
    expect(name([2, 0, -4, 0, 1])).toBe("C4");
    expect(name([12, 8, 0, 0, 1])).toBe("A4");
    expect(name([1, 1, 0, 0, 1])).toBe("S4");
    expect(name([1, 0, 0, 0, "1/2"])).toBe("D4");
  });
  test("quintics through certified resolvents", () => {
    expect(name([-1, -1, 0, 0, 0, 1])).toBe("S5");
    expect(name([-2, 0, 0, 0, 0, 1])).toBe("F20");
    expect(name([16, 20, 0, 0, 0, 1])).toBe("A5");
    expect(name([12, -5, 0, 0, 0, 1])).toBe("D5");
    expect(name([1, 3, -3, -4, 1, 1])).toBe("C5");
    const g = galoisGroup([-1, -1, 0, 0, 0, 1]);
    expect(g.order).toBe(120);
    expect(g.solvable).toBe(false);
  });
  test("reducible input reports each factor", () => {
    const r = galoisGroup([-1, 0, 0, 0, 0, 1]);
    expect(r.reducible).toBe(true);
    expect(r.factors.map((f) => f.group.name)).toEqual(["1", "C4"]);
  });
});
