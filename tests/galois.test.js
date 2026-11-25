import { describe, expect, test } from "bun:test";
import { galoisGroup, galoisAction, PermGroup, cycleString, fromCycles, polys } from "../index.js";

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
  test("sextics", () => {
    const label = (p) => galoisGroup(p).label;
    expect(label([1, 1, 1, 1, 1, 1, 1])).toBe("6T1"); // the 7th cyclotomic polynomial
    expect(label([108, 0, 0, 0, 0, 0, 1])).toBe("6T2"); // the splitting field of x³ − 2
    expect(label([-2, 0, 0, 0, 0, 0, 1])).toBe("6T3");
    expect(label([-1, 0, -3, 0, 0, 0, 1])).toBe("6T4");
    expect(label([3, 0, 0, -3, 0, 0, 1])).toBe("6T5");
    expect(label([1, 0, 0, 0, -3, 0, 1])).toBe("6T6");
    expect(label([-1, 0, -3, 0, -3, 0, 1])).toBe("6T7");
    expect(label([2, 0, 2, 0, -1, 0, 1])).toBe("6T8");
    expect(label([-3, 0, 0, -3, 0, 0, 1])).toBe("6T9");
    expect(label([-2, 0, -3, 0, -3, 0, 1])).toBe("6T11");
    expect(label([25, 170, 175, 140, 55, 10, 1])).toBe("6T12");
    expect(label([2, 8, 4, -4, -4, 0, 1])).toBe("6T13");
    expect(label([7, 6, 5, 4, 3, 2, 1])).toBe("6T14");
    const s6 = galoisGroup([-1, -1, 0, 0, 0, 0, 1]);
    expect([s6.name, s6.order, s6.solvable]).toEqual(["S6", 720, false]);
  }, 30000);
  test("septics", () => {
    const name = (p) => galoisGroup(p).name;
    expect(name([1, -9, 14, 28, -7, -12, 1, 1])).toBe("C7");
    expect(name([-1, 7, 7, 7, 0, 0, 0, 1])).toBe("D7");
    expect(name([-2, -6, 6, 16, -2, -8, 0, 1])).toBe("F21");
    expect(name([-2, 0, 0, 0, 0, 0, 0, 1])).toBe("F42");
    expect(name([3, -7, 0, 0, 0, 0, 0, 1])).toBe("PSL(3,2)");
    expect(name([-1, -1, 0, 0, 0, 0, 0, 1])).toBe("S7");
    expect(() => galoisGroup([-2, 0, 0, 0, 0, 0, 0, 0, 1])).toThrow(RangeError);
  }, 30000);
});

describe("the group acting on the roots", () => {
  test("agrees with the classification by resolvents in degree ≤ 5", () => {
    for (const p of [[-2, 0, 1], [-2, 0, 0, 1], [1, -3, 0, 1], [-2, 0, 0, 0, 1], [1, 0, 0, 0, 1], [2, 0, -4, 0, 1], [12, 8, 0, 0, 1], [1, 1, 0, 0, 1],
      [-1, -1, 0, 0, 0, 1], [-2, 0, 0, 0, 0, 1], [16, 20, 0, 0, 0, 1], [12, -5, 0, 0, 0, 1], [1, 3, -3, -4, 1, 1]]) {
      const a = galoisAction(p);
      const g = galoisGroup(p);
      expect([a.label, a.order, a.solvable]).toEqual([g.label, g.order, g.solvable]);
      expect(a.group.isTransitive).toBe(true);
    }
  });
  test("x⁴ − 2: D4 with complex conjugation as an element", () => {
    const a = galoisAction([-2, 0, 0, 0, 1]);
    expect(a.name).toBe("D4");
    expect(a.roots.length).toBe(4);
    expect(a.permutations().length).toBe(8);
    expect(a.group.has(a.conjugation)).toBe(true);
    expect(cycleString(a.conjugation)).toBe("(3 4)"); // two real roots, then the conjugate pair
    expect(a.derivedSeries()).toEqual([8, 2, 1]);
  });
  test("the Galois correspondence for x⁴ − 2", () => {
    const a = galoisAction([-2, 0, 0, 0, 1]);
    const rows = a.correspondence();
    expect(rows.map((r) => r.order)).toEqual([1, 2, 2, 2, 4, 4, 4, 8]);
    for (const r of rows) {
      expect(r.field.degree).toBe(r.index);
      expect(r.field.polynomial.length - 1).toBe(r.index);
    }
    // the three quadratic subfields of Q(⁴√2, i) are Q(i), Q(√2), Q(√−2)
    const quads = rows.filter((r) => r.index === 2).map((r) => polys.discriminant(r.field.polynomial));
    const squarefree = quads.map((d) => { let n = d.numerator; for (let k = 2n; k * k <= (n < 0n ? -n : n); k += 1n) while (n % (k * k) === 0n) n /= k * k; return Number(n); });
    expect(squarefree.sort((x, y) => x - y)).toEqual([-2, -1, 2]);
    expect(rows.filter((r) => r.normal).length).toBe(6);
    // a root generates the field fixed by its stabiliser
    const stab = a.group.subgroup(a.group.elements.filter((g) => g[0] === 0));
    expect(a.fixedField(stab).polynomial).toEqual([-2n, 0n, 0n, 0n, 1n]);
    expect(a.fixedField(a.group).degree).toBe(1);
  });
  test("x³ − 2: the fixed field of the rotations is Q(√−3)", () => {
    const a = galoisAction([-2, 0, 0, 1]);
    const f = a.fixedField("(1 2 3)");
    expect(f.degree).toBe(2);
    expect(f.normal).toBe(true);
    const d = polys.discriminant(f.polynomial).numerator;
    expect(d < 0n && d % 3n === 0n).toBe(true);
    // the splitting field has degree 6; the generator needs moved roots here
    const top = a.fixedField([]);
    expect(top.degree).toBe(6);
    expect(typeof top.substitution).toBe("string");
  });
  test("a non-monic polynomial", () => {
    const a = galoisAction([-3, 0, 0, 0, 2]);
    expect(a.name).toBe("D4");
    expect(a.fixedField(a.group.subgroup(a.group.elements.filter((g) => g[0] === 0))).substitution).toBe("y = 2x");
  });
  test("factor patterns modulo primes are cycle types (Dedekind, Chebotarev)", () => {
    const a = galoisAction([3, -7, 0, 0, 0, 0, 0, 1]); // PSL(3,2)
    const { primes, patterns } = a.frobenius(300);
    expect(primes).toBe(300);
    expect(patterns.reduce((s, x) => s + x.primes, 0)).toBe(300);
    expect(patterns.map((x) => x.type.join("")).sort()).toEqual(["1111111", "22111", "331", "421", "7"]);
    for (const x of patterns) expect(Math.abs(x.observed - x.expected)).toBeLessThan(0.08);
    // 7-cycles are 2/7 of the group
    expect(patterns.find((x) => x.type[0] === 7).elements).toBe(48);
  });
  test("subgroups of PSL(3,2)", () => {
    const a = galoisAction([3, -7, 0, 0, 0, 0, 0, 1]);
    expect(a.group.subgroups().length).toBe(179);
    const rows = a.correspondence({ maxDegree: 8 });
    expect(rows.length).toBe(15);
    // two classes of index 7: the points and the lines of the Fano plane
    const sevens = rows.filter((r) => r.index === 7);
    expect(sevens.length).toBe(2);
    expect(sevens.map((r) => galoisGroup(r.field.polynomial).name)).toEqual(["PSL(3,2)", "PSL(3,2)"]);
    expect(rows[0].field).toBe(null);
  });
  test("reducible input is refused", () => {
    expect(() => galoisAction([-1, 0, 0, 0, 1])).toThrow(RangeError);
  });
});

describe("PermGroup", () => {
  test("orders and structure", () => {
    expect(PermGroup.symmetric(5).order).toBe(120);
    expect(PermGroup.alternating(6).order).toBe(360);
    const D4 = PermGroup.from(["(1 2 3 4)", "(1 3)"], 4);
    expect(D4.order).toBe(8);
    expect(D4.isAbelian).toBe(false);
    expect(D4.subgroups().length).toBe(10);
    expect(D4.subgroupClasses().length).toBe(8);
    expect(D4.derivedSubgroup().order).toBe(2);
    expect(PermGroup.symmetric(4).derivedSeries().map((H) => H.order)).toEqual([24, 12, 4, 1]);
    expect(PermGroup.symmetric(5).isSolvable).toBe(false);
    expect(PermGroup.symmetric(4).subgroups().length).toBe(30);
    expect(PermGroup.symmetric(5).subgroupClasses().length).toBe(19);
  });
  test("cycle notation", () => {
    expect(cycleString([1, 2, 0, 4, 3])).toBe("(1 2 3)(4 5)");
    expect(fromCycles("(1 2 3)(4 5)", 5)).toEqual([1, 2, 0, 4, 3]);
    expect(cycleString([0, 1, 2])).toBe("()");
    expect(PermGroup.symmetric(4).cycleTypes().find((t) => t.type.join("") === "22").count).toBe(3);
    expect(() => fromCycles("(1 9)", 4)).toThrow(RangeError);
  });
});
