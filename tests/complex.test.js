import { describe, expect, test } from "bun:test";
import {
  fromRational, I, add, sub, mul, pow, equals, sqrt, nthRoot, complex, norm, abs, format, re, im, roots, isReal, conjugate,
} from "../index.js";

const q = (x) => fromRational(x);

describe("principal roots of complex numbers", () => {
  test("square roots", () => {
    expect(String(sqrt(complex(-3, 4)))).toBe("1 + 2i");
    expect(String(sqrt(complex(-3, -4)))).toBe("1 − 2i");
    expect(equals(pow(sqrt(I()), 2), I())).toBe(true);
    expect(format(sqrt(I()))).toBe("√2/2 + i√2/2");
    expect(String(sqrt(q(-4)))).toBe("2i");
  });
  test("higher roots pick the largest real part", () => {
    expect(format(nthRoot(I(), 3))).toBe("√3/2 + i/2");
    expect(format(nthRoot(complex(0, -8), 3))).toBe("√3 − i");
    expect(format(nthRoot(q(-16), 4))).toBe("√2 + i√2");
    expect(format(nthRoot(q(-1), 6))).toBe("√3/2 + i/2");
    expect(nthRoot(q(-8), 3).toRational().toString()).toBe("-2"); // odd roots of negative reals stay real
    const z = complex(1, sqrt(q(2)));
    expect(equals(pow(nthRoot(z, 3), 3), z)).toBe(true);
  });
});

describe("modulus and parts", () => {
  test("abs and norm", () => {
    expect(abs(complex(3, 4)).toRational().toString()).toBe("5");
    expect(format(abs(complex(1, 1)))).toBe("√2");
    expect(format(abs(add(sqrt(q(2)), I())))).toBe("√3");
    expect(norm(complex(1, 2)).toRational().toString()).toBe("5");
    expect(abs(q(-3)).toRational().toString()).toBe("3");
  });
  test("real and imaginary parts of quadratic numbers are read off exactly", () => {
    const w = roots([1, 1, 1], { complex: true })[0].value;
    expect(format(re(w))).toBe("−1/2");
    expect(format(im(w))).toBe("√3/2");
    expect(format(im(conjugate(w)))).toBe("−√3/2");
  });
});

describe("display", () => {
  test("quadratic numbers print exactly", () => {
    expect(String(pow(complex(1, 1), 3))).toBe("−2 + 2i");
    expect(String(complex("1/2", "1/2"))).toBe("1/2 + i/2");
    expect(roots([1, 1, 1], { complex: true }).map((r) => String(r.value))).toEqual(["−1/2 + i√3/2", "−1/2 − i√3/2"]);
    expect(roots([-1, -1, 1]).map((r) => format(r.value))).toEqual(["(1 − √5)/2", "(1 + √5)/2"]);
    expect(format(sub(q(0), mul(sqrt(q(2)), q("1/2"))))).toBe("−√2/2");
  });
  test("other numbers print certified digits", () => {
    const [, z] = roots([-2, 0, 0, 1], { complex: true }).map((r) => r.value);
    expect(z.toString(10)).toBe("−0.6299605249 + 1.091123636i…");
    expect(isReal(complex(2, 0))).toBe(true);
  });
});
