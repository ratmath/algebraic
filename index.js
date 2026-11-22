// @ratmath/algebraic: exact algebraic numbers, polynomials and fields.
//
//   import { sqrt, add, fromRational, roots, factor, NumberField, primeField } from "@ratmath/algebraic";
//   add(sqrt(fromRational(2)), sqrt(fromRational(3))).toString(20)   // "3.14626436994197234232…"

export { AlgebraicReal, isolate, countRealRoots, squareFree, simplestBetween } from "./src/poly.js";
export { ComplexAlgebraic, certifyComplexRoots } from "./src/complex.js";
export {
  fromRational, rationalOf, roots, realRoots, minimalPolynomial, simplify,
  neg, add, sub, mul, inv, div, pow, sign, isZero, compare, equals,
  I, conjugate, re, im, nthRoot, sqrt, isReal, isComplex, MAX_DEGREE,
} from "./src/arith.js";
export { factor, isIrreducible, factorModP } from "./src/factor.js";
export { NumberField, NumberFieldElement } from "./src/numberfield.js";
export { primeField, galoisField, isProbablePrime } from "./src/finite.js";

/** Lower-level pieces: integer and rational polynomials, interval Krawczyk. */
export * as poly from "./src/poly.js";
export * as qpoly from "./src/qpoly.js";
export * as interval from "./src/interval.js";
