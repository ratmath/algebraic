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
  complex, norm, abs, format,
} from "./src/arith.js";
export { factor, isIrreducible, factorModP } from "./src/factor.js";
export { NumberField, NumberFieldElement } from "./src/numberfield.js";
export { factorOver, rootsIn, extend, fieldContaining, splittingField } from "./src/fieldpoly.js";
export { primeField, galoisField, isProbablePrime } from "./src/finite.js";

/** Lower-level pieces: integer and rational polynomials, interval Krawczyk. */
export * as poly from "./src/poly.js";
export * as qpoly from "./src/qpoly.js";
export * as interval from "./src/interval.js";
export { cosPi, sinPi, tanPi, rootOfUnity, rootsOfUnity, cyclotomic, chebyshevT, cosMinimal } from "./src/trig.js";
export { continuedFraction, convergents, fromContinuedFraction, periodicContinuedFraction, pell, fundamentalUnit, fromQuadratic } from "./src/cf.js";

/** Polynomial toolkit over Q (friendly inputs, Core Rational outputs). */
export * as polys from "./src/polys.js";
export { resultant, discriminant, countRoots, interpolate, partialFractions } from "./src/polys.js";
export * as nt from "./src/nt.js";
export { factorInteger, crt, modPow, modInverse, jacobi, sqrtMod, primitiveRoot, multiplicativeOrder, discreteLogMod, totient } from "./src/nt.js";
export { galoisGroup, galoisAction, GaloisAction } from "./src/galois.js";
export { PermGroup, cycleString, fromCycles, cycleType } from "./src/permgroup.js";
export { lll, recognize } from "./src/lll.js";
