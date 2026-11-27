// Type declarations for @ratmath/algebraic.
import type { Rational, Integer } from "@ratmath/core";

/** A rational as an exact pair [numerator, denominator], denominator > 0. */
export type Pair = [bigint, bigint];
/** Anything read as an exact rational: a pair, a Core Rational, a bigint, an integer number or a string such as "3/4". */
export type RationalLike = Pair | Rational | Integer | bigint | number | string;
/** An integer polynomial, constant term first. */
export type IntPoly = bigint[];
/** A polynomial with rational coefficients, constant term first. */
export type PolyLike = RationalLike[];
/** An interval of pairs [lo, hi]. */
export type PairInterval = [Pair, Pair];

/** A real algebraic number: a root of `poly` isolated in [lo, hi]. */
export class AlgebraicReal {
  constructor(poly: IntPoly, lo: Pair, hi: Pair);
  static rational(pair: Pair): AlgebraicReal;
  poly: IntPoly;
  lo: Pair;
  hi: Pair;
  /** Shrink the isolating interval below width 2^-bits (quadratic interval refinement). */
  refine(bits?: number): this;
  /** The exact value as a pair when rational, else null. */
  readonly rational: Pair | null;
  readonly isRational: boolean;
  interval(bits?: number): [Rational, Rational];
  toRational(): Rational | null;
  toNumber(): number;
  /** Decimal digits, the last one certain up to one unit; ends in "…" when irrational. */
  toString(digits?: number): string;
  /** Rounded to k decimal places with a Core rounding mode (default "half-even"), decided exactly. */
  roundTo(k?: number, mode?: RoundingMode): Rational;
  /** k decimal places as a string, rounded with roundTo: every digit is certain. */
  toFixed(k?: number, mode?: RoundingMode): string;
}
export type RoundingMode = "half-even" | "half-up" | "toward-zero" | "floor" | "ceil";

/** A non-real algebraic number: a root of `poly` in a certified box re × im. */
export class ComplexAlgebraic {
  constructor(poly: IntPoly, re: PairInterval, im: PairInterval);
  poly: IntPoly;
  re: PairInterval;
  im: PairInterval;
  conjugate(): ComplexAlgebraic;
  refine(bits?: number): this;
  box(bits?: number): { re: [Rational, Rational]; im: [Rational, Rational] };
  toNumbers(): { re: number; im: number };
  toDecimal(digits?: number): string;
  toString(digits?: number): string;
  imSign(): number;
}

export type Algebraic = AlgebraicReal | ComplexAlgebraic;
/** An algebraic number, or a plain number read as a rational one. */
export type AlgebraicLike = Algebraic | RationalLike;

export const MAX_DEGREE: number;
export function fromRational(x: RationalLike): AlgebraicReal;
/** a as an algebraic number: plain numbers become rational AlgebraicReals. */
export function num(a: AlgebraicLike): Algebraic;
export function rationalOf(a: Algebraic): Rational | null;
export function roots(coeffs: PolyLike, options?: { complex?: boolean }): Array<{ value: Algebraic; multiplicity: number }>;
export function realRoots(coeffs: PolyLike): Array<{ value: AlgebraicReal; multiplicity: number }>;
export function minimalPolynomial(a: AlgebraicLike): IntPoly;
export function simplify<T extends Algebraic>(a: T): T;
export function neg<T extends Algebraic>(a: T): T;
export function add(a: AlgebraicLike, b: AlgebraicLike): Algebraic;
export function sub(a: AlgebraicLike, b: AlgebraicLike): Algebraic;
export function mul(a: AlgebraicLike, b: AlgebraicLike): Algebraic;
export function div(a: AlgebraicLike, b: AlgebraicLike): Algebraic;
export function inv(a: AlgebraicLike): Algebraic;
export function pow(a: AlgebraicLike, k: number | RationalLike): Algebraic;
export function sign(a: AlgebraicReal | RationalLike): -1 | 0 | 1;
export function isZero(a: AlgebraicLike): boolean;
export function compare(a: AlgebraicReal | RationalLike, b: AlgebraicReal | RationalLike): -1 | 0 | 1;
export function equals(a: AlgebraicLike, b: AlgebraicLike): boolean;
export function I(): ComplexAlgebraic;
export function conjugate<T extends Algebraic>(a: T): T;
export function re(a: AlgebraicLike): AlgebraicReal;
export function im(a: AlgebraicLike): AlgebraicReal;
export function nthRoot(a: AlgebraicLike, n: number): Algebraic;
export function sqrt(a: AlgebraicLike): Algebraic;
export function isReal(a: unknown): a is AlgebraicReal;
export function isComplex(a: unknown): a is ComplexAlgebraic;
export function complex(x: Algebraic | RationalLike, y?: Algebraic | RationalLike): Algebraic;
export function norm(a: AlgebraicLike): AlgebraicReal;
export function abs(a: AlgebraicLike): AlgebraicReal;
/** Exact text where one exists (3/4, (1 + √5)/2, √(2 + √2)/2, ∛2, 1 + i√3), else certified decimals ending in "…". */
export function format(a: AlgebraicLike, options?: { digits?: number }): string;

export function isolate(p: IntPoly): AlgebraicReal[];
export function countRealRoots(p: IntPoly): number;
export function squareFree(p: IntPoly): Array<{ factor: IntPoly; multiplicity: number }>;
export function simplestBetween(a: Pair, b: Pair): Pair;
export function certifyComplexRoots(p: IntPoly, count: number): ComplexAlgebraic[] | null;

export interface Factorization { content: bigint; factors: Array<{ factor: IntPoly; multiplicity: number }>; }
export function factor(p: IntPoly): Factorization;
export function isIrreducible(p: IntPoly): boolean;
export function factorModP(f: IntPoly, p: bigint): Array<{ factor: bigint[]; multiplicity: number }> | bigint[][];

// ---- number fields

export class NumberField {
  constructor(minpoly: PolyLike, options?: { name?: string; root?: number | "positive" | Algebraic | (() => Algebraic); check?: boolean });
  static quadratic(d: number | bigint): NumberField;
  static gaussian(): NumberField;
  /** One field holding several exact numbers. */
  static containing(numbers: Algebraic[], options?: { name?: string }): { field: NumberField; elements: NumberFieldElement[] };
  /** The splitting field of a rational polynomial. */
  static splitting(coeffs: PolyLike, options?: { name?: string }): { field: NumberField; roots: NumberFieldElement[]; degree: number };
  /** The monic defining polynomial (pairs, constant first). */
  m: Pair[];
  /** Degree over Q. */
  n: number;
  name: string;
  zero: NumberFieldElement;
  one: NumberFieldElement;
  /** The generator α. */
  gen: NumberFieldElement;
  /** The embedding of α. */
  readonly root: Algebraic;
  of(x: NumberFieldElement | RationalLike | RationalLike[]): NumberFieldElement;
  reduce(coeffs: Pair[]): NumberFieldElement;
  add(a: FieldArg, b: FieldArg): NumberFieldElement;
  sub(a: FieldArg, b: FieldArg): NumberFieldElement;
  mul(a: FieldArg, b: FieldArg): NumberFieldElement;
  div(a: FieldArg, b: FieldArg): NumberFieldElement;
  neg(a: FieldArg): NumberFieldElement;
  inv(a: FieldArg): NumberFieldElement;
  sign(a: FieldArg): 0 | 1;
  /** Factor a polynomial with coefficients in this field into irreducibles over it (Trager's method). */
  factor(coeffs: FieldArg[]): FieldFactorization;
  /** The roots of a polynomial that lie in this field. */
  roots(coeffs: FieldArg[]): Array<{ value: NumberFieldElement; multiplicity: number }>;
  /** Adjoin a root of a polynomial that is irreducible over this field. */
  extend(coeffs: FieldArg[], options?: { name?: string; pick?: Algebraic }): FieldExtension;
  toString(): string;
}
export type FieldArg = NumberFieldElement | RationalLike;
export interface FieldFactorization { unit: NumberFieldElement; factors: Array<{ factor: NumberFieldElement[]; multiplicity: number }>; }
export interface FieldExtension { field: NumberField; embed(x: FieldArg): NumberFieldElement; root: NumberFieldElement; }

export class NumberFieldElement {
  constructor(field: NumberField, coeffs: Pair[]);
  field: NumberField;
  /** Coefficients in 1, α, α², … (pairs, trimmed). */
  c: Pair[];
  add(b: FieldArg): NumberFieldElement;
  sub(b: FieldArg): NumberFieldElement;
  neg(): NumberFieldElement;
  mul(b: FieldArg): NumberFieldElement;
  div(b: FieldArg): NumberFieldElement;
  inv(): NumberFieldElement;
  pow(k: number): NumberFieldElement;
  isZero(): boolean;
  equals(b: FieldArg): boolean;
  norm(): Rational;
  trace(): Rational;
  coefficients(): Rational[];
  toRational(): Rational | null;
  toAlgebraic(): Algebraic;
  toString(): string;
}

export function factorOver(K: NumberField, coeffs: FieldArg[]): FieldFactorization;
export function rootsIn(K: NumberField, coeffs: FieldArg[]): Array<{ value: NumberFieldElement; multiplicity: number }>;
export function extend(K: NumberField, coeffs: FieldArg[], options?: { name?: string; pick?: Algebraic }): FieldExtension;
export function fieldContaining(numbers: Algebraic[], options?: { name?: string }): { field: NumberField; elements: NumberFieldElement[] };
export function splittingField(coeffs: PolyLike, options?: { name?: string }): { field: NumberField; roots: NumberFieldElement[]; degree: number };

// ---- finite fields

/** The operations every field object offers (also what linalg's `field` module consumes). */
export interface Field<E> {
  zero: E;
  one: E;
  add(a: E, b: E): E;
  sub(a: E, b: E): E;
  mul(a: E, b: E): E;
  div(a: E, b: E): E;
  neg(a: E): E;
  inv(a: E): E;
  sign(a: E): number;
}
export interface FiniteField<E> extends Field<E> {
  p: bigint;
  size: bigint;
  pow(a: E, e: bigint | number): E;
  equals(a: E, b: E): boolean;
  order(a: E): bigint;
  primitive(): E;
  log(a: E, g?: E): bigint | null;
  elements(): Iterable<E>;
  format(a: E): string;
  from(x: unknown): E;
}
export interface PrimeField extends FiniteField<bigint> { sqrt(a: bigint): bigint | null; }
export interface GaloisField extends FiniteField<bigint[]> { k: number; modulus: bigint[]; gen: bigint[]; }
export function primeField(p: bigint | number): PrimeField;
export function galoisField(p: bigint | number, k: number, modulus?: bigint[] | null): GaloisField;
/** Primality by the Baillie–PSW test (proven below 2^64; no counterexample is known). */
export function isProbablePrime(n: bigint | number): boolean;

// ---- trigonometry, roots of unity

export function cosPi(r: RationalLike): AlgebraicReal;
export function sinPi(r: RationalLike): AlgebraicReal;
export function tanPi(r: RationalLike): AlgebraicReal;
export function rootOfUnity(n: number, k?: number): Algebraic;
export function rootsOfUnity(n: number): Algebraic[];
export function cyclotomic(n: number): IntPoly;
export function chebyshevT(n: number): IntPoly;
export function cosMinimal(n: number): IntPoly;

// ---- continued fractions


// ---- polynomial toolkit over Q

export function resultant(p: PolyLike, q: PolyLike): Rational;
export function discriminant(p: PolyLike): Rational;
export function countRoots(p: PolyLike, a?: RationalLike, b?: RationalLike): number;
export function interpolate(points: Array<[RationalLike, RationalLike]>): Rational[];
export interface PartialFractions { polynomial: Rational[]; terms: Array<{ factor: IntPoly; power: number; numerator: Rational[] }>; }
export function partialFractions(num: PolyLike, den: PolyLike): PartialFractions;

// ---- number theory

export function factorInteger(n: bigint | number): Array<{ prime: bigint; exponent: number }>;
export function modPow(a: bigint, e: bigint, m: bigint): bigint;
export function modInverse(a: bigint, m: bigint): bigint;
export function jacobi(a: bigint, n: bigint): -1 | 0 | 1;
export function sqrtMod(a: bigint, p: bigint): bigint | null;
export function primitiveRoot(p: bigint): bigint;
export function multiplicativeOrder(a: bigint, n: bigint): bigint;
export function discreteLogMod(g: bigint, h: bigint, p: bigint): bigint | null;
export function totient(n: bigint): bigint;

// ---- Galois groups, lattices

export interface GaloisGroup { degree: number; name: string; order: number; label: string; solvable: boolean; discriminant: Rational; discriminantIsSquare: boolean; }
/** The Galois group of a polynomial of degree ≤ 7; a reducible one gives the group of each factor. */
export function galoisGroup(coeffs: PolyLike): GaloisGroup;

/** A permutation of 0 … n−1: p[i] is the image of i. */
export type Permutation = number[];
/** Cycle notation counted from 1, "(1 2 3)(4 5)"; "()" is the identity. */
export function cycleString(p: Permutation): string;
export function fromCycles(text: string, n: number): Permutation;
/** The cycle lengths, longest first, fixed points included. */
export function cycleType(p: Permutation): number[];
export interface SubgroupClass { group: PermGroup; order: number; index: number; conjugates: number; normal: boolean; }
/** A permutation group stored by its elements (meant for groups of at most a few thousand). */
export class PermGroup {
  constructor(generators: Permutation[], degree?: number, options?: { maxOrder?: number });
  static symmetric(n: number): PermGroup;
  static alternating(n: number): PermGroup;
  static from(generators: Array<string | Permutation>, n: number): PermGroup;
  degree: number;
  generators: Permutation[];
  /** Every element, by increasing order of the element; the identity is first. */
  elements: Permutation[];
  order: number;
  readonly isTransitive: boolean;
  readonly isAbelian: boolean;
  readonly isCyclic: boolean;
  readonly isEven: boolean;
  readonly isSolvable: boolean;
  has(p: Permutation): boolean;
  contains(H: PermGroup): boolean;
  equals(H: PermGroup): boolean;
  orbits(): number[][];
  orbitsOnSets(k: number): number[][][];
  cycleTypes(): Array<{ type: number[]; count: number; share: number }>;
  subgroup(generators: Array<string | Permutation>): PermGroup;
  /** g·H·g⁻¹. */
  conjugate(H: PermGroup, g: Permutation): PermGroup;
  isNormal(H: PermGroup): boolean;
  /** One element of each left coset g·H. */
  cosets(H: PermGroup): Permutation[];
  derivedSubgroup(): PermGroup;
  derivedSeries(): PermGroup[];
  withSmallGenerators(): PermGroup;
  /** Every subgroup; a RangeError above maxOrder (default 400). */
  subgroups(options?: { maxOrder?: number }): PermGroup[];
  subgroupClasses(options?: { maxOrder?: number }): SubgroupClass[];
  toString(): string;
}
export interface FixedField {
  /** The index of the subgroup. */
  degree: number;
  /** The minimal polynomial of the generator: integers, constant first. */
  polynomial: bigint[];
  /** The generator as a sum of products of roots, such as "x₁x₃ + x₂x₄". */
  generator: string;
  /** Present when the generator is written in moved roots: "y = x + x²" or "y = 2x". */
  substitution?: string;
  normal: boolean;
}
export interface CorrespondenceRow {
  subgroup: PermGroup;
  order: number;
  index: number;
  conjugates: number;
  normal: boolean;
  /** null when the degree exceeds maxDegree. */
  field: FixedField | null;
  /** Rows whose subgroup contains a conjugate of this one. */
  above: number[];
}
/** The Galois group of an irreducible polynomial of degree ≤ 7 acting on its numbered roots. */
export class GaloisAction {
  constructor(coeffs: PolyLike);
  polynomial: bigint[];
  degree: number;
  /** The roots in the order the permutations use: real ones ascending, then the others. */
  roots: Algebraic[];
  group: PermGroup;
  name: string;
  label: string;
  order: number;
  solvable: boolean;
  discriminant: Rational;
  discriminantIsSquare: boolean;
  orbitsOnPairs: number[][][] | null;
  orbitsOnTriples: number[][][] | null;
  /** Complex conjugation as an element of the group. */
  readonly conjugation: Permutation;
  permutations(): string[];
  fixedField(H: PermGroup | string | Array<string | Permutation>): FixedField;
  correspondence(options?: { maxDegree?: number; maxOrder?: number }): CorrespondenceRow[];
  frobenius(count?: number): { primes: number; patterns: Array<{ type: number[]; elements: number; expected: number; primes: number; observed: number; first: number[] }> };
  /** Orders of G ⊇ G′ ⊇ G″ ⊇ …; ends in 1 exactly when the group is solvable. */
  derivedSeries(): number[];
  toString(): string;
}
export function galoisAction(coeffs: PolyLike): GaloisAction;
export function lll(rows: bigint[][], delta?: Pair): bigint[][];
export interface Recognized { value: AlgebraicReal; poly: IntPoly; degree: number; height: bigint; digitsUsed: number; digitsNeeded: number; }
export function recognize(x: string | number, options?: { maxDegree?: number; margin?: number }): Recognized | null;

// ---- lower-level namespaces (see the reference page for their members)

export const poly: Record<string, any>;
export const qpoly: Record<string, any>;
export const interval: Record<string, any>;
export const polys: Record<string, any>;
export const nt: Record<string, any>;
