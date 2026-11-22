# @ratmath/algebraic

Exact algebraic numbers for RatMath: roots of polynomials held exactly, with
arithmetic, comparison and as many digits as you ask for. Also polynomial
factoring over the integers, number fields such as Q(√2) and Q(i), and finite
fields.

```js
import { fromRational, sqrt, add, mul, pow, compare, roots, NumberField } from "@ratmath/algebraic";

const r2 = sqrt(fromRational(2));
const r3 = sqrt(fromRational(3));
const s = add(r2, r3);
s.toString(20);           // "3.14626436994197234232…"
s.poly;                   // [1n, 0n, -10n, 0n, 1n]   (x⁴ − 10x² + 1, its minimal polynomial)
pow(s, 2).poly;           // [1n, -10n, 1n]            (5 + 2√6 has degree 2)
mul(r2, r2).toRational(); // 2, exactly
compare(r2, fromRational("7/5")); // 1

roots([-2, 0, 0, 1], { complex: true }); // ∛2 and its two complex conjugates

const K = NumberField.quadratic(2);
K.of([1, 1]).inv().toString(); // "−1 + √2"
```

## What a number is

- **AlgebraicReal**: an irreducible (or at least square-free) integer
  polynomial and an interval holding exactly one of its real roots. `refine(bits)`
  shrinks the interval, `interval(bits)` gives it as Core Rationals,
  `toRational()` is the exact value when rational (decided exactly),
  `toString(digits)` prints certain digits.
- **ComplexAlgebraic**: a polynomial and a rectangle in the complex plane
  certified (Krawczyk test on p(x + iy) = 0) to hold exactly one root.
  `refine(bits)`, `box(bits)`, `conjugate()`, `toNumbers()`.

## Arithmetic

`add`, `sub`, `mul`, `div`, `inv`, `neg`, `pow`, `sqrt` (√ of a negative real
is i·√), `nthRoot`, `I()`, `conjugate`, `re`, `im`, `sign`, `compare`, `equals`,
`isZero`, `minimalPolynomial`, `simplify`.

How: a ± b and a·b are roots of a resultant, the polynomial whose roots are all
sums (products) of conjugates of a and b; it is computed by evaluation and
interpolation, then **factored**, so the result is put on its minimal
polynomial and degrees stay small. The right root among the candidates is
picked by refining enclosures until exactly one fits. Inverses reverse the
coefficients; roots substitute xⁿ. Equality and zero tests are exact: an
irrational or non-real number is never 0, and rationality is decided exactly.

## Factoring

`factor(p)` → `{ content, factors: [{ factor, multiplicity }] }` for integer
polynomials (constant term first, BigInt or number coefficients), with each
factor primitive and irreducible over Q. Zassenhaus: square-free
decomposition, Cantor–Zassenhaus factoring modulo a small prime, Hensel lifting
past a Mignotte bound, then recombination. Also `isIrreducible`,
`factorModP(p, prime)`.

## Number fields

`new NumberField(minpoly, { name, root })`, `NumberField.quadratic(d)`,
`NumberField.gaussian()`. Elements (`K.of(rational)`, `K.of([c0, c1, …])`,
`K.gen`) support `add`, `sub`, `mul`, `div`, `inv`, `pow`, `equals`, `isZero`,
`norm`, `trace`, `coefficients`, `toRational`, `toString` ("1 + √2", "−3 + 4i")
and `toAlgebraic()`, the exact number through the field's embedding of α.
Arithmetic here is closed and cheap: everything stays a polynomial in α of
degree below deg m.

## Finite fields

`primeField(p)` (integers mod a prime) and `galoisField(p, k)` (GF(pᵏ), with
an irreducible modulus found automatically).

## Plugging into other libraries

`NumberField`, `primeField` and `galoisField` are field objects
`{ zero, one, add, sub, mul, div, neg, sign }`, the interface that
@ratmath/linalg's `field` module (and anything else field-generic) runs on.
@ratmath/linalg uses this package for eigenvalues, singular values and its
algebraic field.

## Limits

Complex square roots of non-real numbers are not supported yet. Resultants
grow as (deg a)·(deg b), so long chains of operations on high-degree numbers
get slow; `MAX_DEGREE` (64) guards the size of any result's polynomial.

## Develop

```sh
bun install
bun test
```

MIT licensed.
