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

roots([-2, 0, 0, 1], { complex: true }); // [{ value, multiplicity }, …]: ∛2 and its two complex conjugates

const K = NumberField.quadratic(2);
K.of([1, 1]).inv().toString(); // "−1 + √2"
```

## Install

```sh
npm install @ratmath/algebraic      # or: bun add @ratmath/algebraic
```

ES modules only (Node 22+, Bun, browsers through a bundler). TypeScript
declarations ship with the package. The one dependency is
[`@ratmath/core`](https://ratmath.com/core/).

## What a number is

- **AlgebraicReal**: an irreducible (or at least square-free) integer
  polynomial and an interval holding exactly one of its real roots. `refine(bits)`
  shrinks the interval (quadratic interval refinement: the correct digits
  double per step once the root is close), `interval(bits)` gives it as Core Rationals,
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
sums (products) of conjugates of a and b; it is computed from the power sums
of the two polynomials (Newton's identities), then **factored**, so the result is put on its minimal
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

Polynomials over a field: `K.factor(p)` factors a polynomial with
coefficients in K into irreducibles over K (Trager's method), `K.roots(p)`
lists the roots that lie in K, and `K.extend(g)` adjoins a root of an
irreducible g, returning `{ field, embed, root }` with the bigger field as a
plain NumberField. `NumberField.containing([a, b, …])` builds one field
holding several exact numbers, and `NumberField.splitting(p)` the splitting
field of a rational polynomial with all its roots.

## Finite fields

`primeField(p)` (integers mod a prime) and `galoisField(p, k)` (GF(pᵏ), with
an irreducible modulus found automatically).

## Plugging into other libraries

`NumberField`, `primeField` and `galoisField` are field objects
`{ zero, one, add, sub, mul, div, neg, sign }`, the interface that
@ratmath/linalg's `field` module (and anything else field-generic) runs on.
@ratmath/linalg uses this package for eigenvalues, singular values and its
algebraic field.

## Complex numbers

`complex(x, y)`, `I()`, `abs`, `norm`, `re`, `im`, `conjugate`. `nthRoot(a, n)`
and `sqrt(a)` take principal roots of non-real numbers (the root with the
largest real part) and even roots of negative reals; odd roots of negative
reals stay real. `format(a)` prints exact forms where they are short:
`(1 + √5)/2`, `√(2 + √2)/2`, `(√6 − √2)/4`, `−1/2 + i√3/2`, `√3/2 + i/2`;
otherwise certified digits ending in "…".

## Angles and roots of unity

`cosPi(r)`, `sinPi(r)`, `tanPi(r)` for rational r (exact, no π involved;
`tanPi` reads its value off the real roots of Im (1 + it)^q),
`rootOfUnity(n, k)`, `rootsOfUnity(n)`, `cyclotomic(n)`, `chebyshevT(n)`,
`cosMinimal(n)`.

## Polynomials

`resultant`, `discriminant`, `countRoots(p, a, b)` (Sturm), `interpolate`,
`partialFractions(num, den)`, and the `polys` namespace (`add`, `mul`,
`divmod`, `gcd`, `gcdExt`, `compose`, `derivative`, `evaluate`, `fromRoots`,
`toInteger`, `format`). `galoisGroup(p)` names the Galois group of a
polynomial of degree ≤ 5 (resolvent cubic and Kappe–Warren for quartics;
sextic and quartic resolvents computed exactly from certified root
enclosures for quintics).

## Continued fractions, Pell, units

`continuedFraction(x, terms)` for any real algebraic number, `convergents`,
`periodicContinuedFraction(x)` for quadratic irrationals,
`fromContinuedFraction(terms, period)`, `pell(D)`, `fundamentalUnit(d)`.

## Number theory and recognition

`factorInteger` (Pollard rho), `crt`, `modPow`, `modInverse`, `jacobi`,
`sqrtMod`, `primitiveRoot`, `multiplicativeOrder`, `discreteLogMod`,
`totient`, `isProbablePrime` (Baillie–PSW: proven below 2⁶⁴, no known
counterexample), and the `nt` namespace. Finite fields gain `order`, `primitive`,
`log`, `elements` (and `sqrt` for F_p). `recognize("1.6180339887")` guesses
an algebraic number from its digits with LLL (`lll` is exported) and accepts
it only when the polynomial is much smaller than the digits given.

## Docs

`docs/` is a static site: an overview, a full reference, a free-form explore
page and ten topic pages with live demos. Edit the pages in `docs/src/` (full
HTML documents), not `docs/*.html`: `bun run docs:build` bundles the library
into `docs/assets/algebraic.js` and writes each page to `docs/` with an
editable JavaScript code box at the end of its article (seed code from the
page's `<template data-runbox>`), math in `\( … \)` / `\[ … \]` typeset with
KaTeX, and content-hash `?v=` stamps on asset links. The build runs every box
and fails if one throws (`--verbose` prints what each shows).
`docs/assets/runbox.js` and `scripts/docs-kit.js` are shared verbatim with the
core, analytic and linalg docs; keep the copies in sync. Serve `docs/` with any
static server or GitHub Pages.

## Limits

Resultants grow as (deg a)·(deg b), so long chains of operations on
high-degree numbers get slow (the factoring step dominates); `MAX_DEGREE`
(64) guards the size of any result's polynomial. Galois groups stop at degree
5. Field extensions stop at absolute degree 48. `format` writes square roots,
one level of nesting and pure n-th roots; other radicals print as digits.

## Develop

```sh
bun install
bun test
```

MIT licensed.
