# Changelog

## Unreleased

- `AlgebraicReal#toFixed(k, mode)` and `roundTo(k, mode)`: rounding to k
  decimal places with the Core rounding modes, refined until the rounding is
  decided, with ties of rational values decided exactly (Course Test roadmap,
  batch A6).
- Plain numbers work wherever an algebraic number is expected: `sqrt(48)`,
  `nthRoot(27n, 3)`, `add(sqrt(2), 1)`, `compare(sqrt(2), "1.5")` and the rest
  of the arithmetic read a bigint, number, numeric string, Core `Rational` or
  `Integer` as a rational algebraic number (they used to fail with
  "p.map is not a function"). The coercion is exported as `num`. Found by the
  Course Test (ratmath/shelf#7).
- `@ratmath/core` is now a peer dependency (`>=0.10.0 <0.12.0`): the package
  uses the program's own copy of Core instead of bringing its own, so Core 0.11
  works and there is never a second copy whose `Rational` is not the same class.

### Moved

- `continuedFraction`, `convergents`, `fromContinuedFraction`,
  `periodicContinuedFraction`, `pell`, `fundamentalUnit`, `fromQuadratic` and
  `crt` moved to `@ratmath/numtheory` (roadmap 4.3), where `crt` handles any
  moduli and returns `{ x, modulus }`. The explorer's `cf`, `period`,
  `convergents`, `pell` and `unit` commands go with them. The modular kernel
  (`modPow`, `modInverse`, `jacobi`, `sqrtMod`, `primitiveRoot`,
  `multiplicativeOrder`, `discreteLogMod`, `totient`, `factorInteger`) stays
  here for finite fields.

### Fixed

- `format` wrote the negative of a denested radical with the wrong sign on its
  second term (`−√3 + √2` for −√3 − √2).
- `isProbablePrime` is now Baillie–PSW (a base-2 strong test plus a strong
  Lucas test) and no longer accepts the composites that pass every one of the
  first twelve prime bases. It accepts Numbers as well as BigInts.

### Added

- `galoisGroup` covers degrees 6 and 7 (all 16 and 7 transitive groups).
- `galoisAction(p)`: the Galois group as permutations of the numbered roots,
  with `fixedField(H)`, `correspondence()` (every subgroup class with the field
  it fixes), `conjugation`, `frobenius(n)` (factor patterns modulo primes
  against the group's cycle types) and `derivedSeries()`.
- `PermGroup`, `cycleString`, `fromCycles`, `cycleType` for small permutation
  groups.
- The Galois docs page has an explorer: the group moving the roots, the
  subgroup diagram with fixed fields, and the statistics modulo primes.
- Factoring over a number field: `K.factor(poly)`, `K.roots(poly)`,
  `factorOver`, `rootsIn` (Trager's method).
- Field extensions: `K.extend(poly)` gives K(β) as one `NumberField`,
  `NumberField.containing([a, b, …])` gives a field holding several algebraic
  numbers, and `NumberField.splitting(poly)` gives the splitting field with
  every root as an element.
- Faster arithmetic: sums and products find their polynomial from power sums
  when the degrees are coprime, and root refinement is quadratic.
- Pure roots print as radicals: `∛2`, `∜5`, `1 + ∛2`.
- `tanPi` is computed directly.
- `interpolate` reports duplicate x values.
- TypeScript declarations (`index.d.ts`), a CI workflow, and docs pages that
  link to GitHub, npm and the sibling libraries.
- Requires `@ratmath/core` ^0.10.0.

## 0.1.0

First version: real and complex algebraic numbers, exact trigonometry,
factoring over Q, Galois groups, continued fractions, number fields and finite
fields.
