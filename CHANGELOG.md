# Changelog

## Unreleased

### Fixed

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
