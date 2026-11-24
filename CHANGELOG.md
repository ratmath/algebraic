# Changelog

## Unreleased

### Fixed

- `isProbablePrime` is now Baillie–PSW (a base-2 strong test plus a strong
  Lucas test) and no longer accepts the composites that pass every one of the
  first twelve prime bases. It accepts Numbers as well as BigInts.

### Added

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
