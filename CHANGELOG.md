# Changelog

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). This file only reports changes that increase major and minor versions, as well as deprecations.

## 3.0.0 - 2021/10/21

[Full diff](https://framagit.org/atooms/postprocessing/-/compare/2.7.2...3.0.0)

### Changed
- Total and self intermediate scattering functions have the same time grid by default
- Privatize `FourierSpaceCorrelation.kvector` as `FourierSpaceCorrelation._kvectors`
- Adaptative bin width for `FourierSpaceCorrelation.kgrid`, it is automatically increased until all bins in `kgrid` are non-empty
- Rename some private variables in the `FourierSpaceCorrelation` hierarchy

### New
- Add `FourierSpaceCorrelation.kvectors` as a property to store the lists of k-vectors used for each entry in `FourierSpaceCorrelation.kgrid`
- Add `Partial.output_path` property
