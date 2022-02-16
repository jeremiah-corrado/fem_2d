# Contributing to `fem_2d`

New features and new trait implementations are welcome! In either case, a new branch should be started with a descriptive name.

All tests (run with `cargo test`) must pass! New tests should be included with new features in the appropriate module.

Breaking changes to the public API are acceptable as this library is pre 1.0; however they should be avioded whenever possible, and should have strong justification when enacted.

## New `Trait` Implementations

`fem_2d` is pulbicly generic over `Integral`s and `ShapeFn`s which are generally problem specific. New implementations of either of these traits for alternate PDEs are welcome.
  
* `Integral` implementations should be added to their own module located under *[/src/integration/]*
* `ShapeFn` implementations should also be added to their own module located under *[/src/basis/]*

## New Features

If `fem_2d` is missing a feature that you would like to see, feel free to add it! Be sure that all new code is well documented and thoroughly tested. PR's should include plenty of description of- and justification for- the new feature. 

Take care to include new code in a sensible location within the existing module structure. The base-level modules: *domain*, *basis*, *integration*, and *linalg* provide a nice conseptual separation of the libraries funcitonality. New features will likely fit somewhere in one of these categories; however a new base-level module may be added with sufficient justification.

## Release Process

1. Implement and test a new feature, trait implementation, or bug-fix.
2. Send a PR to this repository.
3. Once merged, the version will be incremented by a Crate-maintainer, and the new version will be published to [crates.io](https://crates.io/)
