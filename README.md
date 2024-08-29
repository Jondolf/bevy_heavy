# `bevy_heavy`

`bevy_heavy` is a crate for computing mass properties (mass, angular inertia, and center of mass)
for the [geometric primitives] in the [Bevy game engine][Bevy]. This is typically required
for things like physics simulations.

> **Warning**: `bevy_heavy` is WIP and not well tested yet. I will release it once there are adequate tests and docs,
> and I've made sure it works integrated into a physics engine.

[geometric primitives]: https://docs.rs/bevy/latest/bevy/math/primitives/index.html
[Bevy]: https://bevyengine.org

## Features

- `MassProperties2d` and `MassProperties3d` structs containing the mass, angular inertia, and local center of mass
- Mass property computation for all of Bevy's [geometric primitives]
- Eigensolver for symmetric 3x3 matrices
- Support for `bevy_reflect` and `serde` through the `bevy_reflect` and `serialize` feature flags

## Example

```rust
use bevy_heavy::{ComputeMassProperties2d, MassProperties2d};
use bevy_math::primitives::Rectangle;

let rectangle = Rectangle::new(2.0, 1.0);
let density = 2.0;

// You can compute mass properties individually.
let mass = rectangle.mass(density);
let angular_inertia = rectangle.angular_inertia(mass);
let center_of_mass = rectangle.center_of_mass();

// You can also compute all mass properties at once, returning `MassProperties2d`.
// This can be more efficient when more than one property is needed.
let mass_props = rectangle.mass_properties(density);

// `MassProperties2d` has several helpers.
let shifted_inertia = mass_props.shifted_angular_inertia(Vec2::new(-3.5, 1.0));
let global_center_of_mass = mass_props.global_center_of_mass(Vec2::new(5.0, 7.5));

// You can also add and subtract mass properties.
let mass_props_2 = MassProperties2d::new(
    Mass::new(1.0),
    AngularInertia2d::new(0.5),
    Vec2::new(0.0, 1.0),
);
let sum = mass_props + mass_props_2;
assert_eq!(sum - mass_props_2, mass_props);
```

## License

`bevy_heavy` is free, open source, and permissively licensed! Except where noted (below and/or in individual files),
all code in this repository is dual-licensed under either:

- MIT License ([LICENSE-MIT](/LICENSE-MIT) or <http://opensource.org/licenses/MIT>)
- Apache License, Version 2.0 ([LICENSE-APACHE](/LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)

at your option. This dual-licensing approach is the de-facto standard in the Rust ecosystem,
and there are [very good reasons](https://github.com/bevyengine/bevy/issues/2373) to include both.
