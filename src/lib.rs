//! `bevy_heavy` is a crate for computing mass properties (mass, angular inertia, and center of mass)
//! for the [geometric primitives] in the [Bevy game engine][Bevy]. This is typically required
//! for things like physics simulations.
//!
//! [geometric primitives]: https://docs.rs/bevy/latest/bevy/math/primitives/index.html
//! [Bevy]: https://bevyengine.org
//!
//! ## Example
//!
//! ```
//! use bevy_heavy::{ComputeMassProperties2d, MassProperties2d};
//! use bevy_math::{primitives::Rectangle, Vec2};
//!
//! let rectangle = Rectangle::new(2.0, 1.0);
//! let density = 2.0;
//!
//! // You can compute mass properties individually.
//! let mass = rectangle.mass(density);
//! let angular_inertia = rectangle.angular_inertia(mass);
//! let center_of_mass = rectangle.center_of_mass();
//!
//! // You can also compute all mass properties at once, returning `MassProperties2d`.
//! // This can be more efficient when more than one property is needed.
//! let mass_props = rectangle.mass_properties(density);
//!
//! // Mass properties have several helpers.
//! let shifted_inertia = mass_props.shifted_angular_inertia(Vec2::new(-3.5, 1.0));
//! let global_center_of_mass = mass_props.global_center_of_mass(Vec2::new(5.0, 7.5), 0.0);
//!
//! // You can also add and subtract mass properties.
//! let mass_props_2 = MassProperties2d::new(mass, angular_inertia, Vec2::new(0.0, 1.0));
//! let sum = mass_props + mass_props_2;
//! approx::assert_relative_eq!(sum - mass_props_2, mass_props);
//! ```

#![warn(missing_docs)]

#[cfg(feature = "2d")]
mod dim2;
#[cfg(feature = "3d")]
mod dim3;
mod math_ext;

#[cfg(feature = "2d")]
pub use dim2::{ComputeMassProperties2d, MassProperties2d};
#[cfg(feature = "3d")]
pub use dim3::{
    AngularInertiaTensor, AngularInertiaTensorError, ComputeMassProperties3d, MassProperties3d,
    SymmetricEigen3,
};
pub use math_ext::{MatExt, RecipOrZero};
