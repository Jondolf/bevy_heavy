//! `bevy_heavy` provides mass property types and computations for the geometric primitives
//! in the [Bevy game engine](https://bevyengine.org).

#![warn(missing_docs)]

mod dim2;
pub use dim2::MassProperties2d;

/// Returns the multiplicative inverse `1.0 / value` if `value` is greater than zero,
/// and `0.0` otherwise.
pub fn recip_or_zero(value: f32) -> f32 {
    let recip = value.recip();
    if recip.abs() > f32::EPSILON {
        recip
    } else {
        0.0
    }
}
