use std::ops::*;

use crate::{Mass, RecipOrZero};
use bevy_math::Vec2;

/// An error returned for an invalid angular inertia in 2D.
pub enum AngularInertia2dError {
    /// The mass is negative.
    Negative,
    /// The mass is NaN.
    Nan,
}

/// The moment of inertia of a 2D object. This represents the torque needed for a desired angular acceleration.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "bevy_reflect", derive(bevy_reflect::Reflect))]
#[cfg_attr(feature = "bevy_reflect", reflect(Debug, Default, PartialEq))]
#[cfg_attr(feature = "serialize", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    all(feature = "bevy_reflect", feature = "serialize"),
    reflect(Serialize, Deserialize)
)]
pub struct AngularInertia2d(f32);

impl Deref for AngularInertia2d {
    type Target = f32;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl AngularInertia2d {
    /// Zero angular inertia.
    pub const ZERO: Self = Self(0.0);

    /// An angular inertia of `1.0`.
    pub const ONE: Self = Self(1.0);

    /// Infinite angular inertia.
    pub const INFINITY: Self = Self(f32::INFINITY);

    /// Creates a new [`AngularInertia2d`] from the given angular inertia.
    ///
    /// # Panics
    ///
    /// Panics if the angular inertia is negative or NaN when `debug_assertions` are enabled.
    #[inline]
    pub fn new(angular_inertia: f32) -> Self {
        debug_assert!(
            angular_inertia >= 0.0,
            "inverse angular inertia must be positive or zero",
        );

        Self(angular_inertia)
    }

    /// Tries to create a new [`AngularInertia2d`] from the given angular inertia.
    ///
    /// Returns [`Err(AngularInertia2dError)`](AngularInertia2dError) if the angular inertia is negative.
    #[inline]
    pub fn try_new(angular_inertia: f32) -> Result<Self, AngularInertia2dError> {
        if angular_inertia < 0.0 {
            Err(AngularInertia2dError::Negative)
        } else if angular_inertia.is_nan() {
            Err(AngularInertia2dError::Nan)
        } else {
            Ok(Self(angular_inertia))
        }
    }

    /// Creates a new [`AngularInertia2d`] from the given inverse angular inertia.
    ///
    /// If the `inverse_angular_inertia` is infinite, the returned angular inertia will be zero.
    ///
    /// # Panics
    ///
    /// Panics if the inverse angular inertia is negative or NaN when `debug_assertions` are enabled.
    #[inline]
    pub fn from_inverse(inverse_angular_inertia: f32) -> Self {
        Self::new(inverse_angular_inertia.recip_or_zero())
    }

    /// Tries to create a new [`AngularInertia2d`] from the given inverse angular inertia.
    ///
    /// Returns [`Err(AngularInertia2dError)`](AngularInertia2dError) if the inverse angular inertia is negative.
    #[inline]
    pub fn try_from_inverse(inverse_angular_inertia: f32) -> Result<Self, AngularInertia2dError> {
        if inverse_angular_inertia < 0.0 {
            Err(AngularInertia2dError::Negative)
        } else if inverse_angular_inertia.is_nan() {
            Err(AngularInertia2dError::Nan)
        } else {
            Ok(Self::new(inverse_angular_inertia.recip_or_zero()))
        }
    }

    /// Returns the angular inertia.
    #[inline]
    pub fn value(self) -> f32 {
        self.0
    }

    /// Returns the inverse angular inertia.
    #[inline]
    pub fn inverse(self) -> Self {
        Self(self.0.recip_or_zero())
    }

    /// Sets the angular inertia.
    #[inline]
    pub fn set(&mut self, angular_inertia: AngularInertia2d) {
        *self = angular_inertia;
    }

    /// Computes the angular inertia shifted by the given offset, taking into account the given mass.
    #[inline]
    pub fn shifted(&self, mass: Mass, offset: Vec2) -> Self {
        if mass.value() > 0.0 && mass.is_finite() && offset != Vec2::ZERO {
            Self::new(self.0 + offset.length_squared() * mass.value())
        } else {
            Self::new(self.0)
        }
    }
}

impl From<f32> for AngularInertia2d {
    fn from(angular_inertia: f32) -> Self {
        Self::new(angular_inertia)
    }
}

impl From<AngularInertia2d> for f32 {
    fn from(angular_inertia: AngularInertia2d) -> Self {
        angular_inertia.0
    }
}

impl Add<AngularInertia2d> for AngularInertia2d {
    type Output = Self;

    #[inline]
    fn add(self, rhs: AngularInertia2d) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl AddAssign<AngularInertia2d> for AngularInertia2d {
    #[inline]
    fn add_assign(&mut self, rhs: AngularInertia2d) {
        self.0 += rhs.0;
    }
}

impl Mul<f32> for AngularInertia2d {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: f32) -> Self {
        Self(self.0 * rhs)
    }
}

impl Mul<AngularInertia2d> for f32 {
    type Output = AngularInertia2d;

    #[inline]
    fn mul(self, rhs: AngularInertia2d) -> AngularInertia2d {
        AngularInertia2d(self * rhs.0)
    }
}

impl MulAssign<f32> for AngularInertia2d {
    #[inline]
    fn mul_assign(&mut self, rhs: f32) {
        self.0 *= rhs;
    }
}

impl Div<f32> for AngularInertia2d {
    type Output = Self;

    #[inline]
    fn div(self, rhs: f32) -> Self {
        Self(self.0 / rhs)
    }
}

impl DivAssign<f32> for AngularInertia2d {
    #[inline]
    fn div_assign(&mut self, rhs: f32) {
        self.0 /= rhs;
    }
}

impl Mul<AngularInertia2d> for Mass {
    type Output = AngularInertia2d;

    #[inline]
    fn mul(self, angular_inertia: AngularInertia2d) -> AngularInertia2d {
        AngularInertia2d(*self * angular_inertia.0)
    }
}

impl MulAssign<Mass> for AngularInertia2d {
    #[inline]
    fn mul_assign(&mut self, mass: Mass) {
        self.0 *= *mass;
    }
}
