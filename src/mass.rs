use std::ops::*;

use crate::RecipOrZero;

/// An error returned for an invalid mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MassError {
    /// The mass is negative.
    Negative,
    /// The mass is NaN.
    Nan,
}

/// The mass of an object. Must be positive or zero.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "bevy_reflect", derive(bevy_reflect::Reflect))]
#[cfg_attr(feature = "bevy_reflect", reflect(Debug, PartialEq))]
#[cfg_attr(feature = "serialize", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    all(feature = "bevy_reflect", feature = "serialize"),
    reflect(Serialize, Deserialize)
)]
pub struct Mass(f32);

impl Deref for Mass {
    type Target = f32;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Mass {
    /// Zero mass.
    pub const ZERO: Self = Self(0.0);

    /// A mass of `1.0`.
    pub const ONE: Self = Self(1.0);

    /// Infinite mass.
    pub const INFINITY: Self = Self(f32::INFINITY);

    /// Creates a new [`Mass`] from the given mass.
    ///
    /// # Panics
    ///
    /// Panics if the mass is negative or NaN when `debug_assertions` are enabled.
    #[inline]
    pub fn new(mass: f32) -> Self {
        debug_assert!(
            mass >= 0.0 && !mass.is_nan(),
            "mass must be positive or zero"
        );

        Self(mass)
    }

    /// Tries to create a new [`Mass`] from the given mass.
    ///
    /// # Errors
    ///
    /// Returns [`Err(MassError)`](MassError) if the mass is negative or NaN.
    #[inline]
    pub fn try_new(mass: f32) -> Result<Self, MassError> {
        if mass < 0.0 {
            Err(MassError::Negative)
        } else if mass.is_nan() {
            Err(MassError::Nan)
        } else {
            Ok(Self(mass))
        }
    }

    /// Returns the mass.
    #[inline]
    pub fn value(self) -> f32 {
        self.0
    }

    /// Returns the inverse mass.
    #[inline]
    pub fn inverse(self) -> Self {
        Self(self.0.recip_or_zero())
    }

    /// Sets the mass to the given value.
    #[inline]
    pub fn set(&mut self, mass: impl Into<Mass>) {
        *self = mass.into();
    }
}

impl From<f32> for Mass {
    #[inline]
    fn from(mass: f32) -> Self {
        Self::new(mass)
    }
}

impl From<Mass> for f32 {
    #[inline]
    fn from(mass: Mass) -> f32 {
        mass.0
    }
}

impl Add<Mass> for Mass {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Mass) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl AddAssign<Mass> for Mass {
    #[inline]
    fn add_assign(&mut self, rhs: Mass) {
        self.0 += rhs.0;
    }
}

impl Mul<f32> for Mass {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: f32) -> Self {
        Self(self.0 * rhs)
    }
}

impl Mul<Mass> for f32 {
    type Output = Mass;

    #[inline]
    fn mul(self, rhs: Mass) -> Mass {
        Mass(self * rhs.0)
    }
}

impl MulAssign<f32> for Mass {
    #[inline]
    fn mul_assign(&mut self, rhs: f32) {
        self.0 *= rhs;
    }
}

impl Div<f32> for Mass {
    type Output = Self;

    #[inline]
    fn div(self, rhs: f32) -> Self {
        Self(self.0 / rhs)
    }
}

impl DivAssign<f32> for Mass {
    #[inline]
    fn div_assign(&mut self, rhs: f32) {
        self.0 /= rhs;
    }
}
