use bevy_math::{Rot2, Vec2};

mod angular_inertia;
pub use angular_inertia::AngularInertia2d;

use crate::{Mass, RecipOrZero};

/// [`ComputeMassProperties2d`] implementations for 2D geometric primitives.
mod impls;

/// A trait for computing [`MassProperties2d`] for 2D objects.
pub trait ComputeMassProperties2d {
    /// Computes the mass of the object with a given `density`.
    fn mass(&self, density: f32) -> Mass;

    /// Computes the angular inertia corresponding to a unit mass.
    #[doc(alias = "unit_moment_of_inertia")]
    fn unit_angular_inertia(&self) -> AngularInertia2d;

    /// Computes the angular inertia corresponding to the given `mass`.
    ///
    /// Equivalent to `mass * shape.unit_angular_inertia()`.
    #[inline]
    #[doc(alias = "moment_of_inertia")]
    fn angular_inertia(&self, mass: impl Into<Mass>) -> AngularInertia2d {
        mass.into() * self.unit_angular_inertia()
    }

    /// Computes the local center of mass.
    fn center_of_mass(&self) -> Vec2;

    /// Computes the [`MassProperties2d`] with a given `density`.
    #[inline]
    fn mass_properties(&self, density: f32) -> MassProperties2d {
        let mass = self.mass(density);
        MassProperties2d::new(mass, self.angular_inertia(mass), self.center_of_mass())
    }
}

/// Mass properties in 2D.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "bevy_reflect", derive(bevy_reflect::Reflect))]
#[cfg_attr(feature = "serialize", derive(serde::Serialize, serde::Deserialize))]
pub struct MassProperties2d {
    /// The inverse mass, `1.0 / mass`.
    pub inverse_mass: Mass,
    /// The inverse of the angular inertia along the principal axis, `(1.0 / inertia).sqrt()`.
    pub inverse_angular_inertia_sqrt: AngularInertia2d,
    /// The local center of mass.
    pub center_of_mass: Vec2,
}

impl Default for MassProperties2d {
    /// Returns the default [`MassProperties2d`], with zero mass and angular inertia.
    fn default() -> Self {
        Self::ZERO
    }
}

impl MassProperties2d {
    /// Zero mass and angular inertia.
    pub const ZERO: Self = Self {
        inverse_mass: Mass::INFINITY,
        inverse_angular_inertia_sqrt: AngularInertia2d::ZERO,
        center_of_mass: Vec2::ZERO,
    };

    /// Creates a new [`MassProperties2d`] from a given mass, principal angular inertia,
    /// and center of mass in local space.
    #[inline]
    pub fn new(
        mass: impl Into<Mass>,
        angular_inertia: impl Into<AngularInertia2d>,
        center_of_mass: Vec2,
    ) -> Self {
        Self {
            center_of_mass,
            inverse_mass: Mass::new(mass.into().recip_or_zero()),
            inverse_angular_inertia_sqrt: AngularInertia2d::new(
                angular_inertia.into().sqrt().recip_or_zero(),
            ),
        }
    }

    /// Returns the mass.
    #[inline]
    pub fn mass(&self) -> Mass {
        self.inverse_mass.inverse()
    }

    /// Returns the principal angular inerta.
    #[inline]
    pub fn angular_inertia(&self) -> AngularInertia2d {
        AngularInertia2d::new(self.inverse_angular_inertia_sqrt.powi(2).recip_or_zero())
    }

    /// Returns the center of mass transformed into global space using `translation` and `rotation`.
    #[inline]
    pub fn global_center_of_mass(&self, translation: Vec2, rotation: Rot2) -> Vec2 {
        translation + rotation * self.center_of_mass
    }

    /// Computes the principal angular inertia at a given `offset`.
    #[inline]
    pub fn shifted_angular_inertia(&self, offset: Vec2) -> AngularInertia2d {
        AngularInertia2d::new(self.inverse_angular_inertia_sqrt.powi(2).recip_or_zero())
            .shifted(self.inverse_mass.inverse(), offset)
    }

    /// Returns the mass properties transformed by `translation` and `rotation`.
    ///
    /// In 2D, this only transforms the center of mass.
    #[inline]
    pub fn transformed_by(mut self, translation: Vec2, rotation: Rot2) -> Self {
        self.transform_by(translation, rotation);
        self
    }

    /// Transforms the mass properties by `translation` and `rotation`.
    ///
    /// In 2D, this only transforms the center of mass.
    #[inline]
    pub fn transform_by(&mut self, translation: Vec2, rotation: Rot2) {
        self.center_of_mass = self.global_center_of_mass(translation, rotation);
    }

    /// Sets the mass to the given `new_mass`. This also affects the angular inertia.
    #[inline]
    pub fn set_mass(&mut self, new_mass: impl Into<Mass>) {
        let new_inverse_mass = new_mass.into().inverse();

        // Adjust angular inertia based on new mass.
        let old_mass = self.inverse_mass.recip_or_zero();
        self.inverse_angular_inertia_sqrt *= new_inverse_mass.sqrt() * old_mass.sqrt();

        self.inverse_mass = new_inverse_mass;
    }
}

impl std::ops::Add for MassProperties2d {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        if self == Self::ZERO {
            return other;
        } else if other == Self::ZERO {
            return self;
        }

        let mass1 = self.mass();
        let mass2 = other.mass();
        let new_inverse_mass = (mass1 + mass2).inverse();

        // The new center of mass is the weighted average of the centers of masses of `self` and `other`.
        let new_center_of_mass = (self.center_of_mass * mass1.value()
            + other.center_of_mass * mass2.value())
            * new_inverse_mass.value();

        // Compute the new principal angular inertia, taking the new center of mass into account.
        let inertia1 = self.shifted_angular_inertia(new_center_of_mass - self.center_of_mass);
        let inertia2 = self.shifted_angular_inertia(new_center_of_mass - other.center_of_mass);
        let new_inertia = inertia1 + inertia2;

        Self {
            inverse_mass: new_inverse_mass,
            inverse_angular_inertia_sqrt: AngularInertia2d::new(new_inertia.sqrt().recip_or_zero()),
            center_of_mass: new_center_of_mass,
        }
    }
}

impl std::ops::AddAssign for MassProperties2d {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl std::ops::Sub for MassProperties2d {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        if self == Self::ZERO || other == Self::ZERO {
            return self;
        }

        let mass1 = self.mass();
        let mass2 = other.mass();
        let new_mass = mass1.value() - mass2.value();

        let new_inverse_mass = if new_mass >= f32::EPSILON {
            (mass1 + mass2).inverse()
        } else {
            Mass::INFINITY
        };

        // The new center of mass is the negated weighted average of the centers of masses of `self` and `other`.
        let new_center_of_mass = (self.center_of_mass * mass1.value()
            - other.center_of_mass * mass2.value())
            * new_inverse_mass.value();

        // Compute the new principal angular inertia, taking the new center of mass into account.
        let inertia1 = self.shifted_angular_inertia(new_center_of_mass - self.center_of_mass);
        let inertia2 = self.shifted_angular_inertia(new_center_of_mass - other.center_of_mass);
        let new_inertia = inertia1.value() - inertia2.value();

        Self {
            inverse_mass: new_inverse_mass,
            inverse_angular_inertia_sqrt: AngularInertia2d::new(new_inertia.sqrt().recip_or_zero()),
            center_of_mass: new_center_of_mass,
        }
    }
}

impl std::ops::SubAssign for MassProperties2d {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}
