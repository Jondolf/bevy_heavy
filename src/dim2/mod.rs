use bevy_math::{Rot2, Vec2};

mod angular_inertia;
pub use angular_inertia::AngularInertia2d;

use crate::Mass;

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

/// The [mass], [angular inertia], and local center of mass of an object in 2D space.
///
/// [mass]: Mass
/// [angular inertia]: AngularInertia2d
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "bevy_reflect", derive(bevy_reflect::Reflect))]
#[cfg_attr(feature = "serialize", derive(serde::Serialize, serde::Deserialize))]
pub struct MassProperties2d {
    /// The mass.
    pub mass: Mass,
    /// The angular inertia along the principal axis.
    pub angular_inertia: AngularInertia2d,
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
        mass: Mass::ZERO,
        angular_inertia: AngularInertia2d::ZERO,
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
            mass: mass.into(),
            angular_inertia: angular_inertia.into(),
            center_of_mass,
        }
    }

    /// Returns the center of mass transformed into global space using the given translation and rotation.
    #[inline]
    pub fn global_center_of_mass(&self, translation: Vec2, rotation: Rot2) -> Vec2 {
        translation + rotation * self.center_of_mass
    }

    /// Returns the mass properties transformed by the given translation and rotation.
    ///
    /// In 2D, this only transforms the center of mass.
    #[inline]
    pub fn transformed_by(mut self, translation: Vec2, rotation: Rot2) -> Self {
        self.transform_by(translation, rotation);
        self
    }

    /// Transforms the mass properties by the given translation and rotation.
    ///
    /// In 2D, this only transforms the center of mass.
    #[inline]
    pub fn transform_by(&mut self, translation: Vec2, rotation: Rot2) {
        self.center_of_mass = self.global_center_of_mass(translation, rotation);
    }

    /// Returns the mass propeorties with the inverse of mass and angular inertia.
    pub fn inverse(&self) -> Self {
        Self {
            mass: self.mass.inverse(),
            angular_inertia: self.angular_inertia.inverse(),
            center_of_mass: self.center_of_mass,
        }
    }

    /// Sets the mass to the given `new_mass`. This also affects the angular inertia.
    #[inline]
    pub fn set_mass(&mut self, new_mass: impl Into<Mass>) {
        let new_mass = new_mass.into();

        // Adjust angular inertia based on new mass.
        self.angular_inertia /= self.mass;
        self.angular_inertia *= new_mass;

        self.mass = new_mass;
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

        let mass1 = self.mass;
        let mass2 = other.mass;
        let new_mass = mass1 + mass2;

        // The new center of mass is the weighted average of the centers of masses of `self` and `other`.
        let new_center_of_mass = (self.center_of_mass * mass1.value()
            + other.center_of_mass * mass2.value())
            / new_mass.value();

        // Compute the new principal angular inertia, taking the new center of mass into account.
        let i1 = self
            .angular_inertia
            .shifted(self.mass, new_center_of_mass - self.center_of_mass);
        let i2 = self
            .angular_inertia
            .shifted(self.mass, new_center_of_mass - other.center_of_mass);
        let new_angular_inertia = i1 + i2;

        Self {
            mass: new_mass,
            angular_inertia: new_angular_inertia,
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

        let mass1 = self.mass;
        let mass2 = other.mass;
        let Ok(new_mass) = Mass::try_new(mass1.value() - mass2.value()) else {
            return Self {
                center_of_mass: self.center_of_mass,
                ..Self::ZERO
            };
        };

        // The new center of mass is the negated weighted average of the centers of masses of `self` and `other`.
        let new_center_of_mass = (self.center_of_mass * mass1.value()
            - other.center_of_mass * mass2.value())
            / new_mass.value();

        // Compute the new principal angular inertia, taking the new center of mass into account.
        let i1 = self
            .angular_inertia
            .shifted(mass1, new_center_of_mass - self.center_of_mass);
        let i2 = self
            .angular_inertia
            .shifted(mass1, new_center_of_mass - other.center_of_mass);
        let new_angular_inertia =
            AngularInertia2d::try_new(i1.value() - i2.value()).unwrap_or_default();

        Self {
            mass: new_mass,
            angular_inertia: new_angular_inertia,
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
