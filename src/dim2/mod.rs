use crate::recip_or_zero;
use bevy_math::{Mat2, Vec2};

/// [`ComputeMassProperties2d`] implementations for 2D geometric primitives.
mod impls;

/// A trait for computing [`MassProperties2d`] for 2D objects.
pub trait ComputeMassProperties2d {
    /// Computes the mass of the object with a given `density`.
    fn mass(&self, density: f32) -> f32;

    /// Computes the angular inertia corresponding to the given `mass`.
    #[doc(alias = "moment_of_inertia")]
    fn angular_inertia(&self, mass: f32) -> f32;

    /// Computes the local center of mass.
    fn center_of_mass(&self) -> Vec2;

    /// Computes the [`MassProperties2d`] with a given `density`.
    fn mass_properties(&self, density: f32) -> MassProperties2d {
        let mass = self.mass(density);
        MassProperties2d::new(mass, self.angular_inertia(mass), self.center_of_mass())
    }
}

/// Mass properties in 2D.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassProperties2d {
    /// The multiplicative inverse of mass, `1.0 / mass`.
    pub inverse_mass: f32,
    /// The square root of the multiplicative inverse of the angular inertia
    /// along the principal axis, `(1.0 / inertia).sqrt()`.
    pub inverse_angular_inertia_sqrt: f32,
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
        inverse_mass: 0.0,
        inverse_angular_inertia_sqrt: 0.0,
        center_of_mass: Vec2::ZERO,
    };

    /// Creates a new [`MassProperties2d`] from a given mass, principal angular inertia,
    /// and center of mass in local space.
    pub fn new(mass: f32, angular_inertia: f32, center_of_mass: Vec2) -> Self {
        Self {
            center_of_mass,
            inverse_mass: recip_or_zero(mass),
            inverse_angular_inertia_sqrt: recip_or_zero(angular_inertia.sqrt()),
        }
    }

    /// Returns the mass.
    pub fn mass(&self) -> f32 {
        recip_or_zero(self.inverse_mass)
    }

    /// Returns the principal angular inerta.
    pub fn angular_inertia(&self) -> f32 {
        recip_or_zero(self.inverse_angular_inertia_sqrt.powi(2))
    }

    /// Returns the center of mass transformed into global space using `translation` and `rotation`.
    pub fn global_center_of_mass(&self, translation: Vec2, rotation: f32) -> Vec2 {
        translation + Mat2::from_angle(rotation) * self.center_of_mass
    }

    /// Computes the principal angular inertia at a given `offset`.
    pub fn shifted_angular_inertia(&self, offset: Vec2) -> f32 {
        let inertia = recip_or_zero(self.inverse_angular_inertia_sqrt.powi(2));

        if self.inverse_mass > f32::EPSILON {
            let mass = 1.0 / self.inverse_mass;
            inertia + offset.length_squared() * mass
        } else {
            inertia
        }
    }

    /// Returns the mass properties transformed by `translation` and `rotation`.
    ///
    /// In 2D, this only transforms the center of mass.
    pub fn transformed_by(mut self, translation: Vec2, rotation: f32) -> Self {
        self.transform_by(translation, rotation);
        self
    }

    /// Transforms the mass properties by `translation` and `rotation`.
    ///
    /// In 2D, this only transforms the center of mass.
    pub fn transform_by(&mut self, translation: Vec2, rotation: f32) {
        self.center_of_mass = self.global_center_of_mass(translation, rotation);
    }

    /// Sets the mass to the given `new_mass`. This also affects the angular inertia.
    pub fn set_mass(&mut self, new_mass: f32) {
        let new_inverse_mass = recip_or_zero(new_mass);

        // Adjust angular inertia based on new mass.
        let old_mass = recip_or_zero(self.inverse_mass);
        self.inverse_angular_inertia_sqrt *= new_inverse_mass.sqrt() * old_mass.sqrt();

        self.inverse_mass = new_inverse_mass;
    }
}

impl std::ops::Add for MassProperties2d {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        if self == Self::ZERO {
            return other;
        } else if other == Self::ZERO {
            return self;
        }

        let mass1 = self.mass();
        let mass2 = other.mass();
        let new_inverse_mass = recip_or_zero(mass1 + mass2);

        // The new center of mass is the weighted average of the centers of masses of `self` and `other`.
        let new_center_of_mass =
            (self.center_of_mass * mass1 + other.center_of_mass * mass2) * new_inverse_mass;

        // Compute the new principal angular inertia, taking the new center of mass into account.
        let inertia1 = self.shifted_angular_inertia(new_center_of_mass - self.center_of_mass);
        let inertia2 = self.shifted_angular_inertia(new_center_of_mass - other.center_of_mass);
        let new_inertia = inertia1 + inertia2;

        Self {
            inverse_mass: new_inverse_mass,
            inverse_angular_inertia_sqrt: recip_or_zero(new_inertia.sqrt()),
            center_of_mass: new_center_of_mass,
        }
    }
}

impl std::ops::AddAssign for MassProperties2d {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl std::ops::Sub for MassProperties2d {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        if self == Self::ZERO || other == Self::ZERO {
            return self;
        }

        let mass1 = self.mass();
        let mass2 = other.mass();
        let new_mass = mass1 - mass2;

        let new_inverse_mass = if new_mass >= f32::EPSILON {
            recip_or_zero(mass1 + mass2)
        } else {
            0.0
        };

        // The new center of mass is the negated weighted average of the centers of masses of `self` and `other`.
        let new_center_of_mass =
            (self.center_of_mass * mass1 - other.center_of_mass * mass2) * new_inverse_mass;

        // Compute the new principal angular inertia, taking the new center of mass into account.
        let inertia1 = self.shifted_angular_inertia(new_center_of_mass - self.center_of_mass);
        let inertia2 = self.shifted_angular_inertia(new_center_of_mass - other.center_of_mass);
        let new_inertia = inertia1 - inertia2;

        Self {
            inverse_mass: new_inverse_mass,
            inverse_angular_inertia_sqrt: recip_or_zero(new_inertia.sqrt()),
            center_of_mass: new_center_of_mass,
        }
    }
}

impl std::ops::SubAssign for MassProperties2d {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}
