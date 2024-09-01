use bevy_math::{DVec2, Rot2, Vec2};

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

    /// Computes approximate mass properties from the given set of points.
    ///
    /// This can be used to estimate mass properties for arbitrary shapes
    /// by providing a set of sample points from inside the shape.
    ///
    /// The more points there are, and the more uniformly distributed they are,
    /// the more accurate the estimation will be.
    pub fn from_point_cloud(points: &[Vec2], mass: impl Into<Mass>) -> Self {
        let mass: Mass = mass.into();
        let points_recip = 1.0 / points.len() as f64;

        let center_of_mass =
            (points.iter().fold(DVec2::ZERO, |acc, p| acc + p.as_dvec2()) * points_recip).as_vec2();
        let unit_angular_inertia = points.iter().fold(0.0, |acc, p| {
            let r = p.as_dvec2().length_squared();
            acc + points_recip * r
        }) as f32;

        Self::new(mass, mass.value() * unit_angular_inertia, center_of_mass)
    }

    /// Returns the center of mass transformed into global space using the given translation and rotation.
    #[inline]
    pub fn global_center_of_mass(&self, translation: Vec2, rotation: impl Into<Rot2>) -> Vec2 {
        translation + rotation.into() * self.center_of_mass
    }

    /// Returns the mass properties transformed by the given translation and rotation.
    ///
    /// In 2D, this only transforms the center of mass.
    #[inline]
    pub fn transformed_by(mut self, translation: Vec2, rotation: impl Into<Rot2>) -> Self {
        self.transform_by(translation, rotation);
        self
    }

    /// Transforms the mass properties by the given translation and rotation.
    ///
    /// In 2D, this only transforms the center of mass.
    #[inline]
    pub fn transform_by(&mut self, translation: Vec2, rotation: impl Into<Rot2>) {
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
            * new_mass.value().recip_or_zero();

        // Compute the new principal angular inertia, taking the new center of mass into account.
        let i1 = self
            .angular_inertia
            .shifted(mass1, new_center_of_mass - self.center_of_mass);
        let i2 = other
            .angular_inertia
            .shifted(mass2, new_center_of_mass - other.center_of_mass);
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
            * new_mass.value().recip_or_zero();

        // Compute the new principal angular inertia, taking the new center of mass into account.
        let i1 = self
            .angular_inertia
            .shifted(mass1, new_center_of_mass - self.center_of_mass);
        let i2 = other
            .angular_inertia
            .shifted(mass2, new_center_of_mass - other.center_of_mass);
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

#[cfg(any(feature = "approx", test))]
impl approx::AbsDiffEq for MassProperties2d {
    type Epsilon = f32;
    fn default_epsilon() -> f32 {
        f32::EPSILON
    }
    fn abs_diff_eq(&self, other: &Self, epsilon: f32) -> bool {
        self.mass.abs_diff_eq(&other.mass, epsilon)
            && self
                .angular_inertia
                .abs_diff_eq(&other.angular_inertia, epsilon)
            && self
                .center_of_mass
                .abs_diff_eq(other.center_of_mass, epsilon)
    }
}

#[cfg(any(feature = "approx", test))]
impl approx::RelativeEq for MassProperties2d {
    fn default_max_relative() -> f32 {
        f32::EPSILON
    }
    fn relative_eq(&self, other: &Self, epsilon: f32, max_relative: f32) -> bool {
        self.mass.relative_eq(&other.mass, epsilon, max_relative)
            && self
                .angular_inertia
                .relative_eq(&other.angular_inertia, epsilon, max_relative)
            && self
                .center_of_mass
                .relative_eq(&other.center_of_mass, epsilon, max_relative)
    }
}

#[cfg(any(feature = "approx", test))]
impl approx::UlpsEq for MassProperties2d {
    fn default_max_ulps() -> u32 {
        4
    }
    fn ulps_eq(&self, other: &Self, epsilon: f32, max_ulps: u32) -> bool {
        self.mass.ulps_eq(&other.mass, epsilon, max_ulps)
            && self
                .angular_inertia
                .ulps_eq(&other.angular_inertia, epsilon, max_ulps)
            && self
                .center_of_mass
                .ulps_eq(&other.center_of_mass, epsilon, max_ulps)
    }
}
