use bevy_math::{Isometry3d, Mat3, Quat, Vec3};

mod angular_inertia;
pub use angular_inertia::AngularInertiaTensor;

/// [`ComputeMassProperties3d`] implementations for 3D geometric primitives.
mod impls;

mod eigen3;
pub use eigen3::SymmetricEigen3;

use crate::Mass;

/// A trait for computing [`MassProperties3d`] for 3D objects.
pub trait ComputeMassProperties3d {
    /// Computes the mass of the object with a given `density`.
    fn mass(&self, density: f32) -> Mass;

    /// Computes the principal angular inertia corresponding to a unit mass.
    #[doc(alias = "unit_principal_moment_of_inertia")]
    fn unit_principal_angular_inertia(&self) -> Vec3;

    /// Computes the principal angular inertia corresponding to the given `mass`.
    #[inline]
    #[doc(alias = "principal_moment_of_inertia")]
    fn principal_angular_inertia(&self, mass: impl Into<Mass>) -> Vec3 {
        let mass: Mass = mass.into();
        mass.value() * self.unit_principal_angular_inertia()
    }

    /// Computes the orientation of the inertial frame used by the principal axes of inertia in local space.
    ///
    /// For most primitive shapes, this returns an identity quaternion, which means that the principal axes
    /// are aligned with the object's XYZ axes.
    #[inline]
    fn local_inertial_frame(&self) -> Quat {
        Quat::IDENTITY
    }

    /// Computes the 3x3 [`AngularInertiaTensor`] corresponding to a unit mass.
    #[inline]
    fn unit_angular_inertia_tensor(&self) -> AngularInertiaTensor {
        AngularInertiaTensor::new_with_local_frame(
            self.principal_angular_inertia(1.0),
            self.local_inertial_frame(),
        )
    }

    /// Computes the 3x3 [`AngularInertiaTensor`] corresponding to the given `mass`.
    #[inline]
    fn angular_inertia_tensor(&self, mass: impl Into<Mass>) -> AngularInertiaTensor {
        AngularInertiaTensor::new_with_local_frame(
            self.principal_angular_inertia(mass),
            self.local_inertial_frame(),
        )
    }

    /// Computes the local center of mass.
    fn center_of_mass(&self) -> Vec3;

    /// Computes the [`MassProperties3d`] with a given `density`.
    #[inline]
    fn mass_properties(&self, density: f32) -> MassProperties3d {
        let mass = self.mass(density);
        MassProperties3d::new_with_local_frame(
            mass,
            self.principal_angular_inertia(mass),
            self.local_inertial_frame(),
            self.center_of_mass(),
        )
    }
}

/// The [mass], [angular inertia], and local center of mass of an object in 3D space.
///
/// [mass]: Mass
/// [angular inertia]: AngularInertiaTensor
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "bevy_reflect", derive(bevy_reflect::Reflect))]
#[cfg_attr(feature = "serialize", derive(serde::Serialize, serde::Deserialize))]
pub struct MassProperties3d {
    /// The mass.
    pub mass: Mass,
    /// The angular inertia along the principal axes.
    pub principal_angular_inertia: Vec3,
    /// The orientation of the local inertial frame, defining the principal axes.
    pub local_inertial_frame: Quat,
    /// The local center of mass.
    pub center_of_mass: Vec3,
}

impl Default for MassProperties3d {
    /// Returns the default [`MassProperties3d`], with zero mass and angular inertia.
    fn default() -> Self {
        Self::ZERO
    }
}

impl MassProperties3d {
    /// Zero mass and angular inertia.
    pub const ZERO: Self = Self {
        mass: Mass::ZERO,
        principal_angular_inertia: Vec3::ZERO,
        local_inertial_frame: Quat::IDENTITY,
        center_of_mass: Vec3::ZERO,
    };

    /// Creates a new [`MassProperties3d`] from a given mass, principal angular inertia,
    /// and center of mass in local space.
    #[inline]
    pub fn new(mass: impl Into<Mass>, angular_inertia: Vec3, center_of_mass: Vec3) -> Self {
        Self::new_with_local_frame(mass, angular_inertia, Quat::IDENTITY, center_of_mass)
    }

    /// Creates a new [`MassProperties3d`] from a given mass, principal angular inertia,
    /// local inertial frame, and center of mass in local space.
    ///
    /// The principal angular inertia is the angular inertia along the coordinate axes defined
    /// by the `local_inertial_frame`, expressed in local space.
    #[inline]
    pub fn new_with_local_frame(
        mass: impl Into<Mass>,
        principal_angular_inertia: Vec3,
        local_inertial_frame: Quat,
        center_of_mass: Vec3,
    ) -> Self {
        Self {
            mass: mass.into(),
            principal_angular_inertia,
            local_inertial_frame,
            center_of_mass,
        }
    }

    /// Creates a new [`MassProperties3d`] from a given mass, angular inertia tensor,
    /// and center of mass in local space.
    ///
    /// The angular inertia tensor will be diagonalized in order to extract the principal inertia
    /// values and local principal inertia frme.
    #[inline]
    pub fn new_with_angular_inertia_tensor(
        mass: impl Into<Mass>,
        tensor: impl Into<AngularInertiaTensor>,
        center_of_mass: Vec3,
    ) -> Self {
        let (principal_inertia, principal_inertia_local_frame) =
            tensor.into().principal_angular_inertia_with_local_frame();

        Self::new_with_local_frame(
            mass,
            principal_inertia,
            principal_inertia_local_frame,
            center_of_mass,
        )
    }

    /// Returns the center of mass transformed into global space using the given [isometry].
    ///
    /// [isometry]: Isometry3d
    #[inline]
    pub fn global_center_of_mass(&self, isometry: Isometry3d) -> Vec3 {
        isometry * self.center_of_mass
    }

    /// Computes the world-space angular inertia tensor from the principal inertia.
    #[inline]
    pub fn angular_inertia_tensor(&self) -> AngularInertiaTensor {
        AngularInertiaTensor::new_with_local_frame(
            self.principal_angular_inertia,
            self.local_inertial_frame,
        )
    }

    /// Computes the inertia tensor at a given `offset`.
    #[inline]
    pub fn shifted_angular_inertia_tensor(&self, offset: Vec3) -> AngularInertiaTensor {
        self.angular_inertia_tensor().shifted(self.mass, offset)
    }

    /// Computes the world-space inverse angular inertia tensor with the square root of each element.
    #[inline]
    pub fn global_angular_inertia_tensor(&self, rotation: Quat) -> AngularInertiaTensor {
        let mut lhs = Mat3::from_quat(rotation * self.local_inertial_frame);
        let rhs = lhs.transpose();

        lhs.x_axis *= self.principal_angular_inertia.x;
        lhs.y_axis *= self.principal_angular_inertia.y;
        lhs.z_axis *= self.principal_angular_inertia.z;

        AngularInertiaTensor::from_mat3(lhs * rhs)
    }

    /// Returns the mass properties transformed by the given [isometry].
    ///
    /// [isometry]: Isometry3d
    #[inline]
    pub fn transformed_by(mut self, isometry: Isometry3d) -> Self {
        self.transform_by(isometry);
        self
    }

    /// Transforms the mass properties by the given [isometry].
    ///
    /// [isometry]: Isometry3d
    #[inline]
    pub fn transform_by(&mut self, isometry: Isometry3d) {
        self.center_of_mass = self.global_center_of_mass(isometry);
        self.local_inertial_frame = isometry.rotation * self.local_inertial_frame;
    }

    /// Sets the mass to the given `new_mass`. This also affects the angular inertia.
    #[inline]
    pub fn set_mass(&mut self, new_mass: impl Into<Mass>) {
        let new_mass = new_mass.into();

        // Adjust angular inertia based on new mass.
        self.principal_angular_inertia /= self.mass.value();
        self.principal_angular_inertia *= new_mass.value();

        self.mass = new_mass;
    }
}

impl std::ops::Add for MassProperties3d {
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
        let i1 = self.shifted_angular_inertia_tensor(new_center_of_mass - self.center_of_mass);
        let i2 = self.shifted_angular_inertia_tensor(new_center_of_mass - other.center_of_mass);
        let new_angular_inertia = AngularInertiaTensor::from_mat3(i1.value() + i2.value());

        Self::new_with_angular_inertia_tensor(new_mass, new_angular_inertia, new_center_of_mass)
    }
}

impl std::ops::AddAssign for MassProperties3d {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl std::ops::Sub for MassProperties3d {
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
        let i1 = self.shifted_angular_inertia_tensor(new_center_of_mass - self.center_of_mass);
        let i2 = self.shifted_angular_inertia_tensor(new_center_of_mass - other.center_of_mass);
        let new_angular_inertia = AngularInertiaTensor::from_mat3(i1.value() - i2.value());

        Self::new_with_angular_inertia_tensor(new_mass, new_angular_inertia, new_center_of_mass)
    }
}

impl std::ops::SubAssign for MassProperties3d {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}
