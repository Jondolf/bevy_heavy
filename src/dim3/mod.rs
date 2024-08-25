use bevy_math::{Mat3, Quat, Vec3};

mod angular_inertia;
pub use angular_inertia::AngularInertiaTensor;

/// [`ComputeMassProperties3d`] implementations for 3D geometric primitives.
mod impls;

mod eigen3;
pub use eigen3::SymmetricEigen3;

use crate::{Mass, RecipOrZero};

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

/// Mass properties in 3D.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "bevy_reflect", derive(bevy_reflect::Reflect))]
#[cfg_attr(feature = "serialize", derive(serde::Serialize, serde::Deserialize))]
pub struct MassProperties3d {
    /// The multiplicative inverse of mass, `1.0 / mass`.
    pub inverse_mass: f32,
    /// The square root of the multiplicative inverse of the angular inertia
    /// along the principal axes.
    pub inverse_principal_angular_inertia_sqrt: Vec3,
    /// The principal angular inertia local frame or orientation as a `Quat`.
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
        inverse_mass: 0.0,
        inverse_principal_angular_inertia_sqrt: Vec3::ZERO,
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
    /// by `angular_inertia_local_frame`, expressed in local space.
    #[inline]
    pub fn new_with_local_frame(
        mass: impl Into<Mass>,
        principal_angular_inertia: Vec3,
        local_inertial_frame: Quat,
        center_of_mass: Vec3,
    ) -> Self {
        Self {
            inverse_mass: mass.into().recip_or_zero(),
            inverse_principal_angular_inertia_sqrt: principal_angular_inertia
                .to_array()
                .map(|val| val.sqrt().recip_or_zero())
                .into(),
            local_inertial_frame,
            center_of_mass,
        }
    }

    /// Creates a new [`MassProperties3d`] from a given mass, angular inertia matrix,
    /// and center of mass in local space.
    ///
    /// The angular inertia matrix will be diagonalized in order to extract the principal inertia
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

    /// Returns the mass.
    #[inline]
    pub fn mass(&self) -> f32 {
        self.inverse_mass.recip_or_zero()
    }

    /// Returns the principal angular inerta.
    #[inline]
    pub fn principal_angular_inertia(&self) -> Vec3 {
        self.inverse_principal_angular_inertia_sqrt
            .to_array()
            .map(|val| val.powi(2).recip_or_zero())
            .into()
    }

    /// Returns the center of mass transformed into global space using `translation` and `rotation`.
    #[inline]
    pub fn global_center_of_mass(&self, translation: Vec3, rotation: Quat) -> Vec3 {
        translation + rotation * self.center_of_mass
    }

    /// Computes the world-space inverse angular inertia tensor with the square root of each element.
    #[inline]
    pub fn global_inverse_angular_inertia_tensor_sqrt(&self, rotation: Quat) -> Mat3 {
        if self.inverse_principal_angular_inertia_sqrt != Vec3::ZERO {
            let mut lhs = Mat3::from_quat(rotation * self.local_inertial_frame);
            let rhs = lhs.transpose();

            lhs.x_axis *= self.inverse_principal_angular_inertia_sqrt.x;
            lhs.y_axis *= self.inverse_principal_angular_inertia_sqrt.y;
            lhs.z_axis *= self.inverse_principal_angular_inertia_sqrt.z;

            lhs * rhs
        } else {
            Mat3::ZERO
        }
    }

    /// Computes the world-space angular inertia tensor from the principal inertia.
    #[inline]
    pub fn angular_inertia_tensor(&self) -> AngularInertiaTensor {
        if self.inverse_principal_angular_inertia_sqrt != Vec3::ZERO {
            let inertia = self
                .inverse_principal_angular_inertia_sqrt
                .powf(2.0)
                .recip();
            AngularInertiaTensor::new(inertia)
        } else {
            AngularInertiaTensor::ZERO
        }
    }

    /// Computes the world-space inverse angular inertia tensor from the principal inertia.
    #[inline]
    pub fn inverse_angular_inertia_tensor(&self) -> AngularInertiaTensor {
        let inverse_inertia = self.inverse_principal_angular_inertia_sqrt.powf(2.0);
        AngularInertiaTensor::new(inverse_inertia)
    }

    /// Computes the inertia tensor at a given `offset`.
    pub fn shifted_angular_inertia_tensor(&self, offset: Vec3) -> AngularInertiaTensor {
        let inertia_tensor = self.angular_inertia_tensor();

        if self.inverse_mass > f32::EPSILON {
            let mass = 1.0 / self.inverse_mass;
            let diagonal_element = offset.length_squared();
            let diagonal_mat = Mat3::from_diagonal(Vec3::splat(diagonal_element));
            let offset_outer_product =
                Mat3::from_cols(offset * offset.x, offset * offset.y, offset * offset.z);
            AngularInertiaTensor::from_mat3(
                inertia_tensor.value() + (diagonal_mat + offset_outer_product) * mass,
            )
        } else {
            inertia_tensor
        }
    }

    /// Returns the mass properties transformed by `translation` and `rotation`.
    #[inline]
    pub fn transformed_by(mut self, translation: Vec3, rotation: Quat) -> Self {
        self.transform_by(translation, rotation);
        self
    }

    /// Transforms the mass properties by `translation` and `rotation`.
    #[inline]
    pub fn transform_by(&mut self, translation: Vec3, rotation: Quat) {
        self.center_of_mass = self.global_center_of_mass(translation, rotation);
        self.local_inertial_frame = rotation * self.local_inertial_frame;
    }

    /// Sets the mass to the given `new_mass`. This also affects the angular inertia.
    #[inline]
    pub fn set_mass(&mut self, new_mass: impl Into<Mass>) {
        let new_inverse_mass = new_mass.into().recip_or_zero();

        // Adjust angular inertia based on new mass.
        let old_mass = self.inverse_mass.recip_or_zero();
        self.inverse_principal_angular_inertia_sqrt *= new_inverse_mass.sqrt() * old_mass.sqrt();

        self.inverse_mass = new_inverse_mass;
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

        let mass1 = self.mass();
        let mass2 = other.mass();
        let new_mass = mass1 + mass2;
        let new_inverse_mass = new_mass.recip_or_zero();

        // The new center of mass is the weighted average of the centers of masses of `self` and `other`.
        let new_center_of_mass =
            (self.center_of_mass * mass1 + other.center_of_mass * mass2) * new_inverse_mass;

        // Compute the new principal angular inertia, taking the new center of mass into account.
        let inertia1 = self
            .shifted_angular_inertia_tensor(new_center_of_mass - self.center_of_mass)
            .value();
        let inertia2 = self
            .shifted_angular_inertia_tensor(new_center_of_mass - other.center_of_mass)
            .value();
        let new_inertia = AngularInertiaTensor::from_mat3(inertia1 + inertia2);

        Self::new_with_angular_inertia_tensor(new_mass, new_inertia, new_center_of_mass)
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

        let mass1 = self.mass();
        let mass2 = other.mass();
        let new_mass = mass1 - mass2;

        let new_inverse_mass = if new_mass >= f32::EPSILON {
            (mass1 + mass2).recip_or_zero()
        } else {
            0.0
        };

        // The new center of mass is the negated weighted average of the centers of masses of `self` and `other`.
        let new_center_of_mass =
            (self.center_of_mass * mass1 - other.center_of_mass * mass2) * new_inverse_mass;

        // Compute the new principal angular inertia, taking the new center of mass into account.
        let inertia1 =
            self.shifted_angular_inertia_tensor(new_center_of_mass - self.center_of_mass);
        let inertia2 =
            self.shifted_angular_inertia_tensor(new_center_of_mass - other.center_of_mass);
        let new_inertia = AngularInertiaTensor::from_mat3(inertia1.value() - inertia2.value());

        Self::new_with_angular_inertia_tensor(new_mass, new_inertia, new_center_of_mass)
    }
}

impl std::ops::SubAssign for MassProperties3d {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}
