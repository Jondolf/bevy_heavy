use bevy_math::{DVec3, Mat3, Quat, Vec3};

mod angular_inertia;
pub use angular_inertia::{AngularInertiaTensor, AngularInertiaTensorError};

/// [`ComputeMassProperties3d`] implementations for 3D geometric primitives.
mod impls;

mod eigen3;
pub use eigen3::SymmetricEigen3;

use crate::RecipOrZero;

/// A trait for computing [`MassProperties3d`] for 3D objects.
pub trait ComputeMassProperties3d {
    /// Computes the mass of the object with a given `density`.
    fn mass(&self, density: f32) -> f32;

    /// Computes the principal angular inertia corresponding to a unit mass.
    #[doc(alias = "unit_principal_moment_of_inertia")]
    fn unit_principal_angular_inertia(&self) -> Vec3;

    /// Computes the principal angular inertia corresponding to the given `mass`.
    #[inline]
    #[doc(alias = "principal_moment_of_inertia")]
    fn principal_angular_inertia(&self, mass: f32) -> Vec3 {
        mass * self.unit_principal_angular_inertia()
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
    fn angular_inertia_tensor(&self, mass: f32) -> AngularInertiaTensor {
        mass * self.unit_angular_inertia_tensor()
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

/// The mass, [angular inertia], and local center of mass of an object in 3D space.
///
/// [angular inertia]: AngularInertiaTensor
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "bevy_reflect", derive(bevy_reflect::Reflect))]
#[cfg_attr(feature = "serialize", derive(serde::Serialize, serde::Deserialize))]
pub struct MassProperties3d {
    /// The mass.
    pub mass: f32,
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
        mass: 0.0,
        principal_angular_inertia: Vec3::ZERO,
        local_inertial_frame: Quat::IDENTITY,
        center_of_mass: Vec3::ZERO,
    };

    /// Creates a new [`MassProperties3d`] from a given mass, principal angular inertia,
    /// and center of mass in local space.
    #[inline]
    pub fn new(mass: f32, principal_angular_inertia: Vec3, center_of_mass: Vec3) -> Self {
        Self::new_with_local_frame(
            mass,
            principal_angular_inertia,
            Quat::IDENTITY,
            center_of_mass,
        )
    }

    /// Creates a new [`MassProperties3d`] from a given mass, principal angular inertia,
    /// local inertial frame, and center of mass in local space.
    ///
    /// The principal angular inertia is the angular inertia along the coordinate axes defined
    /// by the `local_inertial_frame`, expressed in local space.
    #[inline]
    pub fn new_with_local_frame(
        mass: f32,
        principal_angular_inertia: Vec3,
        local_inertial_frame: Quat,
        center_of_mass: Vec3,
    ) -> Self {
        Self {
            mass,
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
        mass: f32,
        tensor: impl Into<AngularInertiaTensor>,
        center_of_mass: Vec3,
    ) -> Self {
        let (principal, local_frame) = tensor.into().principal_angular_inertia_with_local_frame();
        Self::new_with_local_frame(mass, principal, local_frame, center_of_mass)
    }

    /// Computes approximate mass properties from the given set of points representing a shape.
    ///
    /// This can be used to estimate mass properties for arbitrary shapes
    /// by providing a set of sample points from inside the shape.
    ///
    /// The more points there are, and the more uniformly distributed they are,
    /// the more accurate the estimation will be.
    #[inline]
    pub fn from_point_cloud(points: &[Vec3], mass: f32, local_inertial_frame: Quat) -> Self {
        let points_recip = 1.0 / points.len() as f64;

        let center_of_mass =
            (points.iter().fold(DVec3::ZERO, |acc, p| acc + p.as_dvec3()) * points_recip).as_vec3();
        let unit_angular_inertia = points
            .iter()
            .fold(DVec3::ZERO, |acc, p| {
                let p = p.as_dvec3() - center_of_mass.as_dvec3();
                let r_x = p.reject_from_normalized(DVec3::X).length_squared();
                let r_y = p.reject_from_normalized(DVec3::Y).length_squared();
                let r_z = p.reject_from_normalized(DVec3::Z).length_squared();
                acc + DVec3::new(r_x, r_y, r_z) * points_recip
            })
            .as_vec3();

        Self::new_with_local_frame(
            mass,
            mass * unit_angular_inertia,
            local_inertial_frame,
            center_of_mass,
        )
    }

    /// Returns the center of mass transformed into global space using the given translation and rotation.
    #[inline]
    pub fn global_center_of_mass(&self, translation: Vec3, rotation: Quat) -> Vec3 {
        translation + rotation * self.center_of_mass
    }

    /// Computes the principal angular inertia corresponding to a mass of `1.0`.
    ///
    /// If the mass is zero, a zero vector is returned.
    #[inline]
    pub fn unit_principal_angular_inertia(&self) -> Vec3 {
        self.mass.recip_or_zero() * self.principal_angular_inertia
    }

    /// Computes the world-space angular inertia tensor corresponding to a mass of `1.0`.
    ///
    /// If the mass is zero, a zero tensor is returned.
    #[inline]
    pub fn unit_angular_inertia_tensor(&self) -> AngularInertiaTensor {
        self.mass.recip_or_zero() * self.angular_inertia_tensor()
    }

    /// Computes the world-space angular inertia tensor from the principal inertia.
    #[inline]
    pub fn angular_inertia_tensor(&self) -> AngularInertiaTensor {
        AngularInertiaTensor::new_with_local_frame(
            self.principal_angular_inertia,
            self.local_inertial_frame,
        )
    }

    /// Computes the angular inertia tensor at a given `offset`.
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

    /// Returns the mass properties transformed by the given translation and rotation.
    #[inline]
    pub fn transformed_by(mut self, translation: Vec3, rotation: Quat) -> Self {
        self.transform_by(translation, rotation);
        self
    }

    /// Transforms the mass properties by the given translation and rotation.
    #[inline]
    pub fn transform_by(&mut self, translation: Vec3, rotation: Quat) {
        self.center_of_mass = self.global_center_of_mass(translation, rotation);
        self.local_inertial_frame = rotation * self.local_inertial_frame;
    }

    /// Sets the mass to the given `new_mass`. This also affects the angular inertia.
    #[inline]
    pub fn set_mass(&mut self, new_mass: f32) {
        // Adjust angular inertia based on new mass.
        self.principal_angular_inertia /= self.mass;
        self.principal_angular_inertia *= new_mass;
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
        let new_center_of_mass =
            (self.center_of_mass * mass1 + other.center_of_mass * mass2) / new_mass;

        // Compute the new principal angular inertia, taking the new center of mass into account.
        let i1 = self.shifted_angular_inertia_tensor(new_center_of_mass - self.center_of_mass);
        let i2 = other.shifted_angular_inertia_tensor(new_center_of_mass - other.center_of_mass);
        let new_angular_inertia = AngularInertiaTensor::from_mat3(i1.as_mat3() + i2.as_mat3());

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

        if mass1 <= mass2 {
            // The result would have non-positive mass.
            return Self {
                center_of_mass: self.center_of_mass,
                ..Self::ZERO
            };
        }

        let new_mass = mass1 - mass2;

        // The new center of mass is the negated weighted average of the centers of masses of `self` and `other`.
        let new_center_of_mass =
            (self.center_of_mass * mass1 - other.center_of_mass * mass2) / new_mass;

        // Compute the new principal angular inertia, taking the new center of mass into account.
        let i1 = self.shifted_angular_inertia_tensor(new_center_of_mass - self.center_of_mass);
        let i2 = other.shifted_angular_inertia_tensor(new_center_of_mass - other.center_of_mass);
        let new_angular_inertia = AngularInertiaTensor::from_mat3(i1.as_mat3() - i2.as_mat3());

        Self::new_with_angular_inertia_tensor(new_mass, new_angular_inertia, new_center_of_mass)
    }
}

impl std::ops::SubAssign for MassProperties3d {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl std::iter::Sum for MassProperties3d {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut total_mass = 0.0;
        let mut total_angular_inertia = AngularInertiaTensor::ZERO;
        let mut total_center_of_mass = Vec3::ZERO;

        // TODO: Avoid this allocation if possible. This is currently needed because we iterate twice.
        let mut all_properties = Vec::with_capacity(iter.size_hint().1.unwrap_or_default());

        for props in iter {
            total_mass += props.mass;
            total_center_of_mass += props.center_of_mass * props.mass;
            all_properties.push(props);
        }

        if total_mass > 0.0 {
            total_center_of_mass /= total_mass;
        }

        for props in all_properties {
            total_angular_inertia +=
                props.shifted_angular_inertia_tensor(total_center_of_mass - props.center_of_mass);
        }

        Self::new_with_angular_inertia_tensor(
            total_mass,
            total_angular_inertia,
            total_center_of_mass,
        )
    }
}

#[cfg(any(feature = "approx", test))]
impl approx::AbsDiffEq for MassProperties3d {
    type Epsilon = f32;
    fn default_epsilon() -> f32 {
        f32::EPSILON
    }
    fn abs_diff_eq(&self, other: &Self, epsilon: f32) -> bool {
        self.mass.abs_diff_eq(&other.mass, epsilon)
            && self
                .principal_angular_inertia
                .abs_diff_eq(other.principal_angular_inertia, epsilon)
            && self
                .local_inertial_frame
                .abs_diff_eq(other.local_inertial_frame, epsilon)
            && self
                .center_of_mass
                .abs_diff_eq(other.center_of_mass, epsilon)
    }
}

#[cfg(any(feature = "approx", test))]
impl approx::RelativeEq for MassProperties3d {
    fn default_max_relative() -> f32 {
        f32::EPSILON
    }
    fn relative_eq(&self, other: &Self, epsilon: f32, max_relative: f32) -> bool {
        self.mass.relative_eq(&other.mass, epsilon, max_relative)
            && self.principal_angular_inertia.relative_eq(
                &other.principal_angular_inertia,
                epsilon,
                max_relative,
            )
            && self.local_inertial_frame.relative_eq(
                &other.local_inertial_frame,
                epsilon,
                max_relative,
            )
            && self
                .center_of_mass
                .relative_eq(&other.center_of_mass, epsilon, max_relative)
    }
}

#[cfg(any(feature = "approx", test))]
impl approx::UlpsEq for MassProperties3d {
    fn default_max_ulps() -> u32 {
        4
    }
    fn ulps_eq(&self, other: &Self, epsilon: f32, max_ulps: u32) -> bool {
        self.mass.ulps_eq(&other.mass, epsilon, max_ulps)
            && self.principal_angular_inertia.ulps_eq(
                &other.principal_angular_inertia,
                epsilon,
                max_ulps,
            )
            && self
                .local_inertial_frame
                .ulps_eq(&other.local_inertial_frame, epsilon, max_ulps)
            && self
                .center_of_mass
                .ulps_eq(&other.center_of_mass, epsilon, max_ulps)
    }
}
