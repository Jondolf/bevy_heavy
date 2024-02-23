use crate::recip_or_zero;
use bevy_math::{Mat3, Quat, Vec3};

/// [`ComputeMassProperties3d`] implementations for 3D geometric primitives.
mod impls;

mod eigen3;
pub use eigen3::SymmetricEigen3;

/// A trait for computing [`MassProperties3d`] for 3D objects.
pub trait ComputeMassProperties3d {
    /// Computes the mass of the object with a given `density`.
    fn mass(&self, density: f32) -> f32;

    /// Computes the angular inertia corresponding to the given `mass`.
    #[doc(alias = "moment_of_inertia")]
    fn angular_inertia(&self, mass: f32) -> Vec3;

    /// Computes the angular inertia local frame.
    fn angular_inertia_local_frame(&self) -> Quat {
        Quat::IDENTITY
    }

    /// Computes the local center of mass.
    fn center_of_mass(&self) -> Vec3;

    /// Computes the [`MassProperties3d`] with a given `density`.
    fn mass_properties(&self, density: f32) -> MassProperties3d {
        let mass = self.mass(density);
        MassProperties3d::new_with_local_frame(
            mass,
            self.angular_inertia(mass),
            self.angular_inertia_local_frame(),
            self.center_of_mass(),
        )
    }
}

/// Mass properties in 3D.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassProperties3d {
    /// The multiplicative inverse of mass, `1.0 / mass`.
    pub inverse_mass: f32,
    /// The square root of the multiplicative inverse of the angular inertia
    /// along the principal axes.
    pub inverse_angular_inertia_sqrt: Vec3,
    /// The principal angular inertia local frame or orientation as a `Quat`.
    pub angular_inertia_local_frame: Quat,
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
        inverse_angular_inertia_sqrt: Vec3::ZERO,
        angular_inertia_local_frame: Quat::IDENTITY,
        center_of_mass: Vec3::ZERO,
    };

    /// Creates a new [`MassProperties3d`] from a given mass, principal angular inertia,
    /// and center of mass in local space.
    pub fn new(mass: f32, angular_inertia: Vec3, center_of_mass: Vec3) -> Self {
        Self::new_with_local_frame(mass, angular_inertia, Quat::IDENTITY, center_of_mass)
    }

    /// Creates a new [`MassProperties3d`] from a given mass, principal angular inertia,
    /// angular inertia local frame, and center of mass in local space.
    ///
    /// The principal angular inertia is the angular inertia along the coordinate axes defined
    /// by `angular_inertia_local_frame`, expressed in local space.
    pub fn new_with_local_frame(
        mass: f32,
        angular_inertia: Vec3,
        angular_inertia_local_frame: Quat,
        center_of_mass: Vec3,
    ) -> Self {
        Self {
            inverse_mass: recip_or_zero(mass),
            inverse_angular_inertia_sqrt: angular_inertia
                .to_array()
                .map(|val| recip_or_zero(val.sqrt()))
                .into(),
            angular_inertia_local_frame,
            center_of_mass,
        }
    }

    /// Creates a new [`MassProperties3d`] from a given mass, angular inertia matrix,
    /// and center of mass in local space.
    ///
    /// The angular inertia matrix will be diagonalized in order to extract the principal inertia
    /// values and local principal inertia frme.
    pub fn new_with_inertia_matrix(mass: f32, angular_inertia: Mat3, center_of_mass: Vec3) -> Self {
        let mut eigen = SymmetricEigen3::new(angular_inertia).reverse();

        if eigen.eigenvectors.determinant() < 0.0 {
            std::mem::swap(
                &mut eigen.eigenvectors.y_axis,
                &mut eigen.eigenvectors.z_axis,
            );
            std::mem::swap(&mut eigen.eigenvalues.y, &mut eigen.eigenvalues.z);
        }

        let mut principal_inertia_local_frame = Quat::from_mat3(&eigen.eigenvectors).normalize();

        if !principal_inertia_local_frame.is_finite() {
            principal_inertia_local_frame = Quat::IDENTITY;
        }

        // Clamp eigenvalues to be non-negative.
        let principal_inertia = eigen.eigenvalues.max(Vec3::ZERO);

        Self::new_with_local_frame(
            mass,
            principal_inertia,
            principal_inertia_local_frame,
            center_of_mass,
        )
    }

    /// Returns the mass.
    pub fn mass(&self) -> f32 {
        recip_or_zero(self.inverse_mass)
    }

    /// Returns the principal angular inerta.
    pub fn angular_inertia(&self) -> Vec3 {
        self.inverse_angular_inertia_sqrt
            .to_array()
            .map(|val| recip_or_zero(val.powi(2)))
            .into()
    }

    /// Returns the center of mass transformed into global space using `translation` and `rotation`.
    pub fn global_center_of_mass(&self, translation: Vec3, rotation: Quat) -> Vec3 {
        translation + rotation * self.center_of_mass
    }

    /// Computes the world-space inverse angular inertia tensor with the square root of each element.
    pub fn global_inverse_inertia_tensor_sqrt(&self, rotation: Quat) -> Mat3 {
        if self.inverse_angular_inertia_sqrt != Vec3::ZERO {
            let mut lhs = Mat3::from_quat(rotation * self.angular_inertia_local_frame);
            let rhs = lhs.transpose();

            lhs.x_axis *= self.inverse_angular_inertia_sqrt.x;
            lhs.y_axis *= self.inverse_angular_inertia_sqrt.y;
            lhs.z_axis *= self.inverse_angular_inertia_sqrt.z;

            lhs * rhs
        } else {
            Mat3::ZERO
        }
    }

    /// Computes the world-space angular inertia tensor from the principal inertia.
    pub fn angular_inertia_tensor(&self) -> Mat3 {
        if self.inverse_angular_inertia_sqrt != Vec3::ZERO {
            let inertia = self.inverse_angular_inertia_sqrt.powf(2.0).recip();
            Mat3::from_quat(self.angular_inertia_local_frame)
                * Mat3::from_diagonal(inertia)
                * Mat3::from_quat(self.angular_inertia_local_frame.inverse())
        } else {
            Mat3::ZERO
        }
    }

    /// Computes the world-space inverse angular inertia tensor from the principal inertia.
    pub fn inverse_angular_inertia_tensor(&self) -> Mat3 {
        let inverse_inertia = self.inverse_angular_inertia_sqrt.powf(2.0);
        Mat3::from_quat(self.angular_inertia_local_frame)
            * Mat3::from_diagonal(inverse_inertia)
            * Mat3::from_quat(self.angular_inertia_local_frame.inverse())
    }

    /// Computes the inertia tensor at a given `offset`.
    pub fn shifted_angular_inertia(&self, offset: Vec3) -> Mat3 {
        let inertia_tensor = self.angular_inertia_tensor();

        if self.inverse_mass > f32::EPSILON {
            let mass = 1.0 / self.inverse_mass;
            let diagonal_element = offset.length_squared();
            let diagonal_mat = Mat3::from_diagonal(Vec3::splat(diagonal_element));
            let offset_outer_product =
                Mat3::from_cols(offset * offset.x, offset * offset.y, offset * offset.z);
            inertia_tensor + (diagonal_mat + offset_outer_product) * mass
        } else {
            inertia_tensor
        }
    }

    /// Returns the mass properties transformed by `translation` and `rotation`.
    pub fn transformed_by(mut self, translation: Vec3, rotation: Quat) -> Self {
        self.transform_by(translation, rotation);
        self
    }

    /// Transforms the mass properties by `translation` and `rotation`.
    pub fn transform_by(&mut self, translation: Vec3, rotation: Quat) {
        self.center_of_mass = self.global_center_of_mass(translation, rotation);
        self.angular_inertia_local_frame = rotation * self.angular_inertia_local_frame;
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

impl std::ops::Add for MassProperties3d {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        if self == Self::ZERO {
            return other;
        } else if other == Self::ZERO {
            return self;
        }

        let mass1 = self.mass();
        let mass2 = other.mass();
        let new_mass = mass1 + mass2;
        let new_inverse_mass = recip_or_zero(new_mass);

        // The new center of mass is the weighted average of the centers of masses of `self` and `other`.
        let new_center_of_mass =
            (self.center_of_mass * mass1 + other.center_of_mass * mass2) * new_inverse_mass;

        // Compute the new principal angular inertia, taking the new center of mass into account.
        let inertia1 = self.shifted_angular_inertia(new_center_of_mass - self.center_of_mass);
        let inertia2 = self.shifted_angular_inertia(new_center_of_mass - other.center_of_mass);
        let new_inertia = inertia1 + inertia2;

        Self::new_with_inertia_matrix(new_mass, new_inertia, new_center_of_mass)
    }
}

impl std::ops::AddAssign for MassProperties3d {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl std::ops::Sub for MassProperties3d {
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

        Self::new_with_inertia_matrix(new_mass, new_inertia, new_center_of_mass)
    }
}

impl std::ops::SubAssign for MassProperties3d {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}
