use crate::Mass;

use super::{AngularInertiaTensor, ComputeMassProperties3d, MassProperties3d};
use bevy_math::{
    prelude::Tetrahedron,
    primitives::{
        BoxedPolyline3d, Capsule3d, Cone, ConicalFrustum, Cuboid, Cylinder, Line3d, Measured3d,
        Plane3d, Polyline3d, Segment3d, Sphere, Torus,
    },
    Mat3, Quat, Vec3,
};

impl ComputeMassProperties3d for Sphere {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.volume() * density)
    }

    #[inline]
    fn unit_principal_angular_inertia(&self) -> Vec3 {
        Vec3::splat(0.4 * self.radius.powi(2))
    }

    #[inline]
    fn unit_angular_inertia_tensor(&self) -> AngularInertiaTensor {
        AngularInertiaTensor::new(self.unit_principal_angular_inertia())
    }

    #[inline]
    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }
}

impl ComputeMassProperties3d for Cuboid {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.volume() * density)
    }

    #[inline]
    fn unit_principal_angular_inertia(&self) -> Vec3 {
        let ix = self.half_size.x.powi(2) / 3.0;
        let iy = self.half_size.y.powi(2) / 3.0;
        let iz = self.half_size.z.powi(2) / 3.0;
        Vec3::new(iy + iz, ix + iz, ix + iy)
    }

    #[inline]
    fn unit_angular_inertia_tensor(&self) -> AngularInertiaTensor {
        AngularInertiaTensor::new(self.unit_principal_angular_inertia())
    }

    #[inline]
    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }
}

impl ComputeMassProperties3d for Cylinder {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.volume() * density)
    }

    #[inline]
    fn unit_principal_angular_inertia(&self) -> Vec3 {
        let radius_squared = self.radius.powi(2);
        let height_squared = self.half_height.powi(2) * 4.0;
        let principal = radius_squared / 2.0;
        let off_principal = (radius_squared * 3.0 + height_squared) / 12.0;
        Vec3::new(off_principal, principal, off_principal)
    }

    #[inline]
    fn unit_angular_inertia_tensor(&self) -> AngularInertiaTensor {
        AngularInertiaTensor::new(self.unit_principal_angular_inertia())
    }

    #[inline]
    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }
}

impl ComputeMassProperties3d for Capsule3d {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.volume() * density)
    }

    #[inline]
    fn unit_principal_angular_inertia(&self) -> Vec3 {
        // The cylinder and hemisphere parts
        let cylinder = Cylinder {
            radius: self.radius,
            half_height: self.half_length,
        };
        let cylinder_length = self.half_length * 2.0;
        let sphere = Sphere::new(self.radius);

        // Volumes
        let cylinder_volume = cylinder.volume();
        let sphere_volume = sphere.volume();

        // Masses
        let density = 1.0 / (cylinder_volume + sphere_volume);
        let cylinder_mass = cylinder_volume * density;
        let sphere_mass = sphere_volume * density;

        // Principal inertias
        let cylinder_inertia = cylinder.principal_angular_inertia(cylinder_mass);
        let sphere_inertia = sphere.principal_angular_inertia(sphere_mass);

        // Total inertia
        let mut capsule_inertia = cylinder_inertia + sphere_inertia;

        // Compensate for the hemispheres being away from the rotation axis using the parallel axis theorem.
        let extra = (cylinder_length.powi(2) * 0.25 + cylinder_length * self.radius * 3.0 / 8.0)
            * sphere_mass;
        capsule_inertia.x += extra;
        capsule_inertia.z += extra;

        capsule_inertia
    }

    #[inline]
    fn unit_angular_inertia_tensor(&self) -> AngularInertiaTensor {
        AngularInertiaTensor::new(self.unit_principal_angular_inertia())
    }

    #[inline]
    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }

    #[inline]
    fn mass_properties(&self, density: f32) -> MassProperties3d {
        // The cylinder and hemisphere parts
        let cylinder = Cylinder {
            radius: self.radius,
            half_height: self.half_length,
        };
        let cylinder_length = self.half_length * 2.0;
        let sphere = Sphere::new(self.radius);

        // Volumes
        let cylinder_volume = cylinder.volume();
        let sphere_volume = sphere.volume();

        // Masses
        let cylinder_mass = cylinder_volume * density;
        let sphere_mass = sphere_volume * density;

        // Principal inertias
        let cylinder_inertia = cylinder.principal_angular_inertia(cylinder_mass);
        let sphere_inertia = sphere.principal_angular_inertia(sphere_mass);

        // Total inertia
        let mut capsule_inertia = cylinder_inertia + sphere_inertia;

        // Compensate for the hemispheres being away from the rotation axis using the parallel axis theorem.
        let extra = (cylinder_length.powi(2) * 0.25 + cylinder_length * self.radius * 3.0 / 8.0)
            * sphere_mass;
        capsule_inertia.x += extra;
        capsule_inertia.z += extra;

        MassProperties3d::new(cylinder_mass + sphere_mass, capsule_inertia, Vec3::ZERO)
    }
}

impl ComputeMassProperties3d for Cone {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.volume() * density)
    }

    #[inline]
    fn unit_principal_angular_inertia(&self) -> Vec3 {
        let radius_squared = self.radius.powi(2);
        let height_squared = self.height.powi(2);
        let principal = radius_squared * 3.0 / 10.0;
        let off_principal = radius_squared * 3.0 / 20.0 + height_squared * 3.0 / 80.0;
        Vec3::new(off_principal, principal, off_principal)
    }

    #[inline]
    fn unit_angular_inertia_tensor(&self) -> AngularInertiaTensor {
        AngularInertiaTensor::new(self.unit_principal_angular_inertia())
    }

    #[inline]
    fn center_of_mass(&self) -> Vec3 {
        Vec3::new(0.0, -self.height * 0.25, 0.0)
    }
}

impl ComputeMassProperties3d for ConicalFrustum {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        if self.radius_top == self.radius_bottom {
            Cylinder::new(self.radius_top, self.height).mass(density)
        } else {
            let (min_radius, max_radius) = if self.radius_top < self.radius_bottom {
                (self.radius_top, self.radius_bottom)
            } else {
                (self.radius_bottom, self.radius_top)
            };

            // To compute the mass of the conical frustum, we will subtract the volume
            // of the smaller cone (2) from the mass properties of the larger cone (1).

            // The height of the cone formed by tapering the conical frustum to a tip.
            let cone_height =
                max_radius * (self.height / (self.radius_top - self.radius_bottom).abs());

            // Compute the frustum volume as the difference of the volumes of the larger and smaller cone.
            let volume1 = std::f32::consts::PI / 3.0 * max_radius.powi(2) * cone_height;
            let volume2 =
                std::f32::consts::PI / 3.0 * min_radius.powi(2) * (cone_height - self.height);
            let frustum_volume = volume1 - volume2;

            Mass::new(frustum_volume * density)
        }
    }

    #[inline]
    fn unit_principal_angular_inertia(&self) -> Vec3 {
        if self.radius_top == self.radius_bottom {
            Cylinder::new(self.radius_top, self.height).unit_principal_angular_inertia()
        } else {
            let max_radius = self.radius_top.max(self.radius_bottom);
            let radii_squared = self.radius_top.powi(2) + self.radius_bottom.powi(2);

            // The height of the cone formed by tapering the conical frustum to a tip.
            let cone_height =
                max_radius * (self.height / (self.radius_top - self.radius_bottom).abs());

            // Compute the angular inertia.
            // FIXME: `off_principal` is the same as it is for the cone. Is this correct?
            let principal = 1.0 / 12.0 * (self.height.powi(2) + 3.0 * radii_squared);
            let off_principal = max_radius.powi(2) * 3.0 / 20.0 + cone_height.powi(2) * 3.0 / 80.0;

            Vec3::new(off_principal, principal, off_principal)
        }
    }

    #[inline]
    fn unit_angular_inertia_tensor(&self) -> AngularInertiaTensor {
        AngularInertiaTensor::new(self.unit_principal_angular_inertia())
    }

    #[inline]
    fn center_of_mass(&self) -> Vec3 {
        if self.radius_top == self.radius_bottom {
            Vec3::ZERO
        } else {
            let (min_radius, max_radius) = if self.radius_top < self.radius_bottom {
                (self.radius_top, self.radius_bottom)
            } else {
                (self.radius_bottom, self.radius_top)
            };

            // To compute the center of mass of the conical frustum, we will subtract the mass properties
            // of the smaller cone (2) from the mass properties of the larger cone (1).

            let cone1_height =
                max_radius * (self.height / (self.radius_top - self.radius_bottom).abs());
            let cone2_height = cone1_height - self.height;

            // Compute the frustum volume as the difference of the volumes of the cones.
            let volume1 = std::f32::consts::PI / 3.0 * max_radius.powi(2) * cone1_height;
            let volume2 =
                std::f32::consts::PI / 3.0 * min_radius.powi(2) * (cone1_height - self.height);

            // The Y coordinates of the centers of masses
            let (com1, com2) = if self.radius_top < self.radius_bottom {
                (-cone1_height / 4.0, self.height / 2.0 + cone2_height / 4.0)
            } else {
                (
                    cone1_height / 4.0,
                    (cone1_height - self.height) / 4.0 - self.height / 2.0,
                )
            };

            // We can use volume instead of mass for the weighted average because the density is constant.
            Vec3::new(
                0.0,
                (com1 * volume1 - com2 * volume2) / (volume1 - volume2),
                0.0,
            )
        }
    }

    #[inline]
    fn mass_properties(&self, density: f32) -> MassProperties3d {
        if self.radius_top == self.radius_bottom {
            Cylinder::new(self.radius_top, self.height).mass_properties(density)
        } else {
            let (min_radius, max_radius) = if self.radius_top < self.radius_bottom {
                (self.radius_top, self.radius_bottom)
            } else {
                (self.radius_bottom, self.radius_top)
            };
            let radii_squared = self.radius_top.powi(2) + self.radius_bottom.powi(2);

            // To compute the mass properties of the conical frustum, we will subtract the mass properties
            // of the smaller cone (2) from the mass properties of the larger cone (1).

            let cone1_height =
                max_radius * (self.height / (self.radius_top - self.radius_bottom).abs());
            let cone2_height = cone1_height - self.height;

            // Compute the frustum volume as the difference of the volumes of the cones.
            let volume1 = std::f32::consts::PI / 3.0 * max_radius.powi(2) * cone1_height;
            let volume2 =
                std::f32::consts::PI / 3.0 * min_radius.powi(2) * (cone1_height - self.height);
            let frustum_volume = volume1 - volume2;
            let frustum_mass = frustum_volume * density;

            // The Y coordinates of the centers of masses
            let (com1, com2) = if self.radius_top < self.radius_bottom {
                (-cone1_height / 4.0, self.height / 2.0 + cone2_height / 4.0)
            } else {
                (
                    cone1_height / 4.0,
                    (cone1_height - self.height) / 4.0 - self.height / 2.0,
                )
            };

            // We can use volume instead of mass for the weighted average because the density is constant.
            let center_of_mass = Vec3::new(
                0.0,
                (com1 * volume1 - com2 * volume2) / (volume1 - volume2),
                0.0,
            );

            // Compute the angular inertia.
            // FIXME: `off_principal` is the same as it is for the cone. Is this correct?
            let principal = 1.0 / 12.0 * (self.height.powi(2) + 3.0 * radii_squared);
            let off_principal = max_radius.powi(2) * 3.0 / 20.0 + cone1_height.powi(2) * 3.0 / 80.0;
            let inertia = Vec3::new(off_principal, principal, off_principal) * frustum_mass;

            MassProperties3d::new(frustum_mass, inertia, center_of_mass)
        }
    }
}

impl ComputeMassProperties3d for Torus {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.volume() * density)
    }

    #[inline]
    fn unit_principal_angular_inertia(&self) -> Vec3 {
        // Reference: https://en.wikipedia.org/wiki/List_of_moments_of_inertia

        let major_radius_squared = self.major_radius.powi(2);
        let minor_radius_squared = self.minor_radius.powi(2);

        let principal = 0.25 * (4.0 * major_radius_squared + 3.0 * minor_radius_squared);
        let off_principal = 1.0 / 8.0 * (4.0 * major_radius_squared + 5.0 * minor_radius_squared);
        Vec3::new(off_principal, principal, off_principal)
    }

    #[inline]
    fn unit_angular_inertia_tensor(&self) -> AngularInertiaTensor {
        AngularInertiaTensor::new(self.unit_principal_angular_inertia())
    }

    #[inline]
    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }
}

impl ComputeMassProperties3d for Tetrahedron {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.volume() * density)
    }

    #[inline]
    fn unit_principal_angular_inertia(&self) -> Vec3 {
        let tensor = self.angular_inertia_tensor(Mass::ONE);
        tensor.principal_angular_inertia()
    }

    #[inline]
    fn local_inertial_frame(&self) -> Quat {
        let tensor = self.angular_inertia_tensor(Mass::ONE);
        tensor.local_frame()
    }

    #[inline]
    fn unit_angular_inertia_tensor(&self) -> AngularInertiaTensor {
        // References:
        // - F. Tonon. "Explicit Exact Formulas for the 3-D Tetrahedron Inertia Tensor in Terms of its Vertex Coordinates"
        // - Parry: https://github.com/dimforge/parry/blob/837291f1a051dc04e6e4ab8bc2e5b438ce66d257/src/mass_properties/mass_properties_trimesh3d.rs#L38

        let [p1, p2, p3, p4] = self.vertices;

        // For readability
        let x1 = p1.x;
        let y1 = p1.y;
        let z1 = p1.z;
        let x2 = p2.x;
        let y2 = p2.y;
        let z2 = p2.z;
        let x3 = p3.x;
        let y3 = p3.y;
        let z3 = p3.z;
        let x4 = p4.x;
        let y4 = p4.y;
        let z4 = p4.z;

        let diag_x = x1 * x1
            + x1 * x2
            + x2 * x2
            + x1 * x3
            + x2 * x3
            + x3 * x3
            + x1 * x4
            + x2 * x4
            + x3 * x4
            + x4 * x4;
        let diag_y = y1 * y1
            + y1 * y2
            + y2 * y2
            + y1 * y3
            + y2 * y3
            + y3 * y3
            + y1 * y4
            + y2 * y4
            + y3 * y4
            + y4 * y4;
        let diag_z = z1 * z1
            + z1 * z2
            + z2 * z2
            + z1 * z3
            + z2 * z3
            + z3 * z3
            + z1 * z4
            + z2 * z4
            + z3 * z4
            + z4 * z4;

        let a0 = (diag_y + diag_z) * 0.1;
        let b0 = (diag_z + diag_x) * 0.1;
        let c0 = (diag_x + diag_y) * 0.1;

        let a1 = (y1 * z1 * 2.0
            + y2 * z1
            + y3 * z1
            + y4 * z1
            + y1 * z2
            + y2 * z2 * 2.0
            + y3 * z2
            + y4 * z2
            + y1 * z3
            + y2 * z3
            + y3 * z3 * 2.0
            + y4 * z3
            + y1 * z4
            + y2 * z4
            + y3 * z4
            + y4 * z4 * 2.0)
            * 0.05;
        let b1 = (x1 * z1 * 2.0
            + x2 * z1
            + x3 * z1
            + x4 * z1
            + x1 * z2
            + x2 * z2 * 2.0
            + x3 * z2
            + x4 * z2
            + x1 * z3
            + x2 * z3
            + x3 * z3 * 2.0
            + x4 * z3
            + x1 * z4
            + x2 * z4
            + x3 * z4
            + x4 * z4 * 2.0)
            * 0.05;
        let c1 = (x1 * y1 * 2.0
            + x2 * y1
            + x3 * y1
            + x4 * y1
            + x1 * y2
            + x2 * y2 * 2.0
            + x3 * y2
            + x4 * y2
            + x1 * y3
            + x2 * y3
            + x3 * y3 * 2.0
            + x4 * y3
            + x1 * y4
            + x2 * y4
            + x3 * y4
            + x4 * y4 * 2.0)
            * 0.05;

        AngularInertiaTensor::from_mat3(Mat3::from_cols_array(&[
            a0, -b1, -c1, -b1, b0, -a1, -c1, -a1, c0,
        ]))
    }

    #[inline]
    fn center_of_mass(&self) -> Vec3 {
        (self.vertices[0] + self.vertices[1] + self.vertices[2] + self.vertices[3]) / 4.0
    }

    #[inline]
    fn mass_properties(&self, density: f32) -> MassProperties3d {
        let volume = self.volume();
        let center_of_mass = self.center_of_mass();

        if volume < f32::EPSILON {
            return MassProperties3d::new(Mass::ZERO, Vec3::ZERO, center_of_mass);
        }

        let mass = Mass::new(volume * density);
        let tensor = self.angular_inertia_tensor(mass);

        MassProperties3d::new_with_angular_inertia_tensor(mass, tensor, center_of_mass)
    }
}

macro_rules! impl_zero_mass_properties_3d {
    ($($shape:ty),*) => {
        $(
            impl ComputeMassProperties3d for $shape {
                #[inline]
                fn mass(&self, _density: f32) -> Mass {
                    Mass::ZERO
                }

                #[inline]
                fn unit_principal_angular_inertia(&self) -> Vec3 {
                    Vec3::ZERO
                }

                #[inline]
                fn principal_angular_inertia(&self, _mass: impl Into<Mass>) -> Vec3 {
                    Vec3::ZERO
                }

                #[inline]
                fn local_inertial_frame(&self) -> Quat {
                    Quat::IDENTITY
                }

                #[inline]
                fn unit_angular_inertia_tensor(&self) -> AngularInertiaTensor {
                    AngularInertiaTensor::ZERO
                }

                #[inline]
                fn angular_inertia_tensor(&self, _mass: impl Into<Mass>) -> AngularInertiaTensor {
                    AngularInertiaTensor::ZERO
                }

                #[inline]
                fn center_of_mass(&self) -> Vec3 {
                    Vec3::ZERO
                }

                #[inline]
                fn mass_properties(&self, _density: f32) -> MassProperties3d {
                    MassProperties3d::ZERO
                }
            }
        )*
    };
}

impl_zero_mass_properties_3d!(Plane3d);
impl_zero_mass_properties_3d!(Line3d);
impl_zero_mass_properties_3d!(Segment3d);
impl_zero_mass_properties_3d!(BoxedPolyline3d);

impl<const N: usize> ComputeMassProperties3d for Polyline3d<N> {
    #[inline]
    fn mass(&self, _density: f32) -> Mass {
        Mass::ZERO
    }

    #[inline]
    fn unit_principal_angular_inertia(&self) -> Vec3 {
        Vec3::ZERO
    }

    #[inline]
    fn principal_angular_inertia(&self, _mass: impl Into<Mass>) -> Vec3 {
        Vec3::ZERO
    }

    #[inline]
    fn local_inertial_frame(&self) -> bevy_math::Quat {
        Quat::IDENTITY
    }

    #[inline]
    fn unit_angular_inertia_tensor(&self) -> AngularInertiaTensor {
        AngularInertiaTensor::ZERO
    }

    #[inline]
    fn angular_inertia_tensor(&self, _mass: impl Into<Mass>) -> AngularInertiaTensor {
        AngularInertiaTensor::ZERO
    }

    #[inline]
    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }

    #[inline]
    fn mass_properties(&self, _density: f32) -> MassProperties3d {
        MassProperties3d::ZERO
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    macro_rules! test_shape {
        ($test_name:tt, $shape:expr) => {
            #[test]
            fn $test_name() {
                let shape = $shape;

                /*
                // TODO: Add this back when uniform sampling is implemented for all the shapes.
                // Sample enough points to have a close enough point cloud representation of the shape.
                let mut rng = rand_chacha::ChaCha8Rng::from_seed(Default::default());
                let points = (0..2_000_000)
                    .map(|_| shape.sample_interior(&mut rng))
                    .collect::<Vec<_>>();
                */

                // Compute the mass properties to test.
                let density = 2.0;
                let mass = shape.mass(density);
                let principal_angular_inertia = shape.principal_angular_inertia(mass);
                let local_inertial_frame = shape.local_inertial_frame();
                let angular_inertia_tensor = shape.angular_inertia_tensor(mass);
                let center_of_mass = shape.center_of_mass();

                // First, test that the individually computed properties match the full properties.
                let mass_props = shape.mass_properties(density);
                assert_relative_eq!(mass.value(), mass_props.mass.value());
                assert_relative_eq!(
                    principal_angular_inertia,
                    mass_props.principal_angular_inertia
                );
                assert_relative_eq!(local_inertial_frame, mass_props.local_inertial_frame);
                assert_relative_eq!(center_of_mass, mass_props.center_of_mass);

                /*
                // Estimate the expected mass properties using the point cloud.
                // Note: We could also approximate the mass using Monte Carlo integration.
                //       This would require point containment checks.
                let expected =
                MassProperties3d::from_point_cloud(&points, mass, local_inertial_frame);

                assert_relative_eq!(mass.value(), expected.mass.value());
                assert_relative_eq!(
                    principal_angular_inertia,
                    expected.principal_angular_inertia,
                    epsilon = 0.1
                );
                assert_relative_eq!(center_of_mass, expected.center_of_mass, epsilon = 0.01);
                */

                // Check that the computed tensor matches the tensor computed from
                // the principal moments of inertia and the local inertial frame.
                assert_relative_eq!(
                    angular_inertia_tensor,
                    AngularInertiaTensor::new_with_local_frame(
                        principal_angular_inertia,
                        local_inertial_frame
                    ),
                    epsilon = 1e-5
                );

                // Check that the principal moments of inertia and the local inertial frame
                // extracted from the tensor produce the original tensor. This tests that
                // the diagonalization is correct.
                // Note: The principal moments of inertia are *not* unique, and can differ
                //       from the original values while still producing the same tensor
                //       when taking the local inertial frame into account.
                let (principal, frame) =
                    angular_inertia_tensor.principal_angular_inertia_with_local_frame();
                assert_relative_eq!(
                    (frame * principal).abs(),
                    principal_angular_inertia,
                    epsilon = 1e-5
                );
                assert_relative_eq!(
                    AngularInertiaTensor::new_with_local_frame(principal, frame),
                    angular_inertia_tensor,
                    epsilon = 1e-5
                );
            }
        };
    }

    test_shape!(sphere, Sphere::new(2.0));
    test_shape!(cuboid, Cuboid::new(1.0, 2.0, 3.0));
    test_shape!(cylinder, Cylinder::new(1.0, 4.0));
    test_shape!(capsule, Capsule3d::new(1.0, 2.0));
    test_shape!(
        cone,
        Cone {
            radius: 1.0,
            height: 2.0,
        }
    );
    test_shape!(
        conical_frustum,
        ConicalFrustum {
            radius_top: 1.0,
            radius_bottom: 2.0,
            height: 2.0
        }
    );
    test_shape!(torus, Torus::new(1.0, 2.0));
    test_shape!(tetrahedron, Tetrahedron::default());
}
