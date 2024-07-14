use super::{ComputeMassProperties3d, MassProperties3d};
use bevy_math::{
    primitives::{
        BoxedPolyline3d, Capsule3d, Cone, ConicalFrustum, Cuboid, Cylinder, Line3d, Measured3d,
        Plane3d, Polyline3d, Segment3d, Sphere, Torus,
    },
    Vec3,
};

impl ComputeMassProperties3d for Sphere {
    fn mass(&self, density: f32) -> f32 {
        self.volume() * density
    }

    fn angular_inertia(&self, mass: f32) -> Vec3 {
        let inertia = self.radius.powi(2) / 2.0 / 5.0 * mass;
        Vec3::splat(inertia)
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }
}

impl ComputeMassProperties3d for Cuboid {
    fn mass(&self, density: f32) -> f32 {
        self.volume() * density
    }

    fn angular_inertia(&self, mass: f32) -> Vec3 {
        let ix = self.half_size.x.powi(2) / 3.0;
        let iy = self.half_size.y.powi(2) / 3.0;
        let iz = self.half_size.z.powi(2) / 3.0;
        Vec3::new(iy + iz, ix + iz, ix + iy) * mass
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }
}

impl ComputeMassProperties3d for Cylinder {
    fn mass(&self, density: f32) -> f32 {
        self.volume() * density
    }

    fn angular_inertia(&self, mass: f32) -> Vec3 {
        let radius_squared = self.radius.powi(2);
        let height_squared = self.half_height.powi(2) * 4.0;
        let principal = radius_squared / 2.0;
        let off_principal = (radius_squared * 3.0 + height_squared) / 12.0;
        Vec3::new(off_principal, principal, off_principal) * mass
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }
}

impl ComputeMassProperties3d for Capsule3d {
    fn mass(&self, density: f32) -> f32 {
        let volume = self.radius * (std::f32::consts::PI * self.radius + 4.0 * self.half_length);
        volume * density
    }

    fn angular_inertia(&self, mass: f32) -> Vec3 {
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
        let density = mass / (cylinder_volume + sphere_volume);
        let cylinder_mass = cylinder_volume * density;
        let sphere_mass = sphere_volume * density;

        // Principal inertias
        let cylinder_inertia = cylinder.angular_inertia(1.0);
        let sphere_inertia = sphere.angular_inertia(1.0);

        // Total inertia
        let mut capsule_inertia = cylinder_inertia * cylinder_mass + sphere_inertia * sphere_mass;

        // Compensate for the hemispheres being away from the rotation axis using the parallel axis theorem.
        let extra = (cylinder_length.powi(2) * 0.25 + cylinder_length * self.radius * 3.0 / 8.0)
            * sphere_mass;
        capsule_inertia.x += extra;
        capsule_inertia.z += extra;

        capsule_inertia
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }

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
        let cylinder_inertia = cylinder.angular_inertia(1.0);
        let sphere_inertia = sphere.angular_inertia(1.0);

        // Total inertia
        let mut capsule_inertia = cylinder_inertia * cylinder_mass + sphere_inertia * sphere_mass;

        // Compensate for the hemispheres being away from the rotation axis using the parallel axis theorem.
        let extra = (cylinder_length.powi(2) * 0.25 + cylinder_length * self.radius * 3.0 / 8.0)
            * sphere_mass;
        capsule_inertia.x += extra;
        capsule_inertia.z += extra;

        MassProperties3d::new(cylinder_mass + sphere_mass, capsule_inertia, Vec3::ZERO)
    }
}

impl ComputeMassProperties3d for Cone {
    fn mass(&self, density: f32) -> f32 {
        self.volume() * density
    }

    fn angular_inertia(&self, mass: f32) -> Vec3 {
        let radius_squared = self.radius.powi(2);
        let height_squared = self.height.powi(2);
        let principal = radius_squared * 3.0 / 10.0;
        let off_principal = radius_squared * 3.0 / 20.0 + height_squared * 3.0 / 80.0;
        Vec3::new(off_principal, principal, off_principal) * mass
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::new(0.0, -self.height * 0.25, 0.0)
    }
}

impl ComputeMassProperties3d for ConicalFrustum {
    fn mass(&self, density: f32) -> f32 {
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

            frustum_volume * density
        }
    }

    fn angular_inertia(&self, mass: f32) -> Vec3 {
        if self.radius_top == self.radius_bottom {
            Cylinder::new(self.radius_top, self.height).angular_inertia(mass)
        } else {
            let max_radius = self.radius_top.max(self.radius_bottom);
            let radii_squared = self.radius_top.powi(2) + self.radius_bottom.powi(2);

            // The height of the cone formed by tapering the conical frustum to a tip.
            let cone_height =
                max_radius * (self.height / (self.radius_top - self.radius_bottom).abs());

            // Compute the angular inertia.
            // FIXME: `off_principal` is the same as it is for the cone. Is this correct?
            let principal = 1.0 / 12.0 * mass * (self.height.powi(2) + 3.0 * radii_squared);
            let off_principal = max_radius.powi(2) * 3.0 / 20.0 + cone_height.powi(2) * 3.0 / 80.0;

            Vec3::new(off_principal, principal, off_principal) * mass
        }
    }

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
            let principal = 1.0 / 12.0 * frustum_mass * (self.height.powi(2) + 3.0 * radii_squared);
            let off_principal = max_radius.powi(2) * 3.0 / 20.0 + cone1_height.powi(2) * 3.0 / 80.0;
            let inertia = Vec3::new(off_principal, principal, off_principal) * frustum_mass;

            MassProperties3d::new(frustum_mass, inertia, center_of_mass)
        }
    }
}

impl ComputeMassProperties3d for Torus {
    fn mass(&self, density: f32) -> f32 {
        self.volume() * density
    }

    fn angular_inertia(&self, mass: f32) -> Vec3 {
        // Reference: https://en.wikipedia.org/wiki/List_of_moments_of_inertia

        let major_radius_squared = self.major_radius.powi(2);
        let minor_radius_squared = self.minor_radius.powi(2);

        let principal = 0.25 * mass * (4.0 * major_radius_squared + 3.0 * minor_radius_squared);
        let off_principal =
            1.0 / 8.0 * mass * (4.0 * major_radius_squared + 5.0 * minor_radius_squared);
        Vec3::new(off_principal, principal, off_principal)
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }
}

impl ComputeMassProperties3d for Plane3d {
    fn mass(&self, _density: f32) -> f32 {
        0.0
    }

    fn angular_inertia(&self, _mass: f32) -> Vec3 {
        Vec3::ZERO
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }

    fn mass_properties(&self, _density: f32) -> MassProperties3d {
        MassProperties3d::ZERO
    }
}

impl ComputeMassProperties3d for Line3d {
    fn mass(&self, _density: f32) -> f32 {
        0.0
    }

    fn angular_inertia(&self, _mass: f32) -> Vec3 {
        Vec3::ZERO
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }

    fn mass_properties(&self, _density: f32) -> MassProperties3d {
        MassProperties3d::ZERO
    }
}

impl ComputeMassProperties3d for Segment3d {
    fn mass(&self, _density: f32) -> f32 {
        0.0
    }

    fn angular_inertia(&self, _mass: f32) -> Vec3 {
        Vec3::ZERO
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }

    fn mass_properties(&self, _density: f32) -> MassProperties3d {
        MassProperties3d::ZERO
    }
}

impl<const N: usize> ComputeMassProperties3d for Polyline3d<N> {
    fn mass(&self, _density: f32) -> f32 {
        0.0
    }

    fn angular_inertia(&self, _mass: f32) -> Vec3 {
        Vec3::ZERO
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }

    fn mass_properties(&self, _density: f32) -> MassProperties3d {
        MassProperties3d::ZERO
    }
}

impl ComputeMassProperties3d for BoxedPolyline3d {
    fn mass(&self, _density: f32) -> f32 {
        0.0
    }

    fn angular_inertia(&self, _mass: f32) -> Vec3 {
        Vec3::ZERO
    }

    fn center_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }

    fn mass_properties(&self, _density: f32) -> MassProperties3d {
        MassProperties3d::ZERO
    }
}
