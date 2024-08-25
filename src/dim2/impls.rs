use crate::Mass;

use super::{AngularInertia2d, ComputeMassProperties2d, MassProperties2d};
use bevy_math::{
    primitives::{
        Annulus, Arc2d, BoxedPolygon, BoxedPolyline2d, Capsule2d, Circle, CircularSector,
        CircularSegment, Ellipse, Line2d, Measured2d, Plane2d, Polygon, Polyline2d, Rectangle,
        RegularPolygon, Rhombus, Segment2d, Triangle2d,
    },
    Vec2,
};

impl ComputeMassProperties2d for Circle {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.area() * density)
    }

    #[inline]
    fn unit_angular_inertia(&self) -> AngularInertia2d {
        AngularInertia2d::new(self.radius.powi(2) / 2.0)
    }

    #[inline]
    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }
}

impl ComputeMassProperties2d for CircularSector {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.area() * density)
    }

    #[inline]
    fn unit_angular_inertia(&self) -> AngularInertia2d {
        AngularInertia2d::new(0.5 * self.radius().powi(4) * self.angle())
    }

    #[inline]
    fn center_of_mass(&self) -> Vec2 {
        let angle = self.angle();
        let y = 2.0 * self.radius() * angle.sin() / (3.0 * angle);
        Vec2::new(0.0, y)
    }
}

impl ComputeMassProperties2d for CircularSegment {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.area() * density)
    }

    #[inline]
    fn unit_angular_inertia(&self) -> AngularInertia2d {
        let angle = self.angle();
        let (sin, cos) = angle.sin_cos();
        AngularInertia2d::new(
            self.radius().powi(4) / 4.0 * (angle - sin + 2.0 / 3.0 * sin * (1.0 - cos) / 2.0),
        )
    }

    #[inline]
    fn center_of_mass(&self) -> Vec2 {
        let y = self.radius() * self.half_angle().sin().powi(3)
            / (6.0 * self.half_angle() - self.angle().sin());
        Vec2::new(0.0, y)
    }
}

impl ComputeMassProperties2d for Ellipse {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.area() * density)
    }

    #[inline]
    fn unit_angular_inertia(&self) -> AngularInertia2d {
        AngularInertia2d::new(self.half_size.length_squared() / 4.0)
    }

    #[inline]
    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }
}

impl ComputeMassProperties2d for Annulus {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.area() * density)
    }

    #[inline]
    fn unit_angular_inertia(&self) -> AngularInertia2d {
        AngularInertia2d::new(
            0.5 * (self.outer_circle.radius.powi(2) + self.inner_circle.radius.powi(2)),
        )
    }

    #[inline]
    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }
}

impl ComputeMassProperties2d for Triangle2d {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.area() * density)
    }

    #[inline]
    fn unit_angular_inertia(&self) -> AngularInertia2d {
        // Adapted from Box2D: https://github.com/erincatto/box2d/blob/411acc32eb6d4f2e96fc70ddbdf01fe5f9b16230/src/collision/b2_polygon_shape.cpp#L274

        let e1 = self.vertices[1] - self.vertices[0];
        let e2 = self.vertices[2] - self.vertices[0];

        let int_x2 = e1.x * e1.x + e2.x * e1.x + e2.x * e2.x;
        let int_y2 = e1.y * e1.y + e2.y * e1.y + e2.y * e2.y;

        AngularInertia2d::new((int_x2 + int_y2) / 6.0)
    }

    #[inline]
    fn center_of_mass(&self) -> Vec2 {
        (self.vertices[0] + self.vertices[1] + self.vertices[2]) / 3.0
    }

    #[inline]
    fn mass_properties(&self, density: f32) -> MassProperties2d {
        let area = self.area();
        let center_of_mass = self.center_of_mass();

        if area < f32::EPSILON {
            return MassProperties2d::new(Mass::ZERO, AngularInertia2d::ZERO, center_of_mass);
        }

        let mass = Mass::new(area * density);

        MassProperties2d::new(mass, self.angular_inertia(mass), center_of_mass)
    }
}

impl ComputeMassProperties2d for Rectangle {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.area() * density)
    }

    #[inline]
    fn unit_angular_inertia(&self) -> AngularInertia2d {
        AngularInertia2d::new(self.half_size.length_squared() / 3.0)
    }

    #[inline]
    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }
}

impl ComputeMassProperties2d for Rhombus {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.area() * density)
    }

    #[inline]
    fn unit_angular_inertia(&self) -> AngularInertia2d {
        AngularInertia2d::new(self.half_diagonals.length_squared() / 12.0)
    }

    #[inline]
    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }
}

impl ComputeMassProperties2d for RegularPolygon {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        Mass::new(self.area() * density)
    }

    #[inline]
    fn unit_angular_inertia(&self) -> AngularInertia2d {
        let half_external_angle = std::f32::consts::PI / self.sides as f32;
        AngularInertia2d::new(
            self.circumradius().powi(2) / 6.0 * (1.0 + 2.0 * half_external_angle.cos().powi(2)),
        )
    }

    #[inline]
    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }
}

impl ComputeMassProperties2d for Capsule2d {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        let area = self.radius * (std::f32::consts::PI * self.radius + 4.0 * self.half_length);
        Mass::new(area * density)
    }

    #[inline]
    fn unit_angular_inertia(&self) -> AngularInertia2d {
        // The rectangle and hemicircle parts
        let rectangle = Rectangle {
            half_size: Vec2::new(self.radius, self.half_length),
        };
        let rectangle_height = rectangle.half_size.y * 2.0;
        let circle = Circle::new(self.radius);

        // Areas
        let rectangle_area = rectangle.area();
        let circle_area = circle.area();

        // Masses
        let density = Mass::ONE / (rectangle_area + circle_area);
        let rectangle_mass = rectangle_area * density;
        let circle_mass = circle_area * density;

        // Principal inertias
        let rectangle_inertia = rectangle.angular_inertia(Mass::ONE);
        let circle_inertia = circle.angular_inertia(Mass::ONE);

        // Total inertia
        let mut capsule_inertia = rectangle_mass * rectangle_inertia + circle_mass * circle_inertia;

        // Compensate for the hemicircles being away from the rotation axis using the parallel axis theorem.
        capsule_inertia += AngularInertia2d::new(
            (rectangle_height.powi(2) * 0.25 + rectangle_height * self.radius * 3.0 / 8.0)
                * *circle_mass,
        );

        capsule_inertia
    }

    #[inline]
    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }

    #[inline]
    fn mass_properties(&self, density: f32) -> MassProperties2d {
        // The rectangle and hemicircle parts
        let rectangle = Rectangle {
            half_size: Vec2::new(self.radius, self.half_length),
        };
        let rectangle_height = rectangle.half_size.y * 2.0;
        let circle = Circle::new(self.radius);

        // Areas
        let rectangle_area = rectangle.area();
        let circle_area = circle.area();

        // Masses
        let rectangle_mass = rectangle_area * density;
        let circle_mass = circle_area * density;

        // Principal inertias
        let rectangle_inertia = rectangle.angular_inertia(Mass::ONE);
        let circle_inertia = circle.angular_inertia(Mass::ONE);

        // Total inertia
        let mut capsule_inertia = rectangle_inertia * rectangle_mass + circle_inertia * circle_mass;

        // Compensate for the hemicircles being away from the rotation axis using the parallel axis theorem.
        capsule_inertia += AngularInertia2d::new(
            (rectangle_height.powi(2) * 0.25 + rectangle_height * self.radius * 3.0 / 8.0)
                * circle_mass,
        );

        MassProperties2d::new(
            Mass::new(rectangle_mass + circle_mass),
            capsule_inertia,
            Vec2::ZERO,
        )
    }
}

impl<const N: usize> ComputeMassProperties2d for Polygon<N> {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        // The polygon is assumed to be convex.
        let area = convex_polygon_area(&self.vertices);
        Mass::new(area * density)
    }

    #[inline]
    fn unit_angular_inertia(&self) -> AngularInertia2d {
        convex_polygon_unit_angular_inertia(&self.vertices)
    }

    #[inline]
    fn center_of_mass(&self) -> Vec2 {
        convex_polygon_area_and_com(&self.vertices).1
    }

    #[inline]
    fn mass_properties(&self, density: f32) -> MassProperties2d {
        convex_polygon_mass_properties(&self.vertices, density)
    }
}

impl ComputeMassProperties2d for BoxedPolygon {
    #[inline]
    fn mass(&self, density: f32) -> Mass {
        // The polygon is assumed to be convex.
        let area = convex_polygon_area(&self.vertices);
        Mass::new(area * density)
    }

    #[inline]
    fn unit_angular_inertia(&self) -> AngularInertia2d {
        convex_polygon_unit_angular_inertia(&self.vertices)
    }

    #[inline]
    fn center_of_mass(&self) -> Vec2 {
        convex_polygon_area_and_com(&self.vertices).1
    }

    #[inline]
    fn mass_properties(&self, density: f32) -> MassProperties2d {
        convex_polygon_mass_properties(&self.vertices, density)
    }
}

#[inline]
fn convex_polygon_mass_properties(vertices: &[Vec2], density: f32) -> MassProperties2d {
    // The polygon is assumed to be convex.
    let (area, center_of_mass) = convex_polygon_area_and_com(vertices);

    if area < f32::EPSILON {
        return MassProperties2d::new(Mass::ZERO, AngularInertia2d::ZERO, center_of_mass);
    }

    // Initialize polygon inertia.
    let mut inertia = AngularInertia2d::ZERO;

    // Create a peekable iterator over the polygon vertices.
    let mut iter = vertices.iter().peekable();
    let first = **iter.peek().unwrap();

    // Iterate through vertices, computing the sum of the areas of triangles.
    // Each triangle is formed by the current vertex, next vertex, and the geometric center of the polygon.
    while let Some(vertex) = iter.next() {
        let triangle = Triangle2d::new(
            *vertex,
            iter.peek().copied().copied().unwrap_or(first),
            center_of_mass,
        );
        inertia += triangle.unit_angular_inertia() * triangle.area();
    }

    MassProperties2d::new(Mass::new(area * density), inertia * density, center_of_mass)
}

#[inline]
fn convex_polygon_unit_angular_inertia(vertices: &[Vec2]) -> AngularInertia2d {
    // The polygon is assumed to be convex.
    let (area, center_of_mass) = convex_polygon_area_and_com(vertices);

    if area < f32::EPSILON {
        return AngularInertia2d::ZERO;
    }

    // Initialize polygon inertia.
    let mut inertia = AngularInertia2d::ZERO;

    // Create a peekable iterator over the polygon vertices.
    let mut iter = vertices.iter().peekable();
    let first = **iter.peek().unwrap();

    // Iterate through vertices, computing the sum of the areas of triangles.
    // Each triangle is formed by the current vertex, next vertex, and the geometric center of the polygon.
    while let Some(vertex) = iter.next() {
        let triangle = Triangle2d::new(
            *vertex,
            iter.peek().copied().copied().unwrap_or(first),
            center_of_mass,
        );
        inertia += triangle.unit_angular_inertia() * triangle.area();
    }

    inertia / area
}

#[inline]
fn convex_polygon_area(vertices: &[Vec2]) -> f32 {
    let geometric_center =
        vertices.iter().fold(Vec2::ZERO, |acc, vtx| acc + *vtx) / vertices.len() as f32;

    // Initialize polygon area.
    let mut area = 0.0;

    // Create a peekable iterator over the polygon vertices.
    let mut iter = vertices.iter().peekable();
    let Some(first) = iter.peek().copied().copied() else {
        return 0.0;
    };

    // Iterate through vertices, computing the sum of the areas of triangles.
    // Each triangle is formed by the current vertex, next vertex, and the geometric center of the polygon.
    while let Some(vertex) = iter.next() {
        let (a, b, c) = (
            *vertex,
            iter.peek().copied().copied().unwrap_or(first),
            geometric_center,
        );
        let tri_area = Triangle2d::new(a, b, c).area();

        area += tri_area;
    }

    area
}

#[inline]
fn convex_polygon_area_and_com(vertices: &[Vec2]) -> (f32, Vec2) {
    let geometric_center =
        vertices.iter().fold(Vec2::ZERO, |acc, vtx| acc + *vtx) / vertices.len() as f32;

    // Initialize polygon area and center.
    let mut area = 0.0;
    let mut center = Vec2::ZERO;

    // Create a peekable iterator over the polygon vertices.
    let mut iter = vertices.iter().peekable();
    let Some(first) = iter.peek().copied().copied() else {
        return (0.0, Vec2::ZERO);
    };

    // Iterate through vertices, computing the sum of the areas and centers of triangles.
    // Each triangle is formed by the current vertex, next vertex, and the geometric center of the polygon.
    while let Some(vertex) = iter.next() {
        let (a, b, c) = (
            *vertex,
            iter.peek().copied().copied().unwrap_or(first),
            geometric_center,
        );
        let tri_area = Triangle2d::new(a, b, c).area();
        let tri_center = (a + b + c) / 3.0;

        area += tri_area;
        center += tri_center * tri_area;
    }

    if area < f32::EPSILON {
        (area, geometric_center)
    } else {
        (area, center / area)
    }
}

macro_rules! impl_zero_mass_properties_2d {
    ($($shape:ty),*) => {
        $(
            impl ComputeMassProperties2d for $shape {
                #[inline]
                fn mass(&self, _density: f32) -> Mass {
                    Mass::ZERO
                }

                #[inline]
                fn unit_angular_inertia(&self) -> AngularInertia2d {
                    AngularInertia2d::ZERO
                }

                #[inline]
                fn angular_inertia(&self, _mass: impl Into<Mass>) -> AngularInertia2d {
                    AngularInertia2d::ZERO
                }

                #[inline]
                fn center_of_mass(&self) -> Vec2 {
                    Vec2::ZERO
                }

                #[inline]
                fn mass_properties(&self, _density: f32) -> MassProperties2d {
                    MassProperties2d::ZERO
                }
            }
        )*
    };
}

impl_zero_mass_properties_2d!(Arc2d);
impl_zero_mass_properties_2d!(Plane2d);
impl_zero_mass_properties_2d!(Line2d);
impl_zero_mass_properties_2d!(Segment2d);
impl_zero_mass_properties_2d!(BoxedPolyline2d);

impl<const N: usize> ComputeMassProperties2d for Polyline2d<N> {
    #[inline]
    fn mass(&self, _density: f32) -> Mass {
        Mass::ZERO
    }

    #[inline]
    fn unit_angular_inertia(&self) -> AngularInertia2d {
        AngularInertia2d::ZERO
    }

    #[inline]
    fn angular_inertia(&self, _mass: impl Into<Mass>) -> AngularInertia2d {
        AngularInertia2d::ZERO
    }

    #[inline]
    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }

    #[inline]
    fn mass_properties(&self, _density: f32) -> MassProperties2d {
        MassProperties2d::ZERO
    }
}
