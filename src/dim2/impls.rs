use super::{ComputeMassProperties2d, MassProperties2d};
use bevy_math::{
    primitives::{
        BoxedPolygon, BoxedPolyline2d, Capsule2d, Circle, Ellipse, Line2d, Measured2d, Plane2d,
        Polygon, Polyline2d, Rectangle, RegularPolygon, Segment2d, Triangle2d,
    },
    Vec2,
};

impl ComputeMassProperties2d for Circle {
    fn mass(&self, density: f32) -> f32 {
        self.area() * density
    }

    fn angular_inertia(&self, mass: f32) -> f32 {
        self.radius.powi(2) / 2.0 * mass
    }

    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }
}

impl ComputeMassProperties2d for Ellipse {
    fn mass(&self, density: f32) -> f32 {
        self.area() * density
    }

    fn angular_inertia(&self, mass: f32) -> f32 {
        self.half_size.length_squared() / 4.0 * mass
    }

    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }
}

impl ComputeMassProperties2d for Triangle2d {
    fn mass(&self, density: f32) -> f32 {
        self.area() * density
    }

    fn angular_inertia(&self, mass: f32) -> f32 {
        // Adapted from Box2D: https://github.com/erincatto/box2d/blob/411acc32eb6d4f2e96fc70ddbdf01fe5f9b16230/src/collision/b2_polygon_shape.cpp#L274

        let e1 = self.vertices[1] - self.vertices[0];
        let e2 = self.vertices[2] - self.vertices[0];

        let int_x2 = e1.x * e1.x + e2.x * e1.x + e2.x * e2.x;
        let int_y2 = e1.y * e1.y + e2.y * e1.y + e2.y * e2.y;

        (int_x2 + int_y2) / 6.0 * mass
    }

    fn center_of_mass(&self) -> Vec2 {
        (self.vertices[0] + self.vertices[1] + self.vertices[2]) / 3.0
    }

    fn mass_properties(&self, density: f32) -> MassProperties2d {
        let area = self.area();
        let center_of_mass = self.center_of_mass();

        if area < f32::EPSILON {
            return MassProperties2d::new(0.0, 0.0, center_of_mass);
        }

        let mass = area * density;

        MassProperties2d::new(mass, self.angular_inertia(mass), center_of_mass)
    }
}

impl ComputeMassProperties2d for Rectangle {
    fn mass(&self, density: f32) -> f32 {
        self.area() * density
    }

    fn angular_inertia(&self, mass: f32) -> f32 {
        (self.half_size.x.powi(2) + self.half_size.y.powi(2)) / 3.0 * mass
    }

    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }
}

impl ComputeMassProperties2d for RegularPolygon {
    fn mass(&self, density: f32) -> f32 {
        self.area() * density
    }

    fn angular_inertia(&self, mass: f32) -> f32 {
        let half_external_angle = std::f32::consts::PI / self.sides as f32;
        mass * self.circumradius().powi(2) / 6.0 * (1.0 + 2.0 * half_external_angle.cos().powi(2))
    }

    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }
}

impl ComputeMassProperties2d for Capsule2d {
    fn mass(&self, density: f32) -> f32 {
        let area = self.radius * (std::f32::consts::PI * self.radius + 4.0 * self.half_length);
        area * density
    }

    fn angular_inertia(&self, mass: f32) -> f32 {
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
        let density = mass / (rectangle_area + circle_area);
        let rectangle_mass = rectangle_area * density;
        let circle_mass = circle_area * density;

        // Principal inertias
        let rectangle_inertia = rectangle.angular_inertia(1.0);
        let circle_inertia = circle.angular_inertia(1.0);

        // Total inertia
        let mut capsule_inertia = rectangle_inertia * rectangle_mass + circle_inertia * circle_mass;

        // Compensate for the hemicircles being away from the rotation axis using the parallel axis theorem.
        capsule_inertia += (rectangle_height.powi(2) * 0.25
            + rectangle_height * self.radius * 3.0 / 8.0)
            * circle_mass;

        capsule_inertia
    }

    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }

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
        let rectangle_inertia = rectangle.angular_inertia(1.0);
        let circle_inertia = circle.angular_inertia(1.0);

        // Total inertia
        let mut capsule_inertia = rectangle_inertia * rectangle_mass + circle_inertia * circle_mass;

        // Compensate for the hemicircles being away from the rotation axis using the parallel axis theorem.
        capsule_inertia += (rectangle_height.powi(2) * 0.25
            + rectangle_height * self.radius * 3.0 / 8.0)
            * circle_mass;

        MassProperties2d::new(rectangle_mass + circle_mass, capsule_inertia, Vec2::ZERO)
    }
}

impl<const N: usize> ComputeMassProperties2d for Polygon<N> {
    fn mass(&self, density: f32) -> f32 {
        BoxedPolygon::new(self.vertices).mass(density)
    }

    fn angular_inertia(&self, mass: f32) -> f32 {
        BoxedPolygon::new(self.vertices).angular_inertia(mass)
    }

    fn center_of_mass(&self) -> Vec2 {
        BoxedPolygon::new(self.vertices).center_of_mass()
    }

    fn mass_properties(&self, density: f32) -> MassProperties2d {
        BoxedPolygon::new(self.vertices).mass_properties(density)
    }
}

impl ComputeMassProperties2d for BoxedPolygon {
    fn mass(&self, density: f32) -> f32 {
        // The polygon is assumed to be convex.
        let area = convex_polygon_area(self);
        area * density
    }

    fn angular_inertia(&self, mass: f32) -> f32 {
        // The polygon is assumed to be convex.
        let (area, center_of_mass) = convex_polygon_area_and_com(self);

        if area < f32::EPSILON {
            return 0.0;
        }

        // Initialize polygon inertia.
        let mut inertia = 0.0;

        // Create a peekable iterator over the polygon vertices.
        let mut iter = self.vertices.iter().peekable();
        let first = **iter.peek().unwrap();

        // Iterate through vertices, computing the sum of the areas of triangles.
        // Each triangle is formed by the current vertex, next vertex, and the geometric center of the polygon.
        while let Some(vertex) = iter.next() {
            let triangle = Triangle2d::new(
                *vertex,
                iter.peek().copied().copied().unwrap_or(first),
                center_of_mass,
            );
            inertia += triangle.angular_inertia(1.0) * triangle.area();
        }

        let density = mass / area;
        inertia * density
    }

    fn center_of_mass(&self) -> Vec2 {
        convex_polygon_area_and_com(self).1
    }

    fn mass_properties(&self, density: f32) -> MassProperties2d {
        // The polygon is assumed to be convex.
        let (area, center_of_mass) = convex_polygon_area_and_com(self);

        if area < f32::EPSILON {
            return MassProperties2d::new(0.0, 0.0, center_of_mass);
        }

        // Initialize polygon inertia.
        let mut inertia = 0.0;

        // Create a peekable iterator over the polygon vertices.
        let mut iter = self.vertices.iter().peekable();
        let first = **iter.peek().unwrap();

        // Iterate through vertices, computing the sum of the areas of triangles.
        // Each triangle is formed by the current vertex, next vertex, and the geometric center of the polygon.
        while let Some(vertex) = iter.next() {
            let triangle = Triangle2d::new(
                *vertex,
                iter.peek().copied().copied().unwrap_or(first),
                center_of_mass,
            );
            inertia += triangle.angular_inertia(1.0) * triangle.area();
        }

        MassProperties2d::new(area * density, inertia * density, center_of_mass)
    }
}

fn convex_polygon_area(polygon: &BoxedPolygon) -> f32 {
    let geometric_center = polygon
        .vertices
        .iter()
        .fold(Vec2::ZERO, |acc, vtx| acc + *vtx)
        / polygon.vertices.len() as f32;

    // Initialize polygon area.
    let mut area = 0.0;

    // Create a peekable iterator over the polygon vertices.
    let mut iter = polygon.vertices.iter().peekable();
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

fn convex_polygon_area_and_com(polygon: &BoxedPolygon) -> (f32, Vec2) {
    let geometric_center = polygon
        .vertices
        .iter()
        .fold(Vec2::ZERO, |acc, vtx| acc + *vtx)
        / polygon.vertices.len() as f32;

    // Initialize polygon area and center.
    let mut area = 0.0;
    let mut center = Vec2::ZERO;

    // Create a peekable iterator over the polygon vertices.
    let mut iter = polygon.vertices.iter().peekable();
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

impl ComputeMassProperties2d for Plane2d {
    fn mass(&self, _density: f32) -> f32 {
        0.0
    }

    fn angular_inertia(&self, _mass: f32) -> f32 {
        0.0
    }

    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }

    fn mass_properties(&self, _density: f32) -> MassProperties2d {
        MassProperties2d::ZERO
    }
}

impl ComputeMassProperties2d for Line2d {
    fn mass(&self, _density: f32) -> f32 {
        0.0
    }

    fn angular_inertia(&self, _mass: f32) -> f32 {
        0.0
    }

    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }

    fn mass_properties(&self, _density: f32) -> MassProperties2d {
        MassProperties2d::ZERO
    }
}

impl ComputeMassProperties2d for Segment2d {
    fn mass(&self, _density: f32) -> f32 {
        0.0
    }

    fn angular_inertia(&self, _mass: f32) -> f32 {
        0.0
    }

    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }

    fn mass_properties(&self, _density: f32) -> MassProperties2d {
        MassProperties2d::ZERO
    }
}

impl<const N: usize> ComputeMassProperties2d for Polyline2d<N> {
    fn mass(&self, _density: f32) -> f32 {
        0.0
    }

    fn angular_inertia(&self, _mass: f32) -> f32 {
        0.0
    }

    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }

    fn mass_properties(&self, _density: f32) -> MassProperties2d {
        MassProperties2d::ZERO
    }
}

impl ComputeMassProperties2d for BoxedPolyline2d {
    fn mass(&self, _density: f32) -> f32 {
        0.0
    }

    fn angular_inertia(&self, _mass: f32) -> f32 {
        0.0
    }

    fn center_of_mass(&self) -> Vec2 {
        Vec2::ZERO
    }

    fn mass_properties(&self, _density: f32) -> MassProperties2d {
        MassProperties2d::ZERO
    }
}
