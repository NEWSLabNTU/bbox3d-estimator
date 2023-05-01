mod caliper;

use crate::caliper::Caliper;
use geo::{prelude::*, Coord, Line};
use itertools::{izip, Itertools};
use nalgebra as na;
use noisy_float::prelude::*;
use std::{
    borrow::Borrow,
    f64::consts::{FRAC_PI_2, PI},
};

const EPSILON: f64 = 1e-6;

#[derive(Debug, Clone)]
pub struct BBox3D {
    pub center_x: f64,
    pub center_y: f64,
    pub center_z: f64,
    pub size_x: f64,
    pub size_y: f64,
    pub size_z: f64,
    pub rotation: na::UnitQuaternion<f64>,
}

/// Fit a rotated 3D box on a set of points.
pub fn get_rotated_bbox3d(
    points: impl IntoIterator<Item = impl Borrow<na::Point3<f64>>>,
) -> Option<BBox3D> {
    let (points3d, points2d): (Vec<na::Point3<f64>>, Vec<geo::Point>) = points
        .into_iter()
        .map(|point3d| {
            let point3d = *point3d.borrow();
            let point2d = geo::Point::new(point3d.x, point3d.y);
            (point3d, point2d)
        })
        .unzip();

    let convex_hull: Vec<geo::Point<f64>> = {
        let convex_hull = geo::MultiPoint::from(points2d.clone()).convex_hull();
        let exterior = convex_hull.exterior();
        let num_points = exterior.0.len();
        exterior.points().take(num_points - 1).collect()
    };

    let ((bottom, _), (top, _)) = convex_hull
        .iter()
        .map(|point| point.y())
        .map(r64)
        .enumerate()
        .minmax_by_key(|&(_, y)| y)
        .into_option()?;
    let ((left, _), (right, _)) = convex_hull
        .iter()
        .map(|point| point.x())
        .map(r64)
        .enumerate()
        .minmax_by_key(|&(_, x)| x)
        .into_option()?;

    let mut top_caliper = Caliper {
        convex_hull: convex_hull.clone(),
        index: top,
        angle: PI,
    };
    let mut bottom_caliper = Caliper {
        convex_hull: convex_hull.clone(),
        index: bottom,
        angle: 0.0,
    };
    let mut right_caliper = Caliper {
        convex_hull: convex_hull.clone(),
        index: right,
        angle: FRAC_PI_2,
    };
    let mut left_caliper = Caliper {
        convex_hull,
        index: left,
        angle: PI + FRAC_PI_2,
    };

    let mut bboxes: Vec<[(f64, f64); 4]> = vec![];

    while bottom_caliper.angle < FRAC_PI_2 {
        let vertices = [
            bottom_caliper.intersection(&right_caliper),
            top_caliper.intersection(&right_caliper),
            top_caliper.intersection(&left_caliper),
            bottom_caliper.intersection(&left_caliper),
        ];
        bboxes.push(vertices);

        let angle = top_caliper.rotate_to_next();
        let angle = if bottom_caliper.rotate_to_next() < angle {
            bottom_caliper.rotate_to_next()
        } else {
            angle
        };
        let angle = if right_caliper.rotate_to_next() < angle {
            right_caliper.rotate_to_next()
        } else {
            angle
        };
        let angle = if left_caliper.rotate_to_next() < angle {
            left_caliper.rotate_to_next()
        } else {
            angle
        };

        top_caliper.rotate(angle);
        bottom_caliper.rotate(angle);
        right_caliper.rotate(angle);
        left_caliper.rotate(angle);
    }

    let haus_error: Vec<f64> = {
        let loss: Vec<f64> = bboxes
            .iter()
            .map(|bbox| hausdorff_distance(bbox, &points2d))
            .collect();

        let (min_loss, max_loss) = loss.iter().copied().minmax().into_option()?;
        let loss_range = max_loss - min_loss;
        if loss_range <= EPSILON {
            return None;
        }
        loss.iter().map(|x| (x - min_loss) / loss_range).collect()
    };
    let area_error: Vec<f64> = {
        let loss: Vec<f64> = bboxes.iter().map(bbox_area).collect();
        let (min_loss, max_loss) = loss.iter().copied().minmax().into_option()?;
        let loss_range = max_loss - min_loss;
        if loss_range <= EPSILON {
            return None;
        }
        loss.iter().map(|x| (x - min_loss) / loss_range).collect()
    };
    let loss: Vec<f64> = haus_error
        .iter()
        .zip(area_error)
        .map(|(a, b)| a * 0.7 + b * 0.3)
        .collect();

    let (minimal_bbox, _) = izip!(&bboxes, &loss).min_by_key(|(_, &loss)| r64(loss))?;

    let rotation = {
        let from = na::Point2::new(minimal_bbox[0].0, minimal_bbox[0].1);
        let to = na::Point2::new(minimal_bbox[1].0, minimal_bbox[1].1);
        let vector = to - from;
        let yaw = vector.y.atan2(vector.x);
        na::UnitQuaternion::from_euler_angles(0.0, 0.0, yaw)
    };

    let z_minmax = points3d.iter().map(|p| p.z).map(r64).minmax();
    let (z_min, z_max) = z_minmax.into_option()?;
    let z_range = (z_max - z_min).raw();
    if z_range <= EPSILON {
        return None;
    }

    let center_x =
        (minimal_bbox[0].0 + minimal_bbox[1].0 + minimal_bbox[2].0 + minimal_bbox[3].0) / 4.0;
    let center_y =
        (minimal_bbox[0].1 + minimal_bbox[1].1 + minimal_bbox[2].1 + minimal_bbox[3].1) / 4.0;
    let center_z = z_min.raw() + z_range / 2.0;
    let size_x = ((minimal_bbox[1].0 - minimal_bbox[0].0).powf(2.0)
        + (minimal_bbox[1].1 - minimal_bbox[0].1).powf(2.0))
    .sqrt();
    let size_y = ((minimal_bbox[0].0 - minimal_bbox[3].0).powf(2.0)
        + (minimal_bbox[0].1 - minimal_bbox[3].1).powf(2.0))
    .sqrt();
    let size_z = z_range;

    Some(BBox3D {
        center_x,
        center_y,
        center_z,
        size_x,
        size_y,
        size_z,
        rotation,
    })
}

fn bbox_area(vertex: &[(f64, f64); 4]) -> f64 {
    let w = distance(vertex[0], vertex[3]);
    let h = distance(vertex[0], vertex[1]);
    w * h
}

fn distance(lhs: (f64, f64), rhs: (f64, f64)) -> f64 {
    ((lhs.0 - rhs.0).powi(2) + (lhs.1 - rhs.1).powi(2)).sqrt()
}

fn point_to_bbox(vertex: &[(f64, f64); 4], point: &geo::Point<f64>) -> f64 {
    let line: Vec<_> = vertex
        .iter()
        .map(|&(x, y)| Coord { x, y })
        .circular_tuple_windows()
        .map(|(lhs, rhs)| Line::new(lhs, rhs))
        .collect();
    line.iter()
        .map(|line| point.euclidean_distance(line))
        .map(r64)
        .min()
        .unwrap()
        .raw()
}

// Update this from using all including points
fn hausdorff_distance(vertex: &[(f64, f64); 4], covex_hull: &[geo::Point<f64>]) -> f64 {
    covex_hull
        .iter()
        .map(|p| point_to_bbox(vertex, p))
        .map(r64)
        .sum::<R64>()
        .raw()
}
