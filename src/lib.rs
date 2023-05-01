mod caliper;

use crate::caliper::Caliper;
use geo::{prelude::*, Coord, Line};
use itertools::{izip, Itertools, MinMaxResult};
use nalgebra as na;
use noisy_float::prelude::*;
use std::{
    borrow::Borrow,
    f64::consts::{FRAC_PI_2, PI},
};

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

    let convex_hull_points = geo::MultiPoint::from(points2d.clone())
        .convex_hull()
        .exterior()
        .clone();
    let mut ch: Vec<geo::Point<f64>> = convex_hull_points.points().collect();
    ch.pop();

    let top = ch
        .iter()
        .enumerate()
        .fold(0, |idx, (i, p)| if ch[idx].y() < p.y() { i } else { idx });
    let bottom = ch
        .iter()
        .enumerate()
        .fold(0, |idx, (i, p)| if ch[idx].y() > p.y() { i } else { idx });
    let right = ch
        .iter()
        .enumerate()
        .fold(0, |idx, (i, p)| if ch[idx].x() < p.x() { i } else { idx });
    let left = ch
        .iter()
        .enumerate()
        .fold(0, |idx, (i, p)| if ch[idx].x() > p.x() { i } else { idx });

    let mut top_caliper = Caliper {
        convex_hull: ch.clone(),
        current_point: top,
        current_angle: PI,
    };
    let mut bottom_caliper = Caliper {
        convex_hull: ch.clone(),
        current_point: bottom,
        current_angle: 0.0,
    };
    let mut right_caliper = Caliper {
        convex_hull: ch.clone(),
        current_point: right,
        current_angle: PI * 0.5,
    };
    let mut left_caliper = Caliper {
        convex_hull: ch.clone(),
        current_point: left,
        current_angle: PI * 1.5,
    };

    let mut bboxs: Vec<[(f64, f64); 4]> = vec![];

    while bottom_caliper.current_angle < FRAC_PI_2 {
        let vertices = [
            bottom_caliper.intersection(&right_caliper),
            top_caliper.intersection(&right_caliper),
            top_caliper.intersection(&left_caliper),
            bottom_caliper.intersection(&left_caliper),
        ];
        bboxs.push(vertices);

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
        let loss: Vec<f64> = bboxs
            .iter()
            .map(|x| hausdorff_distance(x, &points2d))
            .collect();

        let MinMaxResult::MinMax(min_loss, max_loss) = loss.iter().copied().minmax() else {
            return None;
        };

        loss.iter()
            .map(|x| (x - min_loss) / (max_loss - min_loss))
            .collect()
    };

    let area_error: Vec<f64> = {
        let loss: Vec<f64> = bboxs.iter().map(bbox_area).collect();

        let MinMaxResult::MinMax(min_loss, max_loss) = loss.iter().copied().minmax() else {
            return None;
        };

        loss.iter()
            .map(|x| (x - min_loss) / (max_loss - min_loss))
            .collect()
    };

    let loss: Vec<f64> = haus_error
        .iter()
        .zip(area_error)
        .map(|(a, b)| a * 0.7 + b * 0.3)
        .collect();

    let (minimal_bbox, _) = izip!(&bboxs, &loss).min_by_key(|(_, &loss)| r64(loss))?;

    let rotation = {
        let from = na::Point2::new(minimal_bbox[0].0, minimal_bbox[0].1);
        let to = na::Point2::new(minimal_bbox[1].0, minimal_bbox[1].1);
        let vector = to - from;
        let yaw = vector.y.atan2(vector.x);
        na::UnitQuaternion::from_euler_angles(0.0, 0.0, yaw)
    };

    let z_minmax = points3d.iter().map(|p| p.z).map(r64).minmax();
    let MinMaxResult::MinMax(z_min, z_max) = z_minmax else {
        return None;
    };

    let center_x =
        (minimal_bbox[0].0 + minimal_bbox[1].0 + minimal_bbox[2].0 + minimal_bbox[3].0) / 4.0;
    let center_y =
        (minimal_bbox[0].1 + minimal_bbox[1].1 + minimal_bbox[2].1 + minimal_bbox[3].1) / 4.0;
    let center_z = ((z_max + z_min) / 2.0).raw();
    let size_x = ((minimal_bbox[1].0 - minimal_bbox[0].0).powf(2.0)
        + (minimal_bbox[1].1 - minimal_bbox[0].1).powf(2.0))
    .sqrt();
    let size_y = ((minimal_bbox[0].0 - minimal_bbox[3].0).powf(2.0)
        + (minimal_bbox[0].1 - minimal_bbox[3].1).powf(2.0))
    .sqrt();
    let size_z = (z_max - z_min).raw();

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
    let w = ((vertex[0].0 - vertex[3].0).powf(2.0) + (vertex[0].1 - vertex[3].1).powf(2.0)).sqrt();
    let h = ((vertex[1].0 - vertex[0].0).powf(2.0) + (vertex[1].1 - vertex[0].1).powf(2.0)).sqrt();
    w * h
}

fn point2bbox(vertex: &[(f64, f64); 4], point: &geo::Point<f64>) -> f64 {
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
        .map(|p| point2bbox(vertex, p))
        .map(r64)
        .sum::<R64>()
        .raw()
}
