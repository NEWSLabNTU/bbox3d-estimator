use anyhow::Result;
use common_types as types;
use geo::prelude::*;
use nalgebra as na;
use noisy_float::prelude::*;
use serde::{Deserialize, Serialize};
use std::{borrow::Borrow, fs, io::prelude::*};

unzip_n::unzip_n!(2);

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct BBox3DConfig {
    pub record_cluster_data: bool,
    pub ransac_tolerance_outlier: R64,
    pub ransac_tolerance_filter: R64,
    pub ransac_iteration: usize,
}

/// Fit a rotated 3D box on a set of points.
pub fn get_rotated_bbox3d(
    points: impl IntoIterator<Item = impl Borrow<na::Point3<f64>>>,
    cluster_id: usize,
    config: impl AsRef<BBox3DConfig>,
) -> Result<types::BBox3D> {
    let bbox3d_config = config.as_ref();

    let data = points
        .into_iter()
        .map(|point| {
            let point = point.borrow();
            na::Point3::new(point.x, point.y, point.z)
        })
        .collect::<Vec<_>>();

    let (points3d, points2d) = data
        .into_iter()
        .map(|point3d| {
            let point2d = geo::Point::new(point3d.x, point3d.y);
            (point3d, point2d)
        })
        .unzip_n_vec();

    if bbox3d_config.record_cluster_data {
        let mut file = fs::OpenOptions::new()
            .write(true)
            .append(true)
            .create(true)
            .open("cluster-points.json5")?;

        let points3d: Vec<_> = points3d.iter().map(|p| vec![p.x, p.y, p.z]).collect();
        writeln!(file, "{:#?}: {:#?},", cluster_id, points3d)?;
    }

    let convex_hull_points = geo::MultiPoint::from(points2d.clone())
        .convex_hull()
        .exterior()
        .clone();
    let convex_hull_points: Vec<geo::Point<f64>> = convex_hull_points.points().collect();
    let mut ch = convex_hull_points;
    ch.pop();

    use rotating_caliper::*;

    mod rotating_caliper {
        #[derive(Debug)]
        pub struct Caliper {
            pub convex_hull: Vec<geo::Point<f64>>,
            pub current_point: i32,
            pub current_angle: f64,
        }

        impl Caliper {
            pub fn rotate(&mut self, angle: f64) {
                if angle == self.rotate_to_next() {
                    let size: i32 = self.convex_hull.len() as i32;
                    self.current_point = (self.current_point + 1) % size;
                }
                self.current_angle += angle;
            }

            pub fn x(&self) -> f64 {
                self.convex_hull[self.current_point as usize].x()
            }

            pub fn y(&self) -> f64 {
                self.convex_hull[self.current_point as usize].y()
            }

            pub fn slope(&self) -> f64 {
                self.current_angle.tan()
            }

            // Return the angle needed to be rotated from current angle to next angle
            pub fn rotate_to_next(&self) -> f64 {
                let size = self.convex_hull.len() as i32;
                let next_idx = (self.current_point + 1) % size;
                let angle = {
                    let _t = (self.convex_hull[next_idx as usize].y() - self.y())
                        .atan2(self.convex_hull[next_idx as usize].x() - self.x());
                    match _t > 0.0 {
                        true => _t,
                        false => _t + std::f64::consts::PI * 2.0,
                    }
                };
                match angle > self.current_angle {
                    true => angle - self.current_angle,
                    false => angle - self.current_angle + (std::f64::consts::PI * 2.0),
                }
            }

            pub fn intersection(&self, target: &Caliper) -> (f64, f64) {
                let m1 = self.slope();
                let m2 = target.slope();
                let x = match (m1 == std::f64::INFINITY, m2 == std::f64::INFINITY) {
                    (false, false) => {
                        let _t = self.y() - target.y() - self.slope() * self.x()
                            + target.slope() * target.x();
                        _t / (target.slope() - self.slope())
                    }
                    (true, false) => self.x(),
                    (false, true) => target.x(),
                    _ => self.x(),
                };
                let y = match (m1 == 0.0, m2 == 0.0) {
                    (false, false) => match m1 == std::f64::INFINITY {
                        false => self.y() + self.slope() * (x - self.x()),
                        true => target.y() + target.slope() * (x - target.x()),
                    },
                    (true, false) => self.y(),
                    (false, true) => target.y(),
                    _ => self.y(),
                };
                (x, y)
            }
        }
    }

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
        current_point: top as i32,
        current_angle: std::f64::consts::PI,
    };
    let mut bottom_caliper = Caliper {
        convex_hull: ch.clone(),
        current_point: bottom as i32,
        current_angle: 0.0,
    };
    let mut right_caliper = Caliper {
        convex_hull: ch.clone(),
        current_point: right as i32,
        current_angle: std::f64::consts::PI * 0.5,
    };
    let mut left_caliper = Caliper {
        convex_hull: ch.clone(),
        current_point: left as i32,
        current_angle: std::f64::consts::PI * 1.5,
    };

    let mut bboxs: Vec<[(f64, f64); 4]> = vec![];

    while bottom_caliper.current_angle < std::f64::consts::PI * 0.5 {
        let vertices = [
            bottom_caliper.intersection(&right_caliper),
            top_caliper.intersection(&right_caliper),
            top_caliper.intersection(&left_caliper),
            bottom_caliper.intersection(&left_caliper),
        ];
        bboxs.push(vertices.clone());
        let mut angle = top_caliper.rotate_to_next();
        angle = if bottom_caliper.rotate_to_next() < angle {
            bottom_caliper.rotate_to_next()
        } else {
            angle
        };
        angle = if right_caliper.rotate_to_next() < angle {
            right_caliper.rotate_to_next()
        } else {
            angle
        };
        angle = if left_caliper.rotate_to_next() < angle {
            left_caliper.rotate_to_next()
        } else {
            angle
        };
        top_caliper.rotate(angle);
        bottom_caliper.rotate(angle);
        right_caliper.rotate(angle);
        left_caliper.rotate(angle);
    }

    fn bbox_area(vertex: &[(f64, f64); 4]) -> f64 {
        let w =
            ((vertex[0].0 - vertex[3].0).powf(2.0) + (vertex[0].1 - vertex[3].1).powf(2.0)).sqrt();
        let h =
            ((vertex[1].0 - vertex[0].0).powf(2.0) + (vertex[1].1 - vertex[0].1).powf(2.0)).sqrt();
        w * h
    }

    fn point2bbox(vertex: &[(f64, f64); 4], p: &geo::Point<f64>) -> f64 {
        let line: Vec<_> = (0..4)
            .map(|idx| (vertex[idx], vertex[(idx + 1) % 4]))
            .map(|((x1, y1), (x2, y2))| {
                geo::Line::new(
                    geo::Coordinate { x: x1, y: y1 },
                    geo::Coordinate { x: x2, y: y2 },
                )
            })
            .collect();
        line.iter()
            .map(|l| p.euclidean_distance(l))
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

    let haus_error: Vec<f64> = {
        let loss: Vec<f64> = bboxs
            .clone()
            .iter()
            .map(|x| hausdorff_distance(x, &points2d))
            .collect();
        let max_loss = loss.iter().fold(std::f64::MIN, |a, &b| a.max(b));
        let min_loss = loss.iter().fold(std::f64::MAX, |a, &b| a.min(b));
        if max_loss - min_loss != 0.0 {
            loss.iter()
                .map(|x| (x - min_loss) / (max_loss - min_loss))
                .collect()
        } else {
            loss
        }
    };

    let area_error: Vec<f64> = {
        let loss: Vec<f64> = bboxs.clone().iter().map(|x| bbox_area(x)).collect();
        let max_loss = loss.iter().fold(std::f64::MIN, |a, &b| a.max(b));
        let min_loss = loss.iter().fold(std::f64::MAX, |a, &b| a.min(b));
        if max_loss - min_loss != 0.0 {
            loss.iter()
                .map(|x| (x - min_loss) / (max_loss - min_loss))
                .collect()
        } else {
            loss
        }
    };

    let loss: Vec<f64> = haus_error
        .iter()
        .zip(area_error)
        .map(|(a, b)| a * 0.7 + b * 0.3)
        .collect();

    let (idx, _) = loss.iter().enumerate().fold((0, std::f64::MAX), |min, x| {
        if &min.1 > x.1 {
            (x.0, *x.1)
        } else {
            min
        }
    });

    let minimal_bbox = bboxs[idx];

    let rotation = {
        let from = na::Point2::new(minimal_bbox[0].0, minimal_bbox[0].1);
        let to = na::Point2::new(minimal_bbox[1].0, minimal_bbox[1].1);
        let vector = to - from;
        let yaw = vector.y.atan2(vector.x);
        na::UnitQuaternion::from_euler_angles(0.0, 0.0, yaw)
    };

    let z_max = points3d.iter().map(|p| p.z).map(r64).max().unwrap().raw();
    let z_min = points3d.iter().map(|p| p.z).map(r64).min().unwrap().raw();

    let center_x =
        (minimal_bbox[0].0 + minimal_bbox[1].0 + minimal_bbox[2].0 + minimal_bbox[3].0) / 4.0;
    let center_y =
        (minimal_bbox[0].1 + minimal_bbox[1].1 + minimal_bbox[2].1 + minimal_bbox[3].1) / 4.0;
    let center_z = (z_max + z_min) / 2.0;
    let size_x = ((minimal_bbox[1].0 - minimal_bbox[0].0).powf(2.0)
        + (minimal_bbox[1].1 - minimal_bbox[0].1).powf(2.0))
    .sqrt();
    let size_y = ((minimal_bbox[0].0 - minimal_bbox[3].0).powf(2.0)
        + (minimal_bbox[0].1 - minimal_bbox[3].1).powf(2.0))
    .sqrt();
    let size_z = z_max - z_min;

    Ok(types::BBox3D {
        center_x,
        center_y,
        center_z,
        size_x,
        size_y,
        size_z,
        rotation: Some(types::protobuf_types::UnitQuaternion {
            x: rotation.i,
            y: rotation.j,
            z: rotation.k,
            w: rotation.w,
        }),
    })
}
