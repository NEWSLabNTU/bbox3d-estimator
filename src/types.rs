use nalgebra as na;

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
