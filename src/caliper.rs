use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct Caliper {
    pub convex_hull: Vec<geo::Point<f64>>,
    pub current_index: usize,
    pub current_angle: f64,
}

impl Caliper {
    pub fn rotate(&mut self, angle: f64) {
        if angle == self.rotate_to_next() {
            let size = self.convex_hull.len();
            self.current_index = (self.current_index + 1) % size;
        }
        self.current_angle += angle;
    }

    pub fn x(&self) -> f64 {
        self.current_point().x()
    }

    pub fn y(&self) -> f64 {
        self.current_point().y()
    }

    pub fn slope(&self) -> f64 {
        self.current_angle.tan()
    }

    // Return the angle needed to be rotated from current angle to next angle
    pub fn rotate_to_next(&self) -> f64 {
        let size = self.convex_hull.len();
        let next_idx = (self.current_index + 1) % size;
        let angle = {
            let point = &self.convex_hull[next_idx];
            let t = (point.y() - self.y()).atan2(point.x() - self.x());

            if t > 0.0 {
                t
            } else {
                t + PI * 2.0
            }
        };

        if angle > self.current_angle {
            angle - self.current_angle
        } else {
            angle - self.current_angle + (PI * 2.0)
        }
    }

    pub fn intersection(&self, target: &Caliper) -> (f64, f64) {
        let m1 = self.slope();
        let m2 = target.slope();
        let x = match (m1.is_finite(), m2.is_finite()) {
            (true, true) => {
                let t =
                    self.y() - target.y() - self.slope() * self.x() + target.slope() * target.x();
                t / (target.slope() - self.slope())
            }
            (false, false | true) => self.x(),
            (true, false) => target.x(),
        };
        let y = match (m1 == 0.0, m2 == 0.0) {
            (false, false) => {
                if m1.is_finite() {
                    self.y() + self.slope() * (x - self.x())
                } else {
                    target.y() + target.slope() * (x - target.x())
                }
            }
            (true, false | true) => self.y(),
            (false, true) => target.y(),
        };
        (x, y)
    }

    fn current_point(&self) -> &geo::Point<f64> {
        &self.convex_hull[self.current_index]
    }
}
