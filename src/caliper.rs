#[derive(Debug, Clone)]
pub struct Caliper {
    pub convex_hull: Vec<geo::Point<f64>>,
    pub current_point: usize,
    pub current_angle: f64,
}

impl Caliper {
    pub fn rotate(&mut self, angle: f64) {
        if angle == self.rotate_to_next() {
            let size = self.convex_hull.len();
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
        let size = self.convex_hull.len();
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
                let _t =
                    self.y() - target.y() - self.slope() * self.x() + target.slope() * target.x();
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
