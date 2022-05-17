use crate::*;
use std::ops::{Add, Sub, AddAssign,Mul, Div};

#[derive(Debug, Clone, Default)]
pub struct Vec3 {
    pub x:SimulationFloat,
    pub y:SimulationFloat,
    pub z:SimulationFloat,
}


impl Vec3 {
    pub fn new<T>(x:T, y:T, z:T) -> Vec3 
        where Float:Assign<T>
    {
        Vec3{
            x:SimulationFloat::new(x),
            y:SimulationFloat::new(y),
            z:SimulationFloat::new(z),
        }

    }
    pub fn modulus(&self) -> SimulationFloat {
        let mut sum = SimulationFloat::new(0_f32);
        sum += self.x.clone() * self.x.clone();
        sum += self.y.clone() * self.y.clone();
        sum += self.z.clone() * self.z.clone();
        sum.sqrt_mut();
        sum
    }

    pub fn unit(&self) -> Vec3 {
        let ret = self.clone();
        let modulus = ret.modulus();
        ret / modulus
    }
}

impl Mul<SimulationFloat> for Vec3 {
    type Output = Self;
    fn mul(self, rhs: SimulationFloat )->Self::Output {
        Self {
            x: self.x * rhs.clone(),
            y: self.y * rhs.clone(),
            z: self.z * rhs,
        }
        
    }
}
impl Div<SimulationFloat> for Vec3 {
    type Output = Self;
    fn div(self, rhs: SimulationFloat )->Self::Output {
        Self {
            x: self.x / rhs.clone(),
            y: self.y / rhs.clone(),
            z: self.z / rhs,
        }
        
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}
impl AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }

}
impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}
impl PartialEq for Vec3 {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y && self.z == other.z
    }
}

#[cfg(test)]
mod tests{
    use crate::*;
    #[test]
    fn vector_unit_test() {
        let vec_a = Vec3::new(-1_f32,5_f32,6_f32);
        let vec_b = Vec3::new(-2_f32,10_f32,12_f32);
        assert_eq!(vec_a.unit(), vec_b.unit());
        // This should be [0.127, 0.635, 0.762]
    }
    #[test]
    fn vector_magnitude_test() {
        let vec_a = Vec3::new(3_f32,4_f32,0_f32);
        assert_eq!(vec_a.modulus(), SimulationFloat::new(5_f32));
    }

}

