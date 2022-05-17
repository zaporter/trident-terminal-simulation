use crate::*;
use std::ops::{Mul, Div};
use std::cmp::Ordering;

custom_derive! {
    #[derive(NewtypeFrom, NewtypeDeref, NewtypeDerefMut, NewtypeAdd, NewtypeSub, NewtypeAddAssign, 
        NewtypeMul(rug::Float), NewtypeMul(i32), Debug, Clone)]
    pub struct SimulationFloat(pub rug::Float);
}


impl SimulationFloat {
    pub fn new<T>(val:T) -> SimulationFloat
    where Float: Assign<T>
    {
        SimulationFloat(Float::with_val(256,val))
    }
}
impl Default for SimulationFloat{
    fn default() -> Self {
        SimulationFloat::new(0_f32)
    }
}
impl Mul for SimulationFloat {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self::new(self.0*rhs.0)
    }
}
impl Div for SimulationFloat {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        Self::new(self.0/rhs.0)
    }
}

impl PartialOrd for SimulationFloat {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
    }
    fn lt(&self, other: &Self) -> bool {
        self.0.lt(&other.0)
    }
    fn gt(&self, other: &Self) -> bool {
        self.0.gt(&other.0)
    }
}
impl PartialEq for SimulationFloat {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}
/*
impl From<rug::Float> for SimulationFloat {
    fn from(float : rug::Float) -> Self {
        SimulationFloat(float)
    }
}*/
