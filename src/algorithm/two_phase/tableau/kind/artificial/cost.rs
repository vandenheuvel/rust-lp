use crate::algorithm::two_phase::tableau::inverse_maintenance::number_types as inverse_maintenance;
use std::ops::{Add, Mul};

pub enum Cost {
    One
}
//
// impl<F: inverse_maintenance::Field> Add<F> for Cost {
//     type Output = F;
//
//     fn add(self, rhs: F) -> Self::Output {
//         rhs + 1
//     }
// }
//
// impl<NZ: inverse_maintenance::NonZero> Mul<NZ> for Cost {
//     type Output = NZ;
//
//     fn mul(self, rhs: NZ) -> Self::Output {
//         rhs
//     }
// }
//
// impl<NZ: inverse_maintenance::NonZero> Mul<&NZ> for &Cost {
//     type Output = NZ;
//
//     fn mul(self, rhs: &NZ) -> Self::Output {
//         rhs.clone()
//     }
// }
