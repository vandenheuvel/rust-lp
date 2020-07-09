use crate::data::number_types::rational::Ratio;
use crate::data::number_types::rational::big::integer::signed::Signed;
use crate::data::number_types::rational::big::integer::non_zero_unsigned::NonZeroUnsigned;
use crate::data::number_types::rational::big::integer::non_zero_signed::NonZeroSigned;
use crate::data::number_types::rational::big::integer::unsigned::Unsigned;

pub mod integer;

pub type Big = Ratio<Signed, NonZeroUnsigned>;
pub type NonZero = Ratio<NonZeroSigned, NonZeroUnsigned>;
pub type NonNegative = Ratio<Unsigned, NonZeroUnsigned>;
pub type NonZeroNonNegative = Ratio<NonZeroUnsigned, NonZeroUnsigned>;
//
//
// impl PartialEq for Big {
//     fn eq(&self, other: &Self) -> bool {
//         unimplemented!()
//     }
// }
// impl Eq for Big {}
//
// impl PartialOrd for Big {
//     fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
//         Some(self.cmp(&other))
//     }
// }
// impl Ord for Big {
//     fn cmp(&self, other: &Self) -> Ordering {
//         unimplemented!()
//     }
// }
// impl Add for Big {
//     type Output = Self;
//
//     fn add(self, rhs: Self) -> Self::Output {
//         unimplemented!()
//     }
// }
// impl Add<&Self> for Big {
//     type Output = Self;
//
//     fn add(self, rhs: &Self) -> Self::Output {
//         unimplemented!()
//     }
// }
// impl AddAssign for Big {
//     fn add_assign(&mut self, rhs: Self) {
//         unimplemented!()
//     }
// }
// impl AddAssign<&Self> for Big {
//     fn add_assign(&mut self, rhs: &Self) {
//         unimplemented!()
//     }
// }
// impl Add<Rational64> for Big {
//     type Output = Self;
//
//     fn add(self, rhs: Rational64) -> Self::Output {
//         unimplemented!()
//     }
// }
//
// impl Sum for Big {
//     fn sum<I: Iterator<Item=Self>>(iter: I) -> Self {
//         unimplemented!()
//     }
// }
//
// impl Sub for Big {
//     type Output = Self;
//
//     fn sub(self, rhs: Self) -> Self::Output {
//         unimplemented!()
//     }
// }
// impl Sub<&Self> for Big {
//     type Output = Self;
//
//     fn sub(self, rhs: &Self) -> Self::Output {
//         unimplemented!()
//     }
// }
// impl SubAssign for Big {
//     fn sub_assign(&mut self, rhs: Self) {
//         unimplemented!()
//     }
// }
// impl SubAssign<&Self> for Big {
//     fn sub_assign(&mut self, rhs: &Self) {
//         unimplemented!()
//     }
// }
//
// impl Mul for Big {
//     type Output = Self;
//
//     fn mul(self, rhs: Self) -> Self::Output {
//         unimplemented!()
//     }
// }
// impl Mul<&Self> for Big {
//     type Output = Self;
//
//     fn mul(self, rhs: &Self) -> Self::Output {
//         unimplemented!()
//     }
// }
// impl MulAssign<Self> for Big {
//     fn mul_assign(&mut self, rhs: Self) {
//         unimplemented!()
//     }
// }
// impl MulAssign<&Self> for Big {
//     fn mul_assign(&mut self, rhs: &Self) {
//         unimplemented!()
//     }
// }
// impl MulAssign<Rational64> for Big {
//     fn mul_assign(&mut self, rhs: Rational64) {
//         unimplemented!()
//     }
// }
//
// impl Mul<Rational64> for Big {
//     type Output = Self;
//
//     fn mul(mut self, rhs: Rational64) -> Self::Output {
//         self *= rhs;
//         self
//     }
// }
//
// impl Div<Self> for Big {
//     type Output = Self;
//
//     fn div(self, rhs: Self) -> Self::Output {
//         unimplemented!()
//     }
// }
// impl Div<&Self> for Big {
//     type Output = Self;
//
//     fn div(self, rhs: &Self) -> Self::Output {
//         unimplemented!()
//     }
// }
// impl DivAssign<Self> for Big {
//     fn div_assign(&mut self, rhs: Self) {
//         unimplemented!()
//     }
// }
// impl DivAssign<&Self> for Big {
//     fn div_assign(&mut self, rhs: &Self) {
//         unimplemented!()
//     }
// }
//
// impl Neg for Big {
//     type Output = Self;
//
//     fn neg(self) -> Self::Output {
//         unimplemented!()
//     }
// }
//
// impl Zero for Big {
//     fn zero() -> Self {
//         unimplemented!()
//     }
//
//     fn is_zero(&self) -> bool {
//         unimplemented!()
//     }
// }
//
// impl One for Big {
//     fn one() -> Self {
//         unimplemented!()
//     }
// }
//
// impl Display for Big {
//     fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
//         unimplemented!()
//     }
// }
