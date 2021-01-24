use std::num::{NonZeroI8, NonZeroU8};
use crate::data::number_types::traits::factorization::{NonzeroFactorizable, NonzeroFactorization};



macro_rules! impl_factorization_non_zero {
    ($t:ident, $f:ident) => {
        impl NonzeroFactorizable for $t {
            type Factor = $f;
            type Power = u8;
            // type Power = NonZeroU8;

            fn factorize(&self) -> NonzeroFactorization<Self::Factor, Self::Power> {
                unimplemented!()
            }
        }
    }
}

// impl_factorization_non_zero!(NonZeroI8, NonZeroU8);
impl_factorization_non_zero!(i8, u8);
impl_factorization_non_zero!(i16, u16);
impl_factorization_non_zero!(i32, u32);
impl_factorization_non_zero!(i64, u64);
