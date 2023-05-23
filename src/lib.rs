use std::unreachable;

use glam::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Axis {
    X,
    Y,
    Z,
    Size,
}

pub type MortonCode = u64;

// The number of bits used for just the size component in the morton code
// Assuming we're using a size bit every 7th element, we have floor(64 / 7) = 9
pub const NUMBER_OF_BITS_FOR_SIZE: u32 = 6;
// The number of bits used in the size LUT.
pub const SIZE_LUT_NUMBER_OF_BITS: u32 = NUMBER_OF_BITS_FOR_SIZE * 2;

pub const NUMBER_OF_BITS_BETWEEN_SIZE_BITS: u32 = 6;

pub struct MortonCodeGenerator {
    pub morton_index_scale: f32,
    pub lut: [[MortonCode; 256]; 9],
    // In the case that we're using 9 bits for the size component,
    // this array will be 2^12 large or 4,096
    pub size_lut: [MortonCode; 1 << SIZE_LUT_NUMBER_OF_BITS],
    pub bits_per_axis: UVec3,
    pub shifts: [[u32; 64]; 3],
    pub offset: Vec3,
    pub multiplier: Vec3,
    pub size_multiplier: f32,
    pub axes: [Axis; 64],
}

impl MortonCodeGenerator {
    pub fn new(bounding_box_min: Vec3, bounding_box_max: Vec3) -> MortonCodeGenerator {
        let morton_index_scale = 1.0;
        // Subtact 0.01 to fit inside the bins and avoid an overflow???
        let scale: f32 = (1 << NUMBER_OF_BITS_FOR_SIZE) as f32 - 0.1;

        let mut size_lut: [MortonCode; 1 << SIZE_LUT_NUMBER_OF_BITS] =
            [0; 1 << SIZE_LUT_NUMBER_OF_BITS];
        // from 0 to 2^12
        for i in 0..1 << SIZE_LUT_NUMBER_OF_BITS {
            // This table will be montonically increasing - first 0, then 1, ..., then 2^|s| = 2^6 = 64
            // Supposedly faster to have a table to compute, turns a binary search of comparisons into 1 lookup
            let j: u32 = ((i as f32 / (1 << SIZE_LUT_NUMBER_OF_BITS) as f32).sqrt() * scale) as u32;
            size_lut[i] = expand_size_bits(j);
        }

        let size = bounding_box_max - bounding_box_min;
        let offset = bounding_box_min;
        let size_multiplier = ((1 << SIZE_LUT_NUMBER_OF_BITS) as f32 - 0.001) / size.length();

        let largest_split_bits = 64;

        let mut axis: Axis = Axis::X;
        let mut axes: [Axis; 64] = [Axis::X; 64];

        let mut diagonal = size.clone();

        for i in 0..64 {
            /*
            (1, 6)
            (2, 13)
            (3, 20)
            (4, 27)
            (5, 34)
            (6, 41)
            (7, 48)
            (8, 55)
            */
            // 6, 13, 20, 27, 34, 41, 48, 55
            // 7, 14, 21, 28, 35, 42, 49, 56
            let is_potential_size_bit_location = ((i + 1) % 7) == 0;
            // suppose number of size bits is 6 and i == 41
            // 6 + 1 > (i + 1) / 7
            // 6 + 1 > (41 + 1) / 7
            // 6 + 1 > 42/7
            // 6 + 1 > 6
            let is_not_after_last_size_bit = NUMBER_OF_BITS_FOR_SIZE + 1 > ((i + 1) / 7);
            let is_size_bit = is_potential_size_bit_location && is_not_after_last_size_bit;
            if is_size_bit {
                axes[i as usize] = Axis::Size;
            } else {
                if i < largest_split_bits {
                    axis = major_axis(diagonal);
                    axes[i as usize] = axis;
                    match axis {
                        Axis::X => diagonal.x *= 0.5,
                        Axis::Y => diagonal.y *= 0.5,
                        Axis::Z => diagonal.z *= 0.5,
                        Axis::Size => unreachable!(),
                    }
                } else {
                    // roundrobin axes starting with axis that is "after" the one last used
                    axis = match axis {
                        Axis::X => Axis::Y,
                        Axis::Y => Axis::Z,
                        Axis::Z => Axis::X,
                        Axis::Size => Axis::X,
                    };
                    axes[i as usize] = axis;
                }
            }
        }

        let mut shifts = [[0; 64]; 3];

        let mut bits_per_axis: UVec3 = UVec3::default();
        let mut multiplier = Vec3::default();
        let mut lut = [[0; 256]; 9];
        for axis_index in 0..3 {
            let axis = match axis_index {
                0 => Axis::X,
                1 => Axis::Y,
                2 => Axis::Z,
                _ => unreachable!(),
            };
            let mut last_bit: u32 = 0;
            let mut bits: u32 = 0;
            let axis_shifts = &mut shifts[axis_index];
            for i in 0..64 {
                if axes[i as usize] == axis {
                    axis_shifts[bits as usize] = i - last_bit;
                    last_bit = i;
                    bits += 1;
                }
            }
            match axis {
                Axis::X => bits_per_axis.x = bits,
                Axis::Y => bits_per_axis.y = bits,
                Axis::Z => bits_per_axis.z = bits,
                Axis::Size => unreachable!(),
            }
            axis_shifts[bits as usize] = 63 - last_bit;
            axis_shifts[0] = 0;
            for i in 0..256 {
                lut[axis_index as usize][i] = expand_bits(bits_per_axis, shifts, axis, i as u32);
                lut[axis_index + 3 as usize][i] =
                    expand_bits(bits_per_axis, shifts, axis, (i as u32) << 8);
                lut[axis_index + 6 as usize][i] =
                    expand_bits(bits_per_axis, shifts, axis, (i as u32) << 16);
            }
            let q = (1u64 << bits) as f32 - 0.5;
            match axis {
                Axis::X => multiplier.x = q / size.x,
                Axis::Y => multiplier.y = q / size.y,
                Axis::Z => multiplier.z = q / size.z,
                Axis::Size => unreachable!(),
            }
        }

        MortonCodeGenerator {
            morton_index_scale,
            lut,
            size_lut,
            bits_per_axis,
            shifts,
            offset,
            multiplier,
            size_multiplier,
            axes,
        }
    }

    pub fn expand_bits(&self, axis: Axis, v: u32) -> MortonCode {
        return expand_bits(self.bits_per_axis, self.shifts, axis, v);
    }

    pub fn code(&self, min: Vec3, max: Vec3) -> MortonCode {
        let p = (((min + max) * 0.5 - self.offset) * self.multiplier - 0.1).as_uvec3();
        let j = ((max - min).length() * self.size_multiplier) as u32;
        let b: u32 = 0xFF;
        // Note: Only coding first 3 bytes/24 bits, as such if number of bits for an axis is over 24 there
        // will be bugs (loss of precision).
        return self.lut[0][(p.x & b) as usize]
            | self.lut[1][(p.y & b) as usize]
            | self.lut[2][(p.z & b) as usize]
            | self.lut[3][((p.x >> 8) & b) as usize]
            | self.lut[4][((p.y >> 8) & b) as usize]
            | self.lut[5][((p.z >> 8) & b) as usize]
            | self.lut[6][((p.x >> 16) & b) as usize]
            | self.lut[7][((p.y >> 16) & b) as usize]
            | self.lut[8][((p.z >> 16) & b) as usize]
            | self.size_lut[j as usize];
    }
}

pub fn expand_size_bits(v: u32) -> MortonCode {
    let mut result: MortonCode = 0;
    for i in 0..NUMBER_OF_BITS_FOR_SIZE {
        // Select the 1st 2nd 3rd ... 6th bit in the size value
        // Starting with lowest bit and working our way up
        let a: MortonCode = (v & (1 << i)).into();
        // add that bit, shifted to the left 6/12/18/24/... times
        // 000 000 S 000 000 S 000 000 S 000 000 S 000 000 S 000 000 S 00000 00000 00000 00000 00
        // 6 * 7 = 42 bits
        // 64 - 6*7 = 22, number of bits remaining
        result += (a
            << ((i * NUMBER_OF_BITS_BETWEEN_SIZE_BITS)
                + (MortonCode::BITS
                    - NUMBER_OF_BITS_FOR_SIZE * (NUMBER_OF_BITS_BETWEEN_SIZE_BITS + 1))))
            as MortonCode;
    }
    return result;
}

pub fn expand_bits(bits_per_axis: UVec3, shifts: [[u32; 64]; 3], axis: Axis, v: u32) -> MortonCode {
    let mut morton_code = 0u64;
    let axis_index = match axis {
        Axis::X => 0,
        Axis::Y => 1,
        Axis::Z => 2,
        Axis::Size => unreachable!(),
    };
    let bits_per_axis = bits_per_axis[axis_index];
    let mask = 1 << (bits_per_axis - 1);
    let mut v_mut = v;
    // iterate through bits, lowest to highest
    for i in 0..bits_per_axis {
        morton_code <<= shifts[axis_index][i as usize];
        // If this bit is on, turn the lowest bit on the morton code on
        if v_mut & mask > 0 {
            morton_code += 1;
        }

        v_mut <<= 1;
    }
    morton_code <<= shifts[axis_index][bits_per_axis as usize];
    return morton_code;
}

pub fn major_axis(diagonal: Vec3) -> Axis {
    if diagonal.x >= diagonal.y {
        if diagonal.x >= diagonal.z {
            return Axis::X;
        } else {
            return Axis::Z;
        }
    } else {
        if diagonal.y >= diagonal.z {
            return Axis::Y;
        } else {
            return Axis::Z;
        }
    }
}

#[cfg(test)]
mod tests {
    use glam::*;

    use crate::{expand_size_bits, Axis, MortonCodeGenerator};

    const SCENE_BB_MIN: Vec3 = Vec3::splat(0.0);
    const SCENE_BB_MAX: Vec3 = Vec3::splat(1000.0);

    const SCENE_BB_SKEWED_XY_MIN: Vec3 = Vec3::splat(0.0);
    const SCENE_BB_SKEWED_XY_MAX: Vec3 = Vec3::new(1000.0, 1000.0, 0.001);

    const SCENE_BB_SKEWED_X_MIN: Vec3 = Vec3::splat(0.0);
    const SCENE_BB_SKEWED_X_MAX: Vec3 = Vec3::new(1000.0, 0.001, 0.001);

    const SCENE_BB_MIN_VEC: [f32; 3] = [0.0, -1000.0, 1000.0];
    const SCENE_BB_MAX_VEC: [f32; 3] = [1000.0, 0.0, 2000.0];

    const ALL_SIZE_BITS: u64 = 0b0000001000000100000010000001000000100000010000000000000000000000;

    // Simplest small code case
    #[test]
    fn test_smallest_code() {
        let gen = MortonCodeGenerator::new(SCENE_BB_MIN, SCENE_BB_MAX);
        let x = gen.code(Vec3::splat(0.0), Vec3::splat(0.0));
        assert_eq!(x, 0u64);
    }

    // When box is negative
    #[test]
    fn test_smallest_code_2() {
        let gen = MortonCodeGenerator::new(Vec3::splat(-1000.0), Vec3::splat(0.0));
        let x = gen.code(Vec3::splat(-1000.0), Vec3::splat(-1000.0));
        assert_eq!(x, 0u64);
    }

    //When box is offset positively
    #[test]
    fn test_smallest_code_3() {
        let gen = MortonCodeGenerator::new(Vec3::splat(1000.0), Vec3::splat(2000.0));
        let x = gen.code(Vec3::splat(1000.0), Vec3::splat(1000.0));
        assert_eq!(x, 0u64);
    }

    #[test]
    fn find_size() {
        let gen = MortonCodeGenerator::new(
            Vec3::splat(SCENE_BB_MIN_VEC[0]),
            Vec3::splat(SCENE_BB_MAX_VEC[0]),
        );
        let mut size_val = 63;
        let mut size_code: u64 = expand_size_bits(size_val);
        let mut size_code_change_points: Vec<i32> = vec![];
        for j in 0..500000 {
            let x = gen.code(
                Vec3::splat(SCENE_BB_MIN_VEC[0] + (j as f32) * 0.001),
                Vec3::splat(SCENE_BB_MAX_VEC[0] - (j as f32) * 0.001),
            );
            if x & ALL_SIZE_BITS != size_code {
                size_code_change_points.push(j);
                size_val -= 1;
                size_code = expand_size_bits(size_val);
                assert!(x & ALL_SIZE_BITS == size_code);
            }
        }
        let size_converted: Vec<f32> = size_code_change_points
            .iter()
            .map(|x| (*x as f32) * 0.001)
            .collect();
        // Checked manually by inputting data into a google sheet and finding a polynomial regression,
        // it looks correct but I haven't verified it exactly.
        println!("{:?}", size_converted);
        assert_eq!(size_code_change_points.len(), 63);
        assert_eq!(size_val, 0, "{}", size_val);
    }

    #[test]
    fn test_extremal_codes() {
        // Morton codes at smallest point
        // Morton codes at smallest size
        for i in 0..3 {
            let gen = MortonCodeGenerator::new(
                Vec3::splat(SCENE_BB_MIN_VEC[i]),
                Vec3::splat(SCENE_BB_MAX_VEC[i]),
            );
            let x = gen.code(
                Vec3::splat(SCENE_BB_MIN_VEC[i]),
                Vec3::splat(SCENE_BB_MIN_VEC[i]),
            );
            assert_eq!(x, 0u64);
        }

        // Morton codes at largest point
        // Morton codes at smallest size
        for i in 0..3 {
            let gen = MortonCodeGenerator::new(
                Vec3::splat(SCENE_BB_MIN_VEC[i]),
                Vec3::splat(SCENE_BB_MAX_VEC[i]),
            );
            let x = gen.code(
                Vec3::splat(SCENE_BB_MAX_VEC[i]),
                Vec3::splat(SCENE_BB_MAX_VEC[i]),
            );
            const LARGEST_MORTON_CODE_WITH_SMALL_SIZE: u64 =
                0b1111110111111011111101111110111111011111101111111111111111111111;
            assert_eq!(
                x, LARGEST_MORTON_CODE_WITH_SMALL_SIZE,
                "x was 
                {:#066b} 
                but was supposed to be 
                {:#066b}",
                x, LARGEST_MORTON_CODE_WITH_SMALL_SIZE
            );
        }

        // Morton codes in center
        // Morton codes at largest size
        for i in 0..3 {
            let gen = MortonCodeGenerator::new(
                Vec3::splat(SCENE_BB_MIN_VEC[i]),
                Vec3::splat(SCENE_BB_MAX_VEC[i]),
            );
            let x = gen.code(
                Vec3::splat(SCENE_BB_MIN_VEC[i]),
                Vec3::splat(SCENE_BB_MAX_VEC[i]),
            );
            assert_eq!(x & ALL_SIZE_BITS, ALL_SIZE_BITS);
        }
    }

    #[test]
    fn test_adaptive_axis_order() {
        let gen = MortonCodeGenerator::new(SCENE_BB_MIN, SCENE_BB_MAX);
        let a = gen.axes;
        let b: [Axis; 64] = string_to_axes(
            "
            XYZ XYZ S
            XYZ XYZ S
            XYZ XYZ S
            XYZ XYZ S
            XYZ XYZ S
            XYZ XYZ S
            XYZ XYZ XYZ XYZ XYZ XYZ XYZ X
            ",
        );
        assert_eq!(a, b);
    }

    #[test]
    fn test_adaptive_axis_order_skewed_xy() {
        let gen = MortonCodeGenerator::new(SCENE_BB_SKEWED_XY_MIN, SCENE_BB_SKEWED_XY_MAX);
        let a = gen.axes;
        // 19 XY splits before XYZ split
        // 2 ^ 19 ~= 500,000 (slightly greater)
        // 1000.0 / 0.001 = 1,000,000
        let b = string_to_axes(
            "
            XY XY XY S
            XY XY XY S
            XY XY XY S
            XY XY XY S
            XY XY XY S
            XY XY XY S
            XY XYZ XYZ XYZ XYZ XYZ XYZ XY
            ",
        );
        assert_eq!(a, b, "\n{},\n{}", axes_to_string(a), axes_to_string(b));
    }

    #[test]
    fn test_adaptive_axis_order_skewed_x() {
        let gen = MortonCodeGenerator::new(SCENE_BB_SKEWED_X_MIN, SCENE_BB_SKEWED_X_MAX);
        let a = gen.axes;
        // 19 X splits before XYZ split
        // 2 ^ 19 ~= 500,000 (slightly greater)
        // 1000.0 / 0.001 = 1,000,000
        // we do 1 more X split in XYZ itself so roundrobin is working properly
        let b = string_to_axes(
            "
            XXXXXX S
            XXXXXX S
            XXXXXX S
            X XYZ XY S
            Z XYZ XY S
            Z XYZ XY S
            Z XYZ XYZ XYZ XYZ XYZ XYZ XYZ
            ",
        );
        assert_eq!(a, b, "\n{},\n{}", axes_to_string(a), axes_to_string(b));
    }

    fn string_to_axes(s: &str) -> [Axis; 64] {
        let cleared = clear_string(s);
        assert_eq!(cleared.chars().count(), 64);
        let mut b: [Axis; 64] = [Axis::X; 64];
        for (i, c) in cleared.chars().enumerate() {
            b[i] = match c {
                'X' => Axis::X,
                'Y' => Axis::Y,
                'Z' => Axis::Z,
                'S' => Axis::Size,
                _ => panic!("Unsupported character {} in test!", c.escape_debug()),
            };
        }
        return b;
    }

    fn axes_to_string(a: [Axis; 64]) -> String {
        a.iter()
            .map(|axis| match axis {
                Axis::X => "X",
                Axis::Y => "Y",
                Axis::Z => "Z",
                Axis::Size => "S",
            })
            .collect()
    }

    fn clear_string(s: &str) -> String {
        s.replace("\n", "").replace(" ", "")
    }
}
