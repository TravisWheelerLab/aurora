use phf::phf_map;

pub trait NucleotideByteUtils {
    fn as_utf8_string(&self) -> String;
}

impl NucleotideByteUtils for Vec<u8> {
    fn as_utf8_string(&self) -> String {
        self.iter()
            .map(|&b| ALIGNMENT_ALPHABET_STR[b as usize])
            .collect()
    }
}

pub const ALIGNMENT_ALPHABET_STR: [&str; 15] = [
    "A", "C", "G", "T", "K", "M", "N", "R", "S", "W", "X", "Y", "-", "+", "*",
];

pub const NUCLEOTIDE_ALPHABET_UTF8: [u8; 12] = [
    "A".as_bytes()[0],
    "C".as_bytes()[0],
    "G".as_bytes()[0],
    "T".as_bytes()[0],
    "K".as_bytes()[0],
    "M".as_bytes()[0],
    "N".as_bytes()[0],
    "R".as_bytes()[0],
    "S".as_bytes()[0],
    "W".as_bytes()[0],
    "X".as_bytes()[0],
    "Y".as_bytes()[0],
];

pub const PLUS_UTF8: u8 = "+".as_bytes()[0];
pub const FORWARD_SLASH_UTF8: u8 = "/".as_bytes()[0];
pub const DASH_UTF8: u8 = "-".as_bytes()[0];

pub const A_DIGITAL: u8 = 0;
pub const C_DIGITAL: u8 = 1;
pub const G_DIGITAL: u8 = 2;
pub const T_DIGITAL: u8 = 3;
pub const K_DIGITAL: u8 = 4;
pub const M_DIGITAL: u8 = 5;
pub const N_DIGITAL: u8 = 6;
pub const R_DIGITAL: u8 = 7;
pub const S_DIGITAL: u8 = 8;
pub const W_DIGITAL: u8 = 9;
pub const X_DIGITAL: u8 = 10;
pub const Y_DIGITAL: u8 = 11;
pub const GAP_OPEN_DIGITAL: u8 = 12;
pub const GAP_EXTEND_DIGITAL: u8 = 13;
pub const PAD_DIGITAL: u8 = 14;

pub const UTF8_TO_DIGITAL_NUCLEOTIDE: phf::Map<u8, u8> = phf_map! {
    // core
    65u8 => A_DIGITAL,
    97u8 => A_DIGITAL,

    67u8 => C_DIGITAL,
    99u8 => C_DIGITAL,

    71u8 =>  G_DIGITAL,
    103u8 => G_DIGITAL,

    84u8 =>  T_DIGITAL,
    116u8 => T_DIGITAL,

    // ambiguity
    75u8 =>  K_DIGITAL,
    107u8 => K_DIGITAL,

    77u8 =>  M_DIGITAL,
    109u8 => M_DIGITAL,

    78u8 =>  N_DIGITAL,
    110u8 => N_DIGITAL,

    82u8 =>  R_DIGITAL,
    114u8 => R_DIGITAL,

    83u8 =>  S_DIGITAL,
    115u8 => S_DIGITAL,

    87u8 =>  W_DIGITAL,
    119u8 => W_DIGITAL,

    88u8 =>  X_DIGITAL,
    120u8 => X_DIGITAL,

    89u8 =>  Y_DIGITAL,
    121u8 => Y_DIGITAL,
};

pub const STR_TO_DIGITAL_NUCLEOTIDE: phf::Map<&str, u8> = phf_map! {
    // core
    "A" => A_DIGITAL,
    "a" => A_DIGITAL,

    "C" => C_DIGITAL,
    "c" => C_DIGITAL,

    "G" =>  G_DIGITAL,
    "g" => G_DIGITAL,

    "T" =>  T_DIGITAL,
    "t" => T_DIGITAL,

    // ambiguity
    "K" =>  K_DIGITAL,
    "k" => K_DIGITAL,

    "M" =>  M_DIGITAL,
    "m" => M_DIGITAL,

    "N" =>  N_DIGITAL,
    "n" => N_DIGITAL,

    "R" =>  R_DIGITAL,
    "r" => R_DIGITAL,

    "S" =>  S_DIGITAL,
    "s" => S_DIGITAL,

    "W" =>  W_DIGITAL,
    "w" => W_DIGITAL,

    "X" =>  X_DIGITAL,
    "x" => X_DIGITAL,

    "Y" =>  Y_DIGITAL,
    "y" => Y_DIGITAL,
};
