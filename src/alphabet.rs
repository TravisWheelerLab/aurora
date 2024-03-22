use phf::phf_map;

pub trait NucleotideByteUtils {
    fn to_utf8_string(&self) -> String;
    fn to_debug_utf8_string(&self) -> String;
    fn into_utf8_string(self) -> String;
}

impl NucleotideByteUtils for Vec<u8> {
    fn to_utf8_string(&self) -> String {
        String::from_utf8(
            self.iter()
                .map(|&b| ALIGNMENT_ALPHABET_UTF8[b as usize])
                .collect::<Vec<u8>>(),
        )
        .expect("failed to convert digital nucleotide byte vector to utf8 string")
    }

    fn to_debug_utf8_string(&self) -> String {
        String::from_utf8(
            self.iter()
                .map(|&b| DEBUG_ALIGNMENT_ALPHABET_UTF8[b as usize])
                .collect::<Vec<u8>>(),
        )
        .expect("failed to convert digital nucleotide byte vector to utf8 string")
    }

    fn into_utf8_string(self) -> String {
        String::from_utf8(
            self.into_iter()
                .map(|b| ALIGNMENT_ALPHABET_UTF8[b as usize])
                .collect::<Vec<u8>>(),
        )
        .expect("failed to convert digital nucleotide byte vector to utf8 string")
    }
}

impl NucleotideByteUtils for u8 {
    fn to_utf8_string(&self) -> String {
        ALIGNMENT_ALPHABET_STR[*self as usize].to_string()
    }

    fn to_debug_utf8_string(&self) -> String {
        DEBUG_ALIGNMENT_ALPHABET_STR[*self as usize].to_string()
    }

    fn into_utf8_string(self) -> String {
        ALIGNMENT_ALPHABET_STR[self as usize].to_string()
    }
}

pub const NUCLEOTIDE_WEIGHTS: [[f64; 4]; 15] = [
    //A    C    G    T
    [1.0, 0.0, 0.0, 0.0],     // A
    [0.0, 1.0, 0.0, 0.0],     // C
    [0.0, 0.0, 1.0, 0.0],     // G
    [0.0, 0.0, 0.0, 1.0],     // T
    [0.0, 0.0, 0.5, 0.5],     // K: G | T
    [0.5, 0.5, 0.0, 0.0],     // M: A | C
    [0.25, 0.25, 0.25, 0.25], // N: A | C | G | T
    [0.5, 0.0, 0.5, 0.0],     // R: A | G
    [0.0, 0.5, 0.5, 0.0],     // S: C | G
    [0.5, 0.0, 0.0, 0.5],     // W: A | T
    [0.0, 0.0, 0.0, 0.0],     // X: I think this is a masked base
    [0.0, 0.5, 0.0, 0.5],     // Y: C | T
    [0.0, 0.0, 0.0, 0.0],     // -: gap open
    [0.0, 0.0, 0.0, 0.0],     // -: gap extend
    [0.25, 0.25, 0.25, 0.25], // *: A | C | G | T (this isn't a IUPAC code, just an aurora thing)
];

pub const ALIGNMENT_ALPHABET_STR: [&str; 15] = [
    "A", "C", "G", "T", "K", "M", "N", "R", "S", "W", "X", "Y", "-", "-", "*",
];

pub const DEBUG_ALIGNMENT_ALPHABET_STR: [&str; 15] = [
    "A", "C", "G", "T", "K", "M", "N", "R", "S", "W", "X", "Y", "-", "+", "*",
];

pub const ALIGNMENT_ALPHABET_UTF8: [u8; 15] = [
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
    "-".as_bytes()[0],
    "-".as_bytes()[0],
    "*".as_bytes()[0],
];

pub const DEBUG_ALIGNMENT_ALPHABET_UTF8: [u8; 15] = [
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
    "-".as_bytes()[0],
    "+".as_bytes()[0],
    "*".as_bytes()[0],
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
pub const SPACE_UTF8: u8 = " ".as_bytes()[0];

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

    45u8 => GAP_OPEN_DIGITAL,
    43u8 => GAP_EXTEND_DIGITAL,

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

    "-" => GAP_OPEN_DIGITAL,
    "+" => GAP_EXTEND_DIGITAL,

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
