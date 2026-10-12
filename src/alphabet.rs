use std::str;

use phf::phf_map;

pub trait NucleotideByteUtils {
    fn to_utf8_string(&self) -> String;
    fn to_debug_utf8_string(&self) -> String;
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
}

impl NucleotideByteUtils for [u8] {
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
}

#[derive(Debug)]
pub enum NucleotideAlignmentType {
    Match,
    Transition,
    Transversion,
    Indel,
    Unknown,
}

impl NucleotideAlignmentType {
    pub fn from_pair(nucleotide_a: u8, nucleotide_b: u8) -> Self {
        // Place nucleotides in sorted order...
        if matches!(
            nucleotide_a,
            GAP_EXTEND_DIGITAL | GAP_OPEN_DIGITAL | PAD_DIGITAL
        ) || matches!(
            nucleotide_b,
            GAP_EXTEND_DIGITAL | GAP_OPEN_DIGITAL | PAD_DIGITAL
        ) {
            return Self::Indel;
        }

        if matches!(nucleotide_a, A_DIGITAL..=T_DIGITAL)
            && matches!(nucleotide_b, A_DIGITAL..=T_DIGITAL)
        {
            if nucleotide_a == nucleotide_b {
                return Self::Match;
            }

            return match (nucleotide_a, nucleotide_b) {
                (A_DIGITAL, G_DIGITAL) | (G_DIGITAL, A_DIGITAL) => Self::Transition,
                (C_DIGITAL, T_DIGITAL) | (T_DIGITAL, C_DIGITAL) => Self::Transition,
                _ => Self::Transversion,
            };
        }

        Self::Unknown
    }
}

impl NucleotideByteUtils for u8 {
    fn to_utf8_string(&self) -> String {
        ALIGNMENT_ALPHABET_STR[*self as usize].to_string()
    }

    fn to_debug_utf8_string(&self) -> String {
        DEBUG_ALIGNMENT_ALPHABET_STR[*self as usize].to_string()
    }
}

// Temp const for one third used in table below
const TRD: f64 = 1.0 / 3.0;

pub const NUCLEOTIDE_WEIGHTS: [[f64; 4]; 19] = [
    //A    C    G    T
    [1.0, 0.0, 0.0, 0.0],     // A
    [0.0, 1.0, 0.0, 0.0],     // C
    [0.0, 0.0, 1.0, 0.0],     // G
    [0.0, 0.0, 0.0, 1.0],     // T
    [0.0, TRD, TRD, TRD],     // B: C | G | T
    [TRD, 0.0, TRD, TRD],     // D: A | G | T
    [TRD, TRD, 0.0, TRD],     // H: A | C | T
    [0.0, 0.0, 0.5, 0.5],     // K: G | T
    [0.5, 0.5, 0.0, 0.0],     // M: A | C
    [0.25, 0.25, 0.25, 0.25], // N: A | C | G | T
    [0.5, 0.0, 0.5, 0.0],     // R: A | G
    [0.0, 0.5, 0.5, 0.0],     // S: C | G
    [TRD, TRD, TRD, 0.0],     // V: A | C | G
    [0.5, 0.0, 0.0, 0.5],     // W: A | T
    [0.0, 0.0, 0.0, 0.0],     // X: I think this is a masked base
    [0.0, 0.5, 0.0, 0.5],     // Y: C | T
    [0.0, 0.0, 0.0, 0.0],     // -: gap open
    [0.0, 0.0, 0.0, 0.0],     // -: gap extend
    [0.25, 0.25, 0.25, 0.25], // *: A | C | G | T (this isn't a IUPAC code, just an aurora thing)
];

pub const NUCLEOTIDE_ALPHABET_UTF8: [u8; 16] = [
    b'A', b'C', b'G', b'T', b'B', b'D', b'H', b'K', b'M', b'N', b'R', b'S', b'V', b'W', b'X', b'Y',
];

pub const NUCLEOTIDE_TO_COMPLEMENT: [u8; 16] = [
    T_DIGITAL, // A
    G_DIGITAL, // C
    C_DIGITAL, // G
    A_DIGITAL, // T
    V_DIGITAL, // B
    H_DIGITAL, // D
    D_DIGITAL, // H
    M_DIGITAL, // K
    K_DIGITAL, // M
    N_DIGITAL, // N
    Y_DIGITAL, // R
    S_DIGITAL, // S
    B_DIGITAL, // V
    W_DIGITAL, // W
    X_DIGITAL, // X
    R_DIGITAL, // Y
];

pub const ALIGNMENT_ALPHABET_UTF8: [u8; 19] = [
    b'A', b'C', b'G', b'T', b'B', b'D', b'H', b'K', b'M', b'N', b'R', b'S', b'V', b'W', b'X', b'Y',
    b'-', b'-', b'*',
];

pub const DEBUG_ALIGNMENT_ALPHABET_UTF8: [u8; 19] = [
    b'A', b'C', b'G', b'T', b'B', b'D', b'H', b'K', b'M', b'N', b'R', b'S', b'V', b'W', b'X', b'Y',
    b'-', b'+', b'*',
];

pub const ALIGNMENT_ALPHABET_STR: [&str; 19] = [
    "A", "C", "G", "T", "B", "D", "H", "K", "M", "N", "R", "S", "V", "W", "X", "Y", "-", "-", "*",
];

pub const DEBUG_ALIGNMENT_ALPHABET_STR: [&str; 19] = [
    "A", "C", "G", "T", "B", "D", "H", "K", "M", "N", "R", "S", "V", "W", "X", "Y", "-", "+", "*",
];

pub const PLUS_UTF8: u8 = "+".as_bytes()[0];
pub const FORWARD_SLASH_UTF8: u8 = "/".as_bytes()[0];
pub const DASH_UTF8: u8 = "-".as_bytes()[0];
pub const SPACE_UTF8: u8 = " ".as_bytes()[0];

pub const A_DIGITAL: u8 = 0;
pub const C_DIGITAL: u8 = 1;
pub const G_DIGITAL: u8 = 2;
pub const T_DIGITAL: u8 = 3;
pub const B_DIGITAL: u8 = 4;
pub const D_DIGITAL: u8 = 5;
pub const H_DIGITAL: u8 = 6;
pub const K_DIGITAL: u8 = 7;
pub const M_DIGITAL: u8 = 8;
pub const N_DIGITAL: u8 = 9;
pub const R_DIGITAL: u8 = 10;
pub const S_DIGITAL: u8 = 11;
pub const V_DIGITAL: u8 = 12;
pub const W_DIGITAL: u8 = 13;
pub const X_DIGITAL: u8 = 14;
pub const Y_DIGITAL: u8 = 15;
pub const GAP_OPEN_DIGITAL: u8 = 16;
pub const GAP_EXTEND_DIGITAL: u8 = 17;
pub const PAD_DIGITAL: u8 = 18;

pub const UTF8_TO_DIGITAL_NUCLEOTIDE: phf::Map<u8, u8> = phf_map! {
    // core
    b'A' => A_DIGITAL,
    b'a' => A_DIGITAL,

    b'C' => C_DIGITAL,
    b'c' => C_DIGITAL,

    b'G' =>  G_DIGITAL,
    b'g' => G_DIGITAL,

    b'T' =>  T_DIGITAL,
    b't' => T_DIGITAL,
    b'U' =>  T_DIGITAL,
    b'u' => T_DIGITAL,

    b'-' => GAP_OPEN_DIGITAL,
    b'+' => GAP_EXTEND_DIGITAL,

    // ambiguity
    b'B' => B_DIGITAL,
    b'b' => B_DIGITAL,

    b'D' => D_DIGITAL,
    b'd' => D_DIGITAL,

    b'H' => H_DIGITAL,
    b'h' => H_DIGITAL,

    b'K' =>  K_DIGITAL,
    b'k' => K_DIGITAL,

    b'M' =>  M_DIGITAL,
    b'm' => M_DIGITAL,

    b'N' =>  N_DIGITAL,
    b'n' => N_DIGITAL,

    b'R' =>  R_DIGITAL,
    b'r' => R_DIGITAL,

    b'S' =>  S_DIGITAL,
    b's' => S_DIGITAL,

    b'V' => V_DIGITAL,
    b'v' => V_DIGITAL,

    b'W' =>  W_DIGITAL,
    b'w' => W_DIGITAL,

    b'X' =>  X_DIGITAL,
    b'x' => X_DIGITAL,

    b'Y' =>  Y_DIGITAL,
    b'y' => Y_DIGITAL,
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
    "U" => T_DIGITAL,
    "u" => T_DIGITAL,

    "-" => GAP_OPEN_DIGITAL,
    "+" => GAP_EXTEND_DIGITAL,

    // ambiguity
    "B" => B_DIGITAL,
    "b" => B_DIGITAL,

    "D" => D_DIGITAL,
    "d" => D_DIGITAL,

    "H" => H_DIGITAL,
    "h" => H_DIGITAL,

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

    "V" => V_DIGITAL,
    "v" => V_DIGITAL,

    "W" =>  W_DIGITAL,
    "w" => W_DIGITAL,

    "X" =>  X_DIGITAL,
    "x" => X_DIGITAL,

    "Y" =>  Y_DIGITAL,
    "y" => Y_DIGITAL,
};
