# aurora

AURORA:
Adjudicate
Uncertain
Regions and
Output
Reliable
Annotations

## About

Aurora is a tool that performs biological sequence annotation adjudication.

## Installation

To build Aurora from source, you'll first need to install Rust and Cargo.
The easiest way to do that is to use [rustup](https://rustup.rs/).

Once that's done, you can then build Aurora:

    git clone https://github.com/TravisWheelerLab/aurora
    cd aurora/
    cargo build --release

You'll then find the compiled binary at: `target/release/aurora`

For example, try running:

    target/release/aurora -h

## Usage

The inputs to Aurora are sequence alignments and the matrices that were used to produce the sequence alignments.

For example:

    $ aurora alignments.caf matrices.mat

## License

Aurora is licensed under the BSD-3-Clause license.

See `LICENSE` for details.

## Authors

Jack Roddy - jroddy@arizona.edu
