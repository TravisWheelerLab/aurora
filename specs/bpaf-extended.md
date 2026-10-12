# BPAF — Binary Pairwise Alignment Format, Version 0 with Sequence Letter Extension

This document is the on-disk specification for **BPAF version 0**, with an extension
that supports storing all sequence letters for sequences stored as alignments
in the file. It is implemented in aurora. It is written for implementers who need
to read or write BPAF files in any language, independent of aurora's implementation
in rust.

BPAF stores a stream of pairwise alignment records (query vs. target) plus a
small set of interned string tables, using variable-length integers to keep
files compact. An optional extension flag can be enabled to store 2 additonal
tables that store sequence letters for both the target and query sequences,
allowing for alignment sequences to be reconstructed entirely within the
BPAF format.

---

## 1. Conventions

| Concept             | Rule                                                                          |
| ------------------- | ----------------------------------------------------------------------------- |
| Byte order          | All fixed-width numbers are **little-endian**.                                |
| Integers (variable) | **ULEB128** unsigned varint (see §2.1).                                       |
| Signed integers     | **ZigZag** mapped to a ULEB128 (see §2.2). Used only inside the CIGAR.        |
| Floats              | IEEE-754, little-endian: `float32` (4 bytes) or `float64` (8 bytes) as noted. |
| Strings             | UTF-8, length-prefixed with a ULEB128 byte count. No NUL terminator.          |
| Coordinates         | **1-based, fully closed** intervals.                                          |
| Offsets             | Absolute byte offsets from the start of the file.                             |

### 2.1 ULEB128 (unsigned varint)

An unsigned integer is emitted 7 bits at a time, least-significant group first.
Each byte carries 7 payload bits in its low bits; the high bit (`0x80`) is set
on every byte **except the last**.

```
encode(x):                       decode(buf, i):
  while x >= 0x80:                 shift, val = 0, 0
    emit((x & 0x7F) | 0x80)        loop:
    x >>= 7                          b = buf[i]; i += 1
  emit(x & 0x7F)                     val |= (b & 0x7F) << shift
                                     if (b & 0x80) == 0: return val, i
                                     shift += 7
```

### 2.2 ZigZag (signed → unsigned)

Signed values are folded so small magnitudes stay small, then ULEB128-encoded:

```
zigzag(n)  = (n << 1)         if n >= 0
             (-n << 1) - 1    if n < 0
unzigzag(z) = (z >> 1) ^ -(z & 1)
```

so `0→0, -1→1, 1→2, -2→3, 2→4, …`.

---

## 2. File layout

A BPAF file has four parts, in this order:

```
┌───────────────────────────────────────────────----------┐
│ 1. Header          (fixed, 7 bytes)                     │
├──────────────────────────────────────────────----------─┤
│ 2. Record stream   (variable, 0..N records)             │
├───────────────────────────────────────────────----------┤
│ 3. String tables   (qids, tids, mats, [qseqs], [tseqs]) │  ← starts at `tables_start`
├───────────────────────────────────────────────----------┤
│ 4. Footer          (fixed, 13 bytes)                    │
└───────────────────────────────────────────────----------┘
```

The record stream occupies the bytes between the end of the header (offset `7`)
and `tables_start`. Records have no count field; a reader iterates them until it
reaches `tables_start`, which it learns from the footer.

---

## 3. Header (7 bytes)

| Offset | Size | Field        | Value (v0)                      | Meaning                       |
| ------ | ---- | ------------ | ------------------------------- | ----------------------------- |
| 0      | 5    | `MAGIC`      | `42 50 41 46 01` (`"BPAF\x01"`) | File signature.               |
| 5      | 1    | `version`    | `0x00`                          | Format version (see §7).      |
| 6      | 1    | `extensions` | `0x00`                          | Extensions bitfield (see §7). |

A reader **must** verify `MAGIC`, then read `version`. If `version` is not a
value the reader implements, it **must** refuse the file (see §7). The
`extensions` byte is a bitfield of optional, backward-compatible additions; a
reader may ignore bits it does not recognize.

> The `\x01` inside `MAGIC` is part of the fixed signature and is **not** the
> format version. The format version is the separate `version` byte at offset 5.

# 3.1 Sequence Extension

If the least significant bit is set in the extensions byte (bit `0x01`), the
file is assumed to contain two additional tables storing sequence letters.
These are described in §6.2.

---

## 4. Record stream

Each record is a ULEB128 **length prefix** followed by exactly that many bytes
of record body. The length prefix lets a reader skip a whole record without
decoding it (this is how [`utils/count_alignments.py`](../utils/count_alignments.py)
counts records cheaply).

```
record := uleb(reclen) , body[reclen]
```

### 4.1 Record body

Fields appear in this exact order:

| #   | Field             | Encoding | Notes                                                                |
| --- | ----------------- | -------- | -------------------------------------------------------------------- |
| 1   | `qidx`            | ULEB128  | Index into the **qids** table → `query_id`.                          |
| 2   | `tidx`            | ULEB128  | Index into the **tids** table → `target_id`.                         |
| 3   | `q_start`         | ULEB128  | Query start (1-based).                                               |
| 4   | `q_len`           | ULEB128  | Query span length. `query_end = q_start + q_len - 1`.                |
| 5   | `t_start`         | ULEB128  | Target start (1-based), the **lower** of the two target coordinates. |
| 6   | `t_len`           | ULEB128  | Target span length. `target_end = t_start + t_len - 1`.              |
| 7   | `midx`            | ULEB128  | Index into the **mats** table → `scoring_system`.                    |
| 8   | `flags`           | 1 byte   | Bitfield (see §4.2).                                                 |
| 9   | _optional fields_ | varies   | Present in flag order (see §4.3).                                    |
| 10  | `cigar`           | see §5   | Compact alignment CIGAR.                                             |

Target coordinates are always stored ascending (`t_start ≤ t_end`); the strand
is conveyed separately by the `ORIENT_C` flag. Writers that receive
`target_start > target_end` must swap them before writing and set `ORIENT_C`
according to the aligned strand.

### 4.2 `flags` byte

| Bit | Mask   | Name           | Meaning                                                |
| --- | ------ | -------------- | ------------------------------------------------------ |
| 0   | `0x01` | `ORIENT_C`     | Target aligns on the `-` strand relative to the query. |
| 1   | `0x02` | `HAS_SCORE`    | A `score` field is present (§4.3).                     |
| 2   | `0x04` | `HAS_EVAL`     | An `evalue` field is present.                          |
| 3   | `0x08` | `HAS_DIV`      | A `divergence` field is present.                       |
| 4   | `0x10` | `HAS_BITSCORE` | A `bit_score` field is present.                        |
| 5–7 | —      | reserved       | Must be 0 in v0.                                       |

### 4.3 Optional fields

For each set flag, the corresponding field is appended **in ascending bit
order** immediately after the `flags` byte (i.e. score, then evalue, then
divergence, then bit_score). A field is absent entirely when its flag is clear —
no placeholder bytes are written.

| Flag           | Field        | Encoding     | Size |
| -------------- | ------------ | ------------ | ---- |
| `HAS_SCORE`    | `score`      | `float32` LE | 4    |
| `HAS_EVAL`     | `evalue`     | `float64` LE | 8    |
| `HAS_DIV`      | `divergence` | `uint16` LE  | 2    |
| `HAS_BITSCORE` | `bit_score`  | `float32` LE | 4    |

`divergence` is an unsigned 16-bit integer divergence measure (`div_pm` in the
reference codec) with implementation-defined scaling; producers and consumers
must agree on the scale out of band.

---

## 5. CIGAR encoding

The CIGAR is the tail of the record body and describes the gap structure of the
alignment as an initial ungapped match run followed by _(indel, match)_ pairs.

```
cigar := uleb(npairs) , uleb(match0) , pair[npairs]
pair  := zigzag_uleb(indel) , uleb(match_after)
```

| Field         | Encoding       | Meaning                                                    |
| ------------- | -------------- | ---------------------------------------------------------- |
| `npairs`      | ULEB128        | Number of _(indel, match)_ pairs that follow.              |
| `match0`      | ULEB128        | Length of the leading ungapped match run (≥ 0).            |
| `indel`       | ZigZag→ULEB128 | Signed gap run (see below). Never 0 in a well-formed pair. |
| `match_after` | ULEB128        | Ungapped match run following this indel (≥ 0).             |

Indel sign convention:

- **`indel > 0`** → `indel` gap characters in the **query** (query has `-`,
  target has letters). Under a _query_-referenced CIGAR string this is an
  **insertion** (`I`).
- **`indel < 0`** → `|indel|` gap characters in the **target** (target has `-`,
  query has letters). Under a _query_-referenced CIGAR string this is a
  **deletion** (`D`).

### 5.1 Derived ungapped span lengths

The ungapped letter counts the CIGAR consumes on each side are:

```
query_letters  = match0 + Σ match_after + Σ (−indel for indel < 0)
target_letters = match0 + Σ match_after + Σ ( indel for indel > 0)
```

These must equal the query/target span lengths (§4.1) for reconstruction to
succeed; the reference decoder raises on mismatch.

### 5.2 Example

`match0 = 5`, pairs `[(+2, 5), (-3, 2)]` describes:

```
query : XXXXX--YYYYYZZZ...    (2 gaps in query)
target: XXXXXYYYYYYY---VV     (3 gaps in target)
```

As a query-referenced CIGAR string: `5M 2I 5M 3D 2M`.

---

## 6. String tables and footer

### 6.1 String tables

Beginning at `tables_start`, three ULEB128-counted string tables appear
**in this fixed order**:

1. `qids` — query identifiers (indexed by `qidx`).
2. `tids` — target identifiers (indexed by `tidx`).
3. `mats` — scoring-system names (indexed by `midx`).

Each table is:

```
table  := uleb(count) , entry[count]
entry  := uleb(byte_len) , utf8_bytes[byte_len]
```

Identifiers are interned: the first time a writer sees a given string it appends
it to the appropriate table and reuses the assigned index thereafter. An empty
file still writes all three tables (each with `count = 0`).

### 6.2 Sequence Letter Tables

If the least significant bit is set in the `extensions` byte in the header (bit `0x01`),
the file contains 2 additional tables after `qids`, `tids`, and `mats`:

4. `qseqs` - query nucleotide sequences (indexed by `qidx`).
5. `tseqs` - target nucleotide sequences (indexed by `tidx`).

Each table is:

```
table  := uleb(count) , entry[count]
entry  := uleb(byte_len) , entry_data[byte_len]
```

Each `entry_data` section contains the following fields, heavily inspired by the
twobitplus format. Unlike the rest of the BPAF format, all indexes are **0-based**.
ULEBs are not used to allow for fast random access while sequence is still on disk
(without needing to load the index fields into memory):

| Field          | Type                                         | Meaning                                                                                                                                                                                  |
| -------------- | -------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `runCount`     | `uint32`                                     | number of runs in this sequence (ascending). Starts with 2-bit, then alternates back and forth between 2 and 4 bit.                                                                      |
| `runEnds`      | `runCount` × `uint32`                        | end of each run (ascending). The next run starts at the current end index. The first run starts at 0.                                                                                    |
| `runDataEnds`  | `runCount` × `uint32`                        | end of the base data for this run (evens are ascending and from 2-bit sequence, odds are ascending 4-bit), Next run starts at the end of the prior run. Both first data runs start at 0. |
| `maskCount`    | `uint32`                                     | number of soft-mask (lower-case) runs.                                                                                                                                                   |
| `maskStarts`   | `maskCount` × `uint32`                       | start of each mask run (ascending).                                                                                                                                                      |
| `maskSizes`    | `maskCount` × `uint32`                       | length of each mask run.                                                                                                                                                                 |
| `twoBitBases`  | ⌈`runDataEnds[lastEvenIndex]` / 4⌉ × `uint8` | 2 bits/base (§6.2.1).                                                                                                                                                                    |
| `fourBitBases` | ⌈`runDataEnds[lastOddIndex]` / 2⌉ × `uint8`  | 4 bits/base (§6.2.2).                                                                                                                                                                    |

The length of the data for a run (`runDataStarts[n] - runDataStarts[n-2]`) can be shorter than the actual run length (`runStarts[n] - runStarts[n-1]`), in which case the seqeunce should be repeated again from the start.
Encoders should only use this feature to encode repeat runs of N's.

#### 6.2.1 Two-bit Encoding

The main sequence is encoded as two bits per base, with the **first base in the most-significant bits** of each byte:

```
A = 00
C = 01
G = 10
T = 11
```

Positions covered by an IUPAC or 4-bit block are excluded. Each 2-bit run is stored contiguously after the other.
Positions covered by a mask-block decode to lower case.

#### 6.2.2 Four-bit IUPAC Encoding

IUPAC bases are stored as four bits per base, the **first base in the most-significant bits** of each byte:

```
A = 0000
C = 0001
G = 0010
T = 0011
B = 0100
D = 0101
H = 0110
K = 0111
M = 1000
N = 1001
R = 1010
S = 1011
V = 1100
W = 1101
X = 1110
Y = 1111
```

This only stores letters within the IUPAC runs. Each run is stored contiguously after
the other in the 4-bit sequence.

This should be applied to the original sequence by **inserting** characters
after each 2-bit run.

The IUPAC data for a run can be shorter than the actual length of the run. If this
happens, then the run data should be repeated until the end of the run is reached.
This technique should only be used for encoding sequences of all N's in IUPAC characters.

### 6.3 Footer (13 bytes)

| Size | Field          | Meaning                                             |
| ---- | -------------- | --------------------------------------------------- |
| 8    | `tables_start` | `uint64` LE — absolute offset of the string tables. |
| 5    | `MAGIC`        | Trailing `"BPAF\x01"` signature.                    |

The trailing `MAGIC` and the `tables_start` pointer let a reader validate the
file and locate the tables via a footer walk-back without scanning records
first.

---

## 7. Versioning and extensions

BPAF distinguishes two independent kinds of change, and they have **different
compatibility contracts**:

### 7.1 `version` — breaking, not backward compatible

The `version` byte identifies the record/table layout. A change to `version`
signals a layout that **older parsers are not expected to read**. Bumping the
version is the mechanism for incompatible changes — reordering fields, changing
an encoding, redefining a field's meaning, etc.

A conforming reader **must reject** any file whose `version` it does not
explicitly implement, rather than attempt a best-effort parse. The reference
reader raises `ValueError: unsupported format version …` in this case. Do not
assume forward compatibility across versions.

### 7.2 `extensions` — additive, backward compatible

The `extensions` byte is a bitfield describing optional additions **within a
given version**. Extensions are expected to be **backward compatible**: a parser
that does not understand an extension bit can still read every field it already
knows about.

To preserve this, extensions add data in ways that do not disturb the existing
layout — most commonly by appending **additional payloads at the end of the
file** (after the footer, or in a designated trailing region), or by using
reserved bits/fields that older readers already skip. The guiding rule for
implementers is:

> **Parse only what you need/expect.** A reader should consume exactly the
> fields defined for the `version` it implements and the `extensions` it
> recognizes, and must tolerate the presence of extension data it does not
> recognize (e.g. trailing bytes) without failing.

Only the first, or least significant bit of the `extensions` field is
reserved in version 0 (for the sequence extension). Readers
should not treat unknown `extensions` bits as fatal.

---

## 8. Reading algorithm (reference)

```
1. Read and verify MAGIC at offset 0.  Reject if absent.
2. Ensure the file is long enough for header + footer.
3. Verify the trailing MAGIC (last 5 bytes).
4. version    = byte[5];  reject if version != 0 (or an unimplemented value).
   extensions = byte[6];  ignore unrecognized bits.
5. tables_start = uint64_le(bytes[-13:-5]);  bounds-check
   (7 <= tables_start <= len - 13).
6. Read the qids, tids, mats tables at tables_start.
8. Read qseqs, tseqs tables if the sequence extension bit in the header (0x01) is set.
7. i = 7  (end of header).
   While i < tables_start:
        reclen, i = uleb(data, i)
        decode the record body in [i, i+reclen)   (or skip it)
        i += reclen
```

---

## 9. Minimal file (empty stream)

The smallest valid BPAF file contains no records and three empty tables:

```
"BPAF\x01"          5   header MAGIC
00                  1   version = 0
00                  1   extensions = 0
00                  1   qids: count = 0
00                  1   tids: count = 0
00                  1   mats: count = 0
07 00 00 00 00 00 00 00   8   tables_start = 7
"BPAF\x01"          5   footer MAGIC
```

Total: 23 bytes; `tables_start = 7` points at the first table-count byte, which
directly follows the 7-byte header.
