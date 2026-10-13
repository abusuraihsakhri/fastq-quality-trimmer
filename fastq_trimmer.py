#!/usr/bin/env python3
"""
FASTQ Quality Trimmer
======================

A small, auditable reimplementation of Trimmomatic's read-trimming logic:

  * SLIDINGWINDOW  - trim once the mean quality within a sliding window
                      drops below a threshold
  * LEADING        - clip low-quality bases from the 5' end
  * TRAILING       - clip low-quality bases from the 3' end
  * ILLUMINACLIP    - a simplified seed-and-extend adapter clip
  * MINLEN         - drop reads that are too short after trimming

Standard library only. See README.md for background and usage.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional, Tuple


# ---------------------------------------------------------------------------
# Quality decoding
# ---------------------------------------------------------------------------

PHRED_OFFSETS = {33: 33, 64: 64}


def decode_quality(qual_str: str, phred: int) -> list[int]:
    """Convert a FASTQ quality string into a list of Phred quality scores."""
    offset = PHRED_OFFSETS[phred]
    if any(ord(c) < offset or ord(c) > 126 for c in qual_str):
        raise ValueError(f"Invalid Phred+{offset} quality character; expected ASCII {offset}–126")
    return [ord(c) - offset for c in qual_str]


# ---------------------------------------------------------------------------
# LEADING / TRAILING clipping
# ---------------------------------------------------------------------------

def leading_trim(quals: list[int], threshold: int) -> int:
    """Return the index of the first base whose quality meets the threshold.

    Every base before that index is clipped from the 5' end. If no base
    meets the threshold, returns len(quals) (the whole read is clipped).
    """
    for i, q in enumerate(quals):
        if q >= threshold:
            return i
    return len(quals)


def trailing_trim(quals: list[int], threshold: int) -> int:
    """Return the exclusive end index that keeps trailing low-quality bases
    clipped from the 3' end. If no base meets the threshold, returns 0.
    """
    for i in range(len(quals) - 1, -1, -1):
        if quals[i] >= threshold:
            return i + 1
    return 0


# ---------------------------------------------------------------------------
# SLIDINGWINDOW trimming
# ---------------------------------------------------------------------------

def sliding_window_trim(quals: list[int], window_size: int, threshold: float) -> int:
    """Scan from the 5' end and return the length to keep.

    Once the mean quality of a window of `window_size` consecutive bases
    drops below `threshold`, the read is cut at the start of that window
    (everything from there to the 3' end is discarded). If the read is
    shorter than the window, its overall mean quality is checked instead.
    """
    n = len(quals)
    if n == 0:
        return 0
    if n <= window_size:
        avg = sum(quals) / n
        return n if avg >= threshold else 0

    window_sum = sum(quals[:window_size])
    if window_sum / window_size < threshold:
        return 0

    for start in range(1, n - window_size + 1):
        window_sum += quals[start + window_size - 1] - quals[start - 1]
        if window_sum / window_size < threshold:
            return start

    return n


# ---------------------------------------------------------------------------
# Simplified ILLUMINACLIP: seed-and-extend adapter clipping
# ---------------------------------------------------------------------------

def find_adapter_clip_pos(
    seq: str,
    adapter: str,
    seed_len: int = 8,
    min_match_len: int = 10,
    max_mismatch_rate: float = 0.1,
) -> Optional[int]:
    """Find the earliest position in `seq` where `adapter` begins to align.

    This is a simplified version of Trimmomatic's palindrome/simple
    ILLUMINACLIP alignment: a short exact "seed" of the adapter is looked
    up at every offset of the read (seed step), and any seed hit is then
    extended base-by-base against the adapter, allowing mismatches
    (extend step). A hit is accepted once the aligned overlap is at least
    `min_match_len` bases and its mismatch rate is at or below
    `max_mismatch_rate` - mirroring the way real adapter read-through only
    needs a partial match near the read's 3' end.

    Returns the read offset where the adapter starts (clip from there to
    the end), or None if no acceptable alignment is found.
    """
    if not adapter:
        return None

    n = len(seq)
    eff_seed_len = min(seed_len, len(adapter))
    seed = adapter[:eff_seed_len]

    for i in range(n):
        overlap = min(n - i, len(adapter))
        if overlap < eff_seed_len:
            # Not enough read left to even match the seed.
            break

        # --- seed step: cheap exact match before the expensive extend ---
        if seq[i:i + eff_seed_len] != seed:
            continue

        # --- extend step: full alignment over the available overlap ---
        mismatches = 0
        for a, b in zip(seq[i:i + overlap], adapter[:overlap]):
            if a != b:
                mismatches += 1

        if overlap >= min_match_len and mismatches / overlap <= max_mismatch_rate:
            return i

    return None


# ---------------------------------------------------------------------------
# Per-read trimming pipeline
# ---------------------------------------------------------------------------

@dataclass
class TrimParams:
    phred: int = 33
    window_size: int = 4
    window_quality: float = 15.0
    leading: int = 3
    trailing: int = 3
    minlen: int = 36
    adapter: Optional[str] = None
    adapter_seed_len: int = 8
    adapter_min_match: int = 10
    adapter_max_mismatch_rate: float = 0.1

    def __post_init__(self):
        """Validate trim parameters after initialization."""
        if self.phred not in (33, 64):
            raise ValueError(f"phred must be 33 or 64, got {self.phred}")
        if self.window_size < 1:
            raise ValueError(f"window_size must be >= 1, got {self.window_size}")
        if not math.isfinite(self.window_quality) or self.window_quality < 0:
            raise ValueError(f"window_quality must be >= 0, got {self.window_quality}")
        if self.leading < 0:
            raise ValueError(f"leading must be >= 0, got {self.leading}")
        if self.trailing < 0:
            raise ValueError(f"trailing must be >= 0, got {self.trailing}")
        if self.minlen < 1:
            raise ValueError(f"minlen must be >= 1, got {self.minlen}")
        if self.adapter_seed_len < 1:
            raise ValueError(f"adapter_seed_len must be >= 1, got {self.adapter_seed_len}")
        if self.adapter_min_match < 1:
            raise ValueError(f"adapter_min_match must be >= 1, got {self.adapter_min_match}")
        if not math.isfinite(self.adapter_max_mismatch_rate) or not (0.0 <= self.adapter_max_mismatch_rate <= 1.0):
            raise ValueError(f"adapter_max_mismatch_rate must be between 0.0 and 1.0, got {self.adapter_max_mismatch_rate}")
        if self.adapter is not None and len(self.adapter) == 0:
            raise ValueError("adapter must be None or a non-empty string")


@dataclass
class TrimStats:
    reads_processed: int = 0
    reads_dropped: int = 0
    total_len_before: int = 0
    total_len_after: int = 0

    @property
    def reads_retained(self) -> int:
        return self.reads_processed - self.reads_dropped

    @property
    def mean_len_before(self) -> float:
        return self.total_len_before / self.reads_processed if self.reads_processed else 0.0

    @property
    def mean_len_after(self) -> float:
        return self.total_len_after / self.reads_retained if self.reads_retained else 0.0


def trim_read(seq: str, qual_str: str, params: TrimParams) -> Optional[Tuple[str, str]]:
    """Run the full trimming pipeline on one read.

    Pipeline order (matches Trimmomatic's recommended step order):
    ILLUMINACLIP -> LEADING -> TRAILING -> SLIDINGWINDOW -> MINLEN.

    Returns (trimmed_seq, trimmed_qual), or None if the read is dropped.
    """
    if len(seq) != len(qual_str):
        raise ValueError("Sequence/quality length mismatch")
    # Validate the complete quality string, including any adapter tail.
    decode_quality(qual_str, params.phred)

    # ILLUMINACLIP (simplified)
    if params.adapter:
        clip_pos = find_adapter_clip_pos(
            seq,
            params.adapter,
            seed_len=params.adapter_seed_len,
            min_match_len=params.adapter_min_match,
            max_mismatch_rate=params.adapter_max_mismatch_rate,
        )
        if clip_pos is not None:
            seq = seq[:clip_pos]
            qual_str = qual_str[:clip_pos]

    if not seq:
        return None

    quals = decode_quality(qual_str, params.phred)

    # LEADING
    start = leading_trim(quals, params.leading)
    # TRAILING (computed on the full read; combined below with `start`)
    end = trailing_trim(quals, params.trailing)

    if start >= end:
        return None

    sub_quals = quals[start:end]

    # SLIDINGWINDOW
    keep_len = sliding_window_trim(sub_quals, params.window_size, params.window_quality)
    if keep_len == 0:
        return None

    final_end = start + keep_len

    trimmed_seq = seq[start:final_end]
    trimmed_qual = qual_str[start:final_end]

    # MINLEN
    if len(trimmed_seq) < params.minlen:
        return None

    return trimmed_seq, trimmed_qual


# ---------------------------------------------------------------------------
# FASTQ I/O
# ---------------------------------------------------------------------------

@dataclass
class FastqRecord:
    header: str
    seq: str
    plus: str
    qual: str


def parse_fastq(lines: Iterator[str]) -> Iterator[FastqRecord]:
    """Yield FastqRecord objects from an iterable of lines (4 lines/record)."""
    while True:
        try:
            header = next(lines).rstrip("\r\n")
        except StopIteration:
            return
        if header == "":
            continue
        try:
            seq = next(lines).rstrip("\r\n")
            plus = next(lines).rstrip("\r\n")
            qual = next(lines).rstrip("\r\n")
        except StopIteration:
            raise ValueError("Truncated FASTQ record at end of file") from None

        if not header.startswith("@"):
            raise ValueError(f"Malformed FASTQ header (expected '@'): {header!r}")
        if not plus.startswith("+"):
            raise ValueError(f"Malformed FASTQ separator (expected '+'): {plus!r}")
        if len(seq) != len(qual):
            raise ValueError(
                f"Sequence/quality length mismatch for read {header!r}: "
                f"{len(seq)} vs {len(qual)}"
            )

        yield FastqRecord(header, seq, plus, qual)


def process_fastq(
    in_lines: Iterator[str], out_file, params: TrimParams
) -> TrimStats:
    stats = TrimStats()
    for record in parse_fastq(in_lines):
        stats.reads_processed += 1
        stats.total_len_before += len(record.seq)

        result = trim_read(record.seq, record.qual, params)
        if result is None:
            stats.reads_dropped += 1
            continue

        trimmed_seq, trimmed_qual = result
        stats.total_len_after += len(trimmed_seq)
        out_file.write(f"{record.header}\n{trimmed_seq}\n{record.plus}\n{trimmed_qual}\n")

    return stats


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fastq_trimmer",
        description=(
            "Reimplementation of Trimmomatic-style sliding-window quality "
            "trimming and adapter clipping for FASTQ reads."
        ),
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file")
    parser.add_argument("-o", "--output", required=True, help="Output (trimmed) FASTQ file")
    parser.add_argument(
        "--phred", type=int, choices=(33, 64), default=33,
        help="Quality encoding of the input file (default: 33)",
    )
    parser.add_argument(
        "--window-size", type=int, default=4,
        help="SLIDINGWINDOW: number of bases averaged per window (default: 4)",
    )
    parser.add_argument(
        "--window-quality", type=float, default=15.0,
        help="SLIDINGWINDOW: minimum mean quality required per window (default: 15)",
    )
    parser.add_argument(
        "--leading", type=int, default=3,
        help="LEADING: minimum quality to keep a base at the read start (default: 3)",
    )
    parser.add_argument(
        "--trailing", type=int, default=3,
        help="TRAILING: minimum quality to keep a base at the read end (default: 3)",
    )
    parser.add_argument(
        "--minlen", type=int, default=36,
        help="MINLEN: minimum read length to keep after trimming (default: 36)",
    )
    parser.add_argument(
        "--adapter", default=None,
        help="Adapter sequence to clip (ILLUMINACLIP-style). Omit to skip adapter clipping.",
    )
    parser.add_argument(
        "--adapter-seed-len", type=int, default=8,
        help="Adapter clip: exact-match seed length (default: 8)",
    )
    parser.add_argument(
        "--adapter-min-match", type=int, default=10,
        help="Adapter clip: minimum aligned overlap to accept a match (default: 10)",
    )
    parser.add_argument(
        "--adapter-max-mismatch-rate", type=float, default=0.1,
        help="Adapter clip: maximum fraction of mismatches tolerated in the aligned overlap (default: 0.1)",
    )
    return parser


def params_from_args(args: argparse.Namespace) -> TrimParams:
    return TrimParams(
        phred=args.phred,
        window_size=args.window_size,
        window_quality=args.window_quality,
        leading=args.leading,
        trailing=args.trailing,
        minlen=args.minlen,
        adapter=args.adapter,
        adapter_seed_len=args.adapter_seed_len,
        adapter_min_match=args.adapter_min_match,
        adapter_max_mismatch_rate=args.adapter_max_mismatch_rate,
    )


def print_summary(stats: TrimStats) -> None:
    drop_pct = (stats.reads_dropped / stats.reads_processed * 100) if stats.reads_processed else 0.0
    print("=== FASTQ Quality Trimmer Summary ===")
    print(f"Reads processed:    {stats.reads_processed}")
    print(f"Reads dropped:      {stats.reads_dropped} ({drop_pct:.2f}%)")
    print(f"Reads retained:     {stats.reads_retained}")
    print(f"Mean length before: {stats.mean_len_before:.2f} bp")
    print(f"Mean length after:  {stats.mean_len_after:.2f} bp")


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    params = params_from_args(args)

    # Validate input file exists
    input_path = Path(args.input)
    if not input_path.is_file():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        return 1

    # Never truncate the input in place.
    output_path = Path(args.output)
    if input_path.resolve() == output_path.resolve():
        print("Error: Input and output must be different files", file=sys.stderr)
        return 1
    try:
        # Ensure parent directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        print(f"Error: Cannot create output directory: {e}", file=sys.stderr)
        return 1

    # Write to a temporary sibling; malformed inputs never destroy existing outputs.
    temp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n",
                                         dir=output_path.parent, prefix=".fastq-",
                                         suffix=".tmp", delete=False) as out_f:
            temp_path = Path(out_f.name)
            with open(input_path, "r", encoding="utf-8") as in_f:
                stats = process_fastq(iter(in_f), out_f, params)
        os.replace(temp_path, output_path)
    except (OSError, UnicodeError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)

    print_summary(stats)
    return 0


if __name__ == "__main__":
    sys.exit(main())
