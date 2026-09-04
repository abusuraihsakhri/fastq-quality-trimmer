"""Tests for fastq_trimmer.py — run with: python test_fastq_trimmer.py"""

import io
import unittest

from fastq_trimmer import (
    TrimParams,
    decode_quality,
    find_adapter_clip_pos,
    leading_trim,
    process_fastq,
    sliding_window_trim,
    trailing_trim,
    trim_read,
)

ADAPTER = "AGATCGGAAGAGC"  # Illumina universal adapter (13bp)
HIGH_Q = "I"  # Phred33 'I' -> Q40
LOW_Q = "#"   # Phred33 '#' -> Q2


class TestQualityDecoding(unittest.TestCase):
    def test_phred33(self):
        self.assertEqual(decode_quality("I#", 33), [40, 2])

    def test_phred64(self):
        # 'h' -> 104 - 64 = 40 ; '@' -> 64 - 64 = 0
        self.assertEqual(decode_quality("h@", 64), [40, 0])


class TestLeadingTrailing(unittest.TestCase):
    def test_leading_trim_clips_low_quality_start(self):
        quals = [2, 2, 25, 25, 25]
        self.assertEqual(leading_trim(quals, 3), 2)

    def test_trailing_trim_clips_low_quality_end(self):
        quals = [25, 25, 25, 2, 2]
        self.assertEqual(trailing_trim(quals, 3), 3)

    def test_leading_trim_all_low_quality(self):
        quals = [1, 1, 1]
        self.assertEqual(leading_trim(quals, 3), 3)

    def test_trailing_trim_all_low_quality(self):
        quals = [1, 1, 1]
        self.assertEqual(trailing_trim(quals, 3), 0)


class TestSlidingWindowTrim(unittest.TestCase):
    def test_cuts_at_start_of_failing_window(self):
        # window=4, threshold=20
        # window avgs: [30,30,30,30]=30 [30,30,30,10]=25 [30,30,10,10]=20 [30,10,10,10]=15<20 -> cut at 3
        quals = [30, 30, 30, 30, 10, 10, 10, 10, 30, 30]
        self.assertEqual(sliding_window_trim(quals, 4, 20), 3)

    def test_keeps_whole_read_when_quality_never_drops(self):
        quals = [30] * 20
        self.assertEqual(sliding_window_trim(quals, 4, 20), 20)

    def test_short_read_uses_overall_mean(self):
        # shorter than the window: whole-read mean is checked
        quals = [10, 10]
        self.assertEqual(sliding_window_trim(quals, 4, 20), 0)
        quals_ok = [25, 25]
        self.assertEqual(sliding_window_trim(quals_ok, 4, 20), 2)


class TestAdapterClip(unittest.TestCase):
    def test_exact_adapter_match_found(self):
        seq = "ACGTACGTAC" + ADAPTER  # adapter starts at index 10
        self.assertEqual(find_adapter_clip_pos(seq, ADAPTER), 10)

    def test_partial_match_within_mismatch_tolerance(self):
        # last base of the adapter mutated -> 1/13 mismatches (~7.7%), seed (first 8bp) untouched
        mutated_tail = ADAPTER[:-1] + ("T" if ADAPTER[-1] != "T" else "A")
        seq = "ACGTACGTAC" + mutated_tail
        self.assertEqual(find_adapter_clip_pos(seq, ADAPTER, max_mismatch_rate=0.1), 10)

    def test_no_match_when_adapter_absent(self):
        seq = "ACGTACGTACGTACGTACGTACGT"
        self.assertIsNone(find_adapter_clip_pos(seq, ADAPTER))


class TestTrimReadPipeline(unittest.TestCase):
    def test_high_quality_read_passes_unchanged(self):
        seq = "ACGT" * 12 + "AC"  # 50bp
        qual = HIGH_Q * len(seq)
        params = TrimParams(minlen=36)
        result = trim_read(seq, qual, params)
        self.assertEqual(result, (seq, qual))

    def test_adapter_clipped_then_kept_above_minlen(self):
        genomic = "ACGT" * 10  # 40bp
        seq = genomic + ADAPTER  # 53bp
        qual = HIGH_Q * len(seq)
        params = TrimParams(minlen=30, adapter=ADAPTER)
        result = trim_read(seq, qual, params)
        self.assertIsNotNone(result)
        trimmed_seq, trimmed_qual = result
        self.assertEqual(trimmed_seq, genomic)
        self.assertEqual(trimmed_qual, HIGH_Q * 40)

    def test_low_quality_read_is_dropped(self):
        seq = "ACGT" * 10  # 40bp
        qual = LOW_Q * len(seq)
        params = TrimParams(minlen=30)
        self.assertIsNone(trim_read(seq, qual, params))

    def test_sliding_window_trim_then_minlen_drop(self):
        # 30bp high quality followed by 20bp low quality; sliding window
        # should cut the tail, leaving 30bp, which is below minlen=35.
        seq = "A" * 30 + "C" * 20
        qual = HIGH_Q * 30 + LOW_Q * 20
        params = TrimParams(minlen=35, window_size=4, window_quality=15)
        self.assertIsNone(trim_read(seq, qual, params))

    def test_sliding_window_trim_then_minlen_kept(self):
        seq = "A" * 30 + "C" * 20
        qual = HIGH_Q * 30 + LOW_Q * 20
        params = TrimParams(minlen=25, window_size=4, window_quality=15)
        result = trim_read(seq, qual, params)
        self.assertIsNotNone(result)
        trimmed_seq, _ = result
        self.assertEqual(trimmed_seq, "A" * 30)


class TestProcessFastqEndToEnd(unittest.TestCase):
    def test_summary_stats_and_output_content(self):
        genomic = "ACGT" * 10  # 40bp
        read1_seq = "ACGT" * 12 + "AC"  # 50bp, high quality, kept as-is
        read2_seq = genomic + ADAPTER  # 53bp, adapter clipped to 40bp
        read3_seq = genomic  # 40bp, all low quality -> dropped

        fastq_text = (
            f"@read1\n{read1_seq}\n+\n{HIGH_Q * len(read1_seq)}\n"
            f"@read2\n{read2_seq}\n+\n{HIGH_Q * len(read2_seq)}\n"
            f"@read3\n{read3_seq}\n+\n{LOW_Q * len(read3_seq)}\n"
        )

        params = TrimParams(minlen=30, adapter=ADAPTER)
        out = io.StringIO()
        stats = process_fastq(iter(io.StringIO(fastq_text)), out, params)

        self.assertEqual(stats.reads_processed, 3)
        self.assertEqual(stats.reads_dropped, 1)
        self.assertEqual(stats.reads_retained, 2)
        self.assertAlmostEqual(stats.total_len_before, 143)
        self.assertAlmostEqual(stats.mean_len_before, 143 / 3)
        self.assertAlmostEqual(stats.total_len_after, 90)
        self.assertAlmostEqual(stats.mean_len_after, 45.0)

        output_lines = out.getvalue().splitlines()
        self.assertEqual(output_lines[0], "@read1")
        self.assertEqual(output_lines[1], read1_seq)
        self.assertEqual(output_lines[4], "@read2")
        self.assertEqual(output_lines[5], genomic)
        # read3 must not appear at all
        self.assertNotIn("@read3", output_lines)


class TestTrimParamsValidation(unittest.TestCase):
    def test_default_params_valid(self):
        params = TrimParams()
        self.assertEqual(params.phred, 33)
        self.assertEqual(params.window_size, 4)

    def test_invalid_phred_raises(self):
        with self.assertRaises(ValueError):
            TrimParams(phred=42)

    def test_invalid_window_size_raises(self):
        with self.assertRaises(ValueError):
            TrimParams(window_size=0)

    def test_negative_window_quality_raises(self):
        with self.assertRaises(ValueError):
            TrimParams(window_quality=-1)

    def test_negative_leading_raises(self):
        with self.assertRaises(ValueError):
            TrimParams(leading=-1)

    def test_negative_trailing_raises(self):
        with self.assertRaises(ValueError):
            TrimParams(trailing=-1)

    def test_invalid_minlen_raises(self):
        with self.assertRaises(ValueError):
            TrimParams(minlen=0)

    def test_invalid_adapter_seed_len_raises(self):
        with self.assertRaises(ValueError):
            TrimParams(adapter_seed_len=0)

    def test_invalid_adapter_min_match_raises(self):
        with self.assertRaises(ValueError):
            TrimParams(adapter_min_match=0)

    def test_invalid_mismatch_rate_raises(self):
        with self.assertRaises(ValueError):
            TrimParams(adapter_max_mismatch_rate=1.5)

    def test_negative_mismatch_rate_raises(self):
        with self.assertRaises(ValueError):
            TrimParams(adapter_max_mismatch_rate=-0.1)

    def test_empty_adapter_raises(self):
        with self.assertRaises(ValueError):
            TrimParams(adapter="")

    def test_valid_adapter_accepted(self):
        params = TrimParams(adapter="AGATCGGAAGAGC")
        self.assertEqual(params.adapter, "AGATCGGAAGAGC")


class TestCLIFileValidation(unittest.TestCase):
    def test_missing_input_file_returns_error(self):
        from fastq_trimmer import main
        result = main(["-i", "/nonexistent/path.fastq", "-o", "/tmp/out.fastq"])
        self.assertEqual(result, 1)


if __name__ == "__main__":
    unittest.main()
