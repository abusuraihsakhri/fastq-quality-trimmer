"""Regression tests for FASTQ integrity, validation, and ledger tampering."""
import io

import pytest

from agents.base import AuditTrail
from cli import parse_boolean
from fastq_trimmer import TrimParams, decode_quality, main, process_fastq, trim_read


def test_crlf_input_and_separator_metadata_preserved():
    source = "@r1\r\nACGT\r\n+r1\r\nIIII\r\n"
    out = io.StringIO()
    stats = process_fastq(iter(io.StringIO(source)), out, TrimParams(minlen=1))
    assert stats.reads_retained == 1
    assert out.getvalue() == "@r1\nACGT\n+r1\nIIII\n"


def test_invalid_quality_characters_are_rejected():
    with pytest.raises(ValueError, match="Invalid Phred"):
        decode_quality(" ", 33)
    with pytest.raises(ValueError, match="Invalid Phred"):
        decode_quality("#", 64)
    with pytest.raises(ValueError, match="Invalid Phred"):
        trim_read("A" * 40, "I" * 39 + " ", TrimParams(minlen=1))


@pytest.mark.parametrize("params", [
    {"window_quality": float("nan")},
    {"window_quality": float("inf")},
    {"adapter_max_mismatch_rate": float("nan")},
])
def test_nonfinite_parameters_rejected(params):
    with pytest.raises(ValueError):
        TrimParams(**params)


def test_input_cannot_be_overwritten(tmp_path):
    source = tmp_path / "reads.fastq"
    original = "@read\nACGT\n+\nIIII\n"
    source.write_text(original)
    assert main(["-i", str(source), "-o", str(source), "--minlen", "1"]) == 1
    assert source.read_text() == original


def test_malformed_input_does_not_overwrite_output(tmp_path):
    source = tmp_path / "reads.fastq"
    target = tmp_path / "trimmed.fastq"
    source.write_text("@broken\nACGT\n+\n")
    target.write_text("existing output")
    assert main(["-i", str(source), "-o", str(target), "--minlen", "1"]) == 1
    assert target.read_text() == "existing output"
    assert not list(tmp_path.glob(".fastq-*.tmp"))


def test_hmac_verification_detects_tampering_and_snapshot_is_immutable():
    ledger = AuditTrail(secret_key="unit-test-secret")
    ledger.log("tester", "worker", "TASK", {"status": "ok"})
    assert ledger.verify_integrity()
    snapshot = ledger.get_trail()
    snapshot[0]["actor"] = "changed via response"
    assert ledger.verify_integrity()
    ledger.logs[0]["actor"] = "forged"
    assert not ledger.verify_integrity()


def test_csv_boolean_parser():
    assert parse_boolean("false") is False
    assert parse_boolean("0") is False
    assert parse_boolean("TRUE") is True
    assert parse_boolean("yes") is True
    with pytest.raises(ValueError):
        parse_boolean("perhaps")
