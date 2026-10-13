import test from "node:test";
import assert from "node:assert/strict";
import { parseFastq, trimFastq, trimRead, validateOptions } from "./trimmer.mjs";

const adapter = "AGATCGGAAGAGC";
const record = (id, seq, qual, plus = "+") => [`@${id}`, seq, plus, qual].join("\n") + "\n";

test("retains a high-quality read with original separator", () => {
  const seq = "ACGT".repeat(10);
  const input = record("read1", seq, "I".repeat(seq.length), "+read1");
  const result = trimFastq(input);
  assert.equal(result.fastq, input);
  assert.deepEqual(result.stats, {
    processed: 1, retained: 1, dropped: 0,
    basesBefore: 40, basesAfter: 40
  });
});

test("clips adapters and drops low-quality reads", () => {
  const seq = "ACGT".repeat(10);
  const input = record("adapter", seq + adapter, "I".repeat(53)) +
                record("low", seq, "#".repeat(40));
  const result = trimFastq(input, { adapter });
  assert.equal(result.stats.processed, 2);
  assert.equal(result.stats.retained, 1);
  assert.equal(result.stats.dropped, 1);
  assert.equal(result.fastq, record("adapter", seq, "I".repeat(40)));
});

test("sliding-window and short-read semantics match Python", () => {
  const seq = "A".repeat(30) + "C".repeat(20);
  assert.equal(trimRead(seq, "I".repeat(30) + "#".repeat(20),
    { minLength: 25 })?.seq, "A".repeat(30));
  assert.equal(trimRead(seq, "I".repeat(30) + "#".repeat(20),
    { minLength: 35 }), null);
  assert.equal(trimRead("AC", "II", { minLength: 1, windowSize: 4 })?.seq, "AC");
  assert.equal(trimRead("AC", "##", { minLength: 1, windowSize: 4 }), null);
});

test("CRLF and missing final newline are accepted", () => {
  const sample = record("windows", "ACGT", "IIII", "+windows");
  const records = Array.from(parseFastq(sample.replace(/\n/g, "\r\n").trimEnd()));
  assert.equal(records.length, 1);
  assert.equal(trimFastq(sample.replace(/\n/g, "\r\n").trimEnd(),
    { minLength: 1 }).stats.retained, 1);
});

test("malformed FASTQ and invalid quality are rejected", () => {
  assert.throws(() => trimFastq("@only\nA\n+"), /Incomplete/);
  assert.throws(() => trimFastq(record("bad", "AAA", "II")), /length mismatch/);
  assert.throws(() => trimFastq(record("bad", "AA", "  "), { minLength: 1 }), /Invalid Phred/);
  assert.throws(() => trimFastq(record("bad", "AA", "##"), { phred: 64, minLength: 1 }), /Invalid Phred/);
  assert.throws(() => trimFastq("not-fastq\nAA\n+\nII"), /header/);
});

test("numeric parameters reject non-finite and fractional invalid inputs", () => {
  assert.throws(() => validateOptions({ windowSize: 1.5 }), /windowSize/);
  assert.throws(() => validateOptions({ windowQuality: NaN }), /Window quality/);
  assert.throws(() => validateOptions({ adapterMismatchRate: Infinity }), /mismatch rate/);
  assert.throws(() => validateOptions({ phred: 42 }), /Phred/);
});

test("adapter seed permits a single tail mismatch within rate", () => {
  const mutant = adapter.slice(0, -1) + "T";
  assert.equal(trimRead("A".repeat(40) + mutant, "I".repeat(53),
    { adapter, minLength: 30 })?.seq, "A".repeat(40));
});
