/**
 * Local FASTQ trimmer. The order and default parameters mirror fastq_trimmer.py.
 * This is a simplified single-end, seed-and-extend implementation, not Trimmomatic.
 * No network access, external dependencies, or server-side processing.
 */
export const defaults = Object.freeze({
  phred: 33, windowSize: 4, windowQuality: 15, leading: 3,
  trailing: 3, minLength: 36, adapter: "",
  adapterSeedLength: 8, adapterMinMatch: 10, adapterMismatchRate: 0.1
});

export function validateOptions(options = {}) {
  const p = { ...defaults, ...options };
  for (const [name, minimum] of [
    ["windowSize", 1], ["leading", 0], ["trailing", 0],
    ["minLength", 1], ["adapterSeedLength", 1], ["adapterMinMatch", 1]
  ]) {
    if (!Number.isSafeInteger(p[name]) || p[name] < minimum) {
      throw new Error(`${name} must be an integer >= ${minimum}`);
    }
  }
  if (p.phred !== 33 && p.phred !== 64) {
    throw new Error("Phred offset must be 33 or 64");
  }
  if (!Number.isFinite(p.windowQuality) || p.windowQuality < 0) {
    throw new Error("Window quality must be a finite nonnegative number");
  }
  if (!Number.isFinite(p.adapterMismatchRate) ||
      p.adapterMismatchRate < 0 || p.adapterMismatchRate > 1) {
    throw new Error("Adapter mismatch rate must be between 0 and 1");
  }
  if (typeof p.adapter !== "string") {
    throw new Error("Adapter must be a string");
  }
  return p;
}

export function* parseFastq(text) {
  if (typeof text !== "string") throw new Error("FASTQ input must be text");
  const lines = text.replace(/\r\n?/g, "\n").split("\n");
  let i = 0;
  while (i < lines.length) {
    if (lines[i] === "") { i++; continue; }
    const recordNumber = Math.floor(i / 4) + 1;
    if (i + 3 >= lines.length) {
      throw new Error(`Incomplete FASTQ record near line ${i + 1}`);
    }
    const [header, seq, plus, qual] = lines.slice(i, i + 4);
    if (!header.startsWith("@")) {
      throw new Error(`Invalid FASTQ header at line ${i + 1}`);
    }
    if (!plus.startsWith("+")) {
      throw new Error(`Invalid FASTQ separator at line ${i + 3}`);
    }
    if (seq.length !== qual.length) {
      throw new Error(`Sequence/quality length mismatch at line ${i + 1}`);
    }
    yield { header, seq, plus, qual };
    i += 4;
  }
}

function clipAdapter(seq, adapter, seedLength, minMatch, mismatchRate) {
  if (!adapter) return null;
  const seed = adapter.slice(0, Math.min(seedLength, adapter.length));
  for (let i = 0; i < seq.length; i++) {
    const overlap = Math.min(seq.length - i, adapter.length);
    if (overlap < seed.length) break;
    if (seq.slice(i, i + seed.length) !== seed) continue;
    let mismatches = 0;
    for (let j = 0; j < overlap; j++) {
      if (seq[i + j] !== adapter[j]) mismatches++;
    }
    if (overlap >= minMatch && mismatches / overlap <= mismatchRate) return i;
  }
  return null;
}

function clipWindow(qualities, size, threshold) {
  const n = qualities.length;
  if (n === 0) return 0;
  if (n <= size) {
    return qualities.reduce((sum, q) => sum + q, 0) / n >= threshold ? n : 0;
  }
  let sum = 0;
  for (let i = 0; i < size; i++) sum += qualities[i];
  if (sum / size < threshold) return 0;
  for (let start = 1; start <= n - size; start++) {
    sum += qualities[start + size - 1] - qualities[start - 1];
    if (sum / size < threshold) return start;
  }
  return n;
}

export function trimRead(seq, qual, options = {}) {
  const p = validateOptions(options);
  if (seq.length !== qual.length) throw new Error("Sequence/quality length mismatch");
  // Always validate the entire input, even if an adapter tail will be discarded.
  for (const char of qual) {
    const value = char.charCodeAt(0);
    if (value < p.phred || value > 126) {
      throw new Error(`Invalid Phred+${p.phred} quality character; expected ASCII ${p.phred}–126`);
    }
  }
  const clip = clipAdapter(seq, p.adapter, p.adapterSeedLength,
                           p.adapterMinMatch, p.adapterMismatchRate);
  if (clip !== null) {
    seq = seq.slice(0, clip);
    qual = qual.slice(0, clip);
  }
  if (!seq.length) return null;
  const qualities = Array.from(qual, c => c.charCodeAt(0) - p.phred);
  let start = 0;
  while (start < qualities.length && qualities[start] < p.leading) start++;
  let end = qualities.length;
  while (end > 0 && qualities[end - 1] < p.trailing) end--;
  if (start >= end) return null;
  const kept = clipWindow(qualities.slice(start, end),
                          p.windowSize, p.windowQuality);
  if (kept === 0 || kept < p.minLength) return null;
  return { seq: seq.slice(start, start + kept),
           qual: qual.slice(start, start + kept) };
}

export function trimFastq(text, options = {}) {
  const p = validateOptions(options);
  const lines = [];
  const stats = { processed: 0, retained: 0, dropped: 0, basesBefore: 0, basesAfter: 0 };
  for (const record of parseFastq(text)) {
    stats.processed++;
    stats.basesBefore += record.seq.length;
    const result = trimRead(record.seq, record.qual, p);
    if (result === null) {
      stats.dropped++;
      continue;
    }
    stats.retained++;
    stats.basesAfter += result.seq.length;
    lines.push(record.header, result.seq, record.plus, result.qual);
  }
  return { fastq: lines.length ? lines.join("\n") + "\n" : "", stats };
}
