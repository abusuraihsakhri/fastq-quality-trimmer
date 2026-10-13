# FASTQ Quality Trimmer

### [Open the Live Application →](https://abusuraihsakhri.github.io/fastq-quality-trimmer/)

A small single-end FASTQ quality and adapter trimmer with a standard-library Python CLI and a dependency-free browser interface. The algorithms are intentionally simple and auditable; this is **not** the full Trimmomatic implementation.

## Features

- **ILLUMINACLIP (simplified):** exact adapter-prefix seed followed by mismatch-tolerant extension.
- **LEADING / TRAILING:** remove bases below a minimum Phred score at each end.
- **SLIDINGWINDOW:** cut at the first window whose mean quality falls below the threshold.
- **MINLEN:** discard short trimmed reads.
- Preserve read headers and FASTQ separator annotations; reject truncated records, invalid quality encodings, and sequence/quality length mismatches.
- Summarize processed/retained/dropped reads and sequence lengths. Write Python CLI output atomically to avoid destroying valid output when an input fails validation.

## Browser application

The `web/` directory is a static, responsive FASTQ trimmer. It runs entirely in your browser with **JavaScript**, not server-side Python, Pyodide, or a remote API. Reads stay on your device; the app makes no external network requests. Select an uncompressed `.fastq` / `.fq` file (up to 20 MB), or paste FASTQ content, adjust the parameters, run trimming, and download a trimmed FASTQ file. A three-read example is provided.

For local use:

```bash
python -m http.server 8080 --directory web
# Open http://localhost:8080
```

Browsers supporting standard ES modules, File APIs and Blob downloads are required. Compressed `.fastq.gz` input, paired-end synchronization and streaming multi-gigabyte datasets are not implemented. The browser version uses JavaScript that mirrors the Python single-end algorithm, covered by Node tests.

## Python installation

Python 3.10 or later. Core FASTQ trimming uses only the Python standard library.

```bash
git clone https://github.com/abusuraihsakhri/fastq-quality-trimmer.git
cd fastq-quality-trimmer
python -m pip install -e .
```

The install brings in FastAPI, Uvicorn and Pydantic for optional auxiliary supervisor commands. Direct use of `python fastq_trimmer.py` does not require these packages.

## Command-line usage

```bash
python fastq_trimmer.py -i reads.fastq -o reads.trimmed.fastq \
  --phred 33 --window-size 4 --window-quality 15 \
  --leading 3 --trailing 3 --minlen 36 \
  --adapter AGATCGGAAGAGC

# Or, after installation:
fastq-trimmer -i reads.fastq -o reads.trimmed.fastq --minlen 36
```

Defaults: Phred+33, window of 4 bases, window quality 15, leading/trailing quality 3, minimum length 36. Adapter clipping is disabled unless `--adapter` is specified. Optional adapter settings: `--adapter-seed-len 8`, `--adapter-min-match 10`, and `--adapter-max-mismatch-rate 0.1`. Input and output must be different files. The input must be uncompressed FASTQ text.

### Python API

```python
from fastq_trimmer import TrimParams, process_fastq

params = TrimParams(minlen=30, adapter="AGATCGGAAGAGC")
with open("reads.fastq", encoding="utf-8") as src, \
     open("trimmed.fastq", "w", encoding="utf-8") as dest:
    stats = process_fastq(iter(src), dest, params)
print(stats.reads_processed, stats.reads_retained, stats.mean_len_after)
```

Unlike the CLI, `process_fastq()` writes directly to its supplied output stream; callers are responsible for output transaction safety.

## Optional supervisor and REST API

`cli.py`, `agents/`, `enrichment.py` and `simulator.py` provide a separate **mock demonstration supervisor** with rule-based alerts, synthetic metrics and an in-memory HMAC audit trail. They do **not** trim reads and are not clinical decision support or validated CAP/CLSI compliance tools. The mock model backend never performs genuine LLM inference.

```bash
python cli.py serve --host 127.0.0.1 --port 8000
# GET http://127.0.0.1:8000/health
```

To use Docker Compose, supply a strong `AUDIT_SECRET_KEY` environment value first. Compose binds port 8000 on localhost by default. The separate API has no user authentication and should not be published directly to the internet. The identifier-pattern filter is **not** sufficient to guarantee PHI de-identification; do not send clinical or sensitive personal records to it. Ledger integrity checking verifies HMAC signatures within the running process, but the ledger is not persisted between restarts.

## Testing

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
node --test web/test_trimmer.mjs   # Node.js 22
```

GitHub Actions executes the Python test suite on 3.10–3.13 and the JavaScript engine tests. On `master`, a successful run also attempts publication of the `web/` app to GitHub Pages.

## License

MIT — see [LICENSE](LICENSE).
