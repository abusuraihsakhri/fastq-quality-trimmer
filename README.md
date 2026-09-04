# Fastq Quality Trimmer

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)

A small, auditable reimplementation of Trimmomatic's read-trimming logic for FASTQ files.

## What It Does

FASTQ Quality Trimmer provides the following trimming operations:

- **SLIDINGWINDOW** - Trim once the mean quality within a sliding window drops below a threshold
- **LEADING** - Clip low-quality bases from the 5' end
- **TRAILING** - Clip low-quality bases from the 3' end
- **ILLUMINACLIP** - A simplified seed-and-extend adapter clip
- **MINLEN** - Drop reads that are too short after trimming

## Installation

```bash
# Clone the repository
git clone https://github.com/abusuraihsakhri/fastq-quality-trimmer.git
cd fastq-quality-trimmer

# Install dependencies
pip install -e .

# Or install with development dependencies
pip install -e ".[dev]"
```

## Usage

### Command Line Interface

Trim a FASTQ file:

```bash
python fastq_trimmer.py -i input.fastq -o output.fastq
```

With custom parameters:

```bash
python fastq_trimmer.py -i input.fastq -o output.fastq \
    --phred 33 \
    --window-size 4 \
    --window-quality 15 \
    --leading 3 \
    --trailing 3 \
    --minlen 36 \
    --adapter AGATCGGAAGAGC
```

### CLI Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-i, --input` | (required) | Input FASTQ file |
| `-o, --output` | (required) | Output (trimmed) FASTQ file |
| `--phred` | 33 | Quality encoding (33 or 64) |
| `--window-size` | 4 | SLIDINGWINDOW: bases per window |
| `--window-quality` | 15.0 | SLIDINGWINDOW: minimum mean quality |
| `--leading` | 3 | LEADING: minimum quality at read start |
| `--trailing` | 3 | TRAILING: minimum quality at read end |
| `--minlen` | 36 | MINLEN: minimum read length to keep |
| `--adapter` | None | Adapter sequence to clip |
| `--adapter-seed-len` | 8 | Adapter clip: exact-match seed length |
| `--adapter-min-match` | 10 | Adapter clip: minimum aligned overlap |
| `--adapter-max-mismatch-rate` | 0.1 | Adapter clip: max mismatch fraction |

### Python API

```python
from fastq_trimmer import TrimParams, process_fastq

params = TrimParams(
    phred=33,
    window_size=4,
    window_quality=15.0,
    leading=3,
    trailing=3,
    minlen=36,
    adapter="AGATCGGAAGAGC"
)

with open("input.fastq") as in_f, open("output.fastq", "w") as out_f:
    stats = process_fastq(iter(in_f), out_f, params)

print(f"Processed: {stats.reads_processed}, Retained: {stats.reads_retained}")
```

## Testing

Run the full test suite:

```bash
pytest -v
```

Run specific test files:

```bash
pytest test_fastq_trimmer.py -v
pytest tests/ -v
```

## Security

### PHI Guard
The system includes a Zero-PHI Outbound Interceptor that blocks SSNs, MRNs, phone numbers, and patient identifiers from being logged or output.

### Audit Trail
HMAC-SHA256 tamper-evident audit trail for all operations. Set the `AUDIT_SECRET_KEY` environment variable for production use:

```bash
# Generate a secure key
python -c "import os; print(os.urandom(32).hex())"

# Set it
export AUDIT_SECRET_KEY="your-secure-key-hello"
```

See `.env.example` for configuration options.

## Docker Deployment

```bash
# Build and run with Docker Compose
AUDIT_SECRET_KEY=$(python -c "import os; print(os.urandom(32).hex())") docker-compose up --build

# Or build and run manually
docker build -t fastq-quality-trimmer .
docker run -p 8000:8000 -e AUDIT_SECRET_KEY="your-key" fastq-quality-trimmer
```

## Project Structure

```
fastq-quality-trimmer/
├── fastq_trimmer.py      # Core trimming algorithms and CLI
├── cli.py                # Enterprise supervisor CLI
├── enrichment.py         # Enrichment feature engines
├── simulator.py          # High-throughput simulation
├── agents/               # Multi-agent supervisor system
│   ├── base.py           # Security, PHI guard, audit trail
│   ├── models.py         # Pydantic data models
│   ├── supervisor.py     # Supervisor orchestrator
│   ├── workers.py        # Specialized worker agents
│   ├── api.py            # FastAPI REST server
│   ├── metrics.py        # Prometheus metrics
│   ├── learning.py       # Bayesian calibration engine
│   └── streamer.py       # WebSocket telemetry
├── tests/                # Pytest test suite
├── web/                  # Operations console (HTML)
├── pyproject.toml        # Project metadata and dependencies
├── Dockerfile            # Container build
├── docker-compose.yml    # Container orchestration
└── .github/workflows/    # CI/CD pipeline
```

## License

MIT License - see [LICENSE](LICENSE) for details.
