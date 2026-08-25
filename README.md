# Fastq Quality Trimmer

> **Clinical Decision Support & Biomedical Computing**  
> Reference Standards: `CAP, CLSI & ISO Standards`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)

</div>

---

## Overview

**Fastq Quality Trimmer** is a production-grade analytical platform designed for high-throughput evaluation, deterministic verification, and automated compliance checking. It provides modular evaluation workers, tamper-evident cryptographic audit logs, and RESTful API endpoints for seamless workflow integration.

---

## Architecture

```
                    +--------------------------------------------------+
                    |             Air-Gapped Telemetry Stream          |
                    +--------------------------------------------------+
                                             |
                                             v
                    +--------------------------------------------------+
                    |                 SystemSupervisor                 |
                    |       (Zero-PHI Memory & HMAC Audit Trail)       |
                    +--------------------------------------------------+
                               /             |             \
                              /              |              \
                             v               v               v
               +-------------------+ +---------------+ +-------------------+
               | InvariantQCWorker | | SafetyWorker  | | ProtocolWorker    |
               | (Boundary Auditor)| | (Safety Alert)| | (Spec Conformance)|
               +-------------------+ +---------------+ +-------------------+
```

---

## Features

* **Zero-PHI Outbound Interceptors**: AST-level pattern matching preventing sensitive identifier leaks.
* **Tamper-Evident Audit Logging**: Cryptographically linked HMAC-SHA256 records securing every transaction.
* **Multi-Worker Event Loops**: Dedicated verification workers for quality control, safety bounds, and protocol conformance.
* **REST & CLI Interfaces**: Complete FastAPI application and interactive command-line interface.
* **Automated Test Coverage**: Comprehensive test suites verifying boundary conditions and operational stability.

---

## Quick Start (CLI)

```bash
# Run task evaluation
python cli.py audit --task-id TASK-2026-001 --primary 28.5

# System configuration and status query
python cli.py chat "Explain standard reference protocols and calibration limits"

# Verify HMAC-SHA256 audit trail integrity
python cli.py verify-audit

# Launch FastAPI REST Server
python cli.py serve --port 8000
```

---

## API Reference

| Endpoint | Method | Description |
|:---------|:------:|:------------|
| `/health` | `GET` | System health check and metadata |
| `/metrics` | `GET` | Operational metrics exporter |
| `/api/audit` | `POST` | Dispatches task payload across workers and compiles consensus dossier |
| `/api/chat` | `POST` | System query interface |
| `/api/audit/logs` | `GET` | Cryptographic HMAC-SHA256 audit trail log with integrity verification |

---

## Python API Usage

```python
from enrichment import enrichment_suite

# Execute the module suite
results = enrichment_suite.execute_all(primary_val=2.5, secondary_val=1.2)
for module_name, res in results.items():
    print(f"[{res.status}] {res.feature_name} -> Score: {res.score}")
```

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
