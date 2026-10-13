"""
FastAPI REST API Server for Fastq Quality Trimmer.
"""
from typing import Dict, Any, List
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from .base import AuditLogger, PHIGuard, SecurityException
from .models import SystemTaskPayload, ConsensusDossier
from .supervisor import SystemSupervisor

supervisor = SystemSupervisor(model_provider="mock")

app = FastAPI(
    title="Fastq Quality Trimmer API",
    description="Auxiliary mock supervisor endpoints (not FASTQ read trimming)",
    version="1.0.0",
)


@app.exception_handler(SecurityException)
async def handle_sensitive_input(request: Request, exc: SecurityException):
    return JSONResponse(status_code=422, content={
        "detail": "Input contains a potentially sensitive identifier and was rejected."
    })


class ChatRequest(BaseModel):
    query: str


@app.get("/health")
def health():
    return {"status": "healthy", "service": "fastq-quality-trimmer", "version": "1.0.0"}


@app.get("/metrics")
def metrics():
    return {
        "dossiers_processed_total": len(supervisor.dossier_registry),
        "audit_blocks_total": len(AuditLogger.get_trail()),
        "system_status": "ready"
    }


@app.post("/api/audit")
def api_audit(payload: SystemTaskPayload):
    dossier = supervisor.process_task(payload)
    return dossier.to_dict()


@app.post("/api/chat")
def api_chat(req: ChatRequest):
    try:
        ans = supervisor.query_supervisory_chat(req.query)
        return {"response": ans}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/audit/logs")
def api_audit_logs():
    return {"audit_trail": AuditLogger.get_trail(), "verified": AuditLogger.verify_integrity()}
