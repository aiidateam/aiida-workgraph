# -*- coding: utf-8 -*-
"""Declaration of FastAPI router for daemon endpoints."""
from __future__ import annotations

import typing as t

from aiida.cmdline.utils.decorators import with_dbenv
from aiida.engine.daemon.client import DaemonException
from aiida_workgraph.engine.scheduler.client import get_scheduler_client
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from aiida_workgraph.engine.scheduler.client import start_scheduler_process


router = APIRouter()


class DaemonStatusModel(BaseModel):
    """Response model for daemon status."""

    running: bool = Field(description="Whether the daemon is running or not.")
    num_workers: t.Optional[int] = Field(
        description="The number of workers if the daemon is running."
    )


@router.get("/api/daemon/scheduler/status", response_model=DaemonStatusModel)
@with_dbenv()
async def get_daemon_status() -> DaemonStatusModel:
    """Return the daemon status."""
    client = get_scheduler_client()

    if not client.is_daemon_running:
        return DaemonStatusModel(running=False, num_workers=None)

    response = client.get_numprocesses()

    return DaemonStatusModel(running=True, num_workers=response["numprocesses"])


@router.get("/api/daemon/scheduler/worker")
@with_dbenv()
async def get_daemon_worker():
    """Return the daemon status."""
    client = get_scheduler_client()

    if not client.is_daemon_running:
        return {}

    response = client.get_worker_info()

    return response["info"]


@router.post("/api/daemon/scheduler/start", response_model=DaemonStatusModel)
@with_dbenv()
async def get_daemon_start() -> DaemonStatusModel:
    """Start the daemon."""
    client = get_scheduler_client()

    if client.is_daemon_running:
        raise HTTPException(status_code=400, detail="The daemon is already running.")

    try:
        client.start_daemon()
    except DaemonException as exception:
        raise HTTPException(status_code=500, detail=str(exception)) from exception

    response = client.get_numprocesses()
    start_scheduler_process()

    return DaemonStatusModel(running=True, num_workers=response["numprocesses"])


@router.post("/api/daemon/scheduler/stop", response_model=DaemonStatusModel)
@with_dbenv()
async def get_daemon_stop() -> DaemonStatusModel:
    """Stop the daemon."""
    client = get_scheduler_client()

    if not client.is_daemon_running:
        raise HTTPException(status_code=400, detail="The daemon is not running.")

    try:
        client.stop_daemon()
    except DaemonException as exception:
        raise HTTPException(status_code=500, detail=str(exception)) from exception

    return DaemonStatusModel(running=False, num_workers=None)


@router.post("/api/daemon/scheduler/increase", response_model=DaemonStatusModel)
@with_dbenv()
async def increase_daemon_worker() -> DaemonStatusModel:
    """increase the daemon worker."""
    client = get_scheduler_client()

    if not client.is_daemon_running:
        raise HTTPException(status_code=400, detail="The daemon is not running.")

    try:
        client.increase_workers(1)
    except DaemonException as exception:
        raise HTTPException(status_code=500, detail=str(exception)) from exception

    response = client.get_numprocesses()
    print(response)
    start_scheduler_process(response["numprocesses"])

    return DaemonStatusModel(running=False, num_workers=None)


@router.post("/api/daemon/scheduler/decrease", response_model=DaemonStatusModel)
@with_dbenv()
async def decrease_daemon_worker() -> DaemonStatusModel:
    """decrease the daemon worker."""
    client = get_scheduler_client()

    if not client.is_daemon_running:
        raise HTTPException(status_code=400, detail="The daemon is not running.")

    try:
        client.decrease_workers(1)
    except DaemonException as exception:
        raise HTTPException(status_code=500, detail=str(exception)) from exception

    return DaemonStatusModel(running=False, num_workers=None)
