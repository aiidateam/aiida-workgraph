# -*- coding: utf-8 -*-
"""Declaration of FastAPI router for daemon endpoints."""
from __future__ import annotations

import typing as t

from aiida.cmdline.utils.decorators import with_dbenv
from aiida.engine.daemon.client import DaemonException, get_daemon_client
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field


router = APIRouter()


class DaemonStatusModel(BaseModel):
    """Response model for daemon status."""

    running: bool = Field(description="Whether the daemon is running or not.")
    num_workers: t.Optional[int] = Field(
        description="The number of workers if the daemon is running."
    )


@router.get("/api/daemon/status", response_model=DaemonStatusModel)
@with_dbenv()
async def get_daemon_status() -> DaemonStatusModel:
    """Return the daemon status."""
    client = get_daemon_client()

    if not client.is_daemon_running:
        return DaemonStatusModel(running=False, num_workers=None)

    response = client.get_numprocesses()

    return DaemonStatusModel(running=True, num_workers=response["numprocesses"])


@router.get("/api/daemon/worker")
@with_dbenv()
async def get_daemon_worker():
    """Return the daemon status."""
    client = get_daemon_client()

    if not client.is_daemon_running:
        return {}

    response = client.get_worker_info()

    return response["info"]


@router.post("/api/daemon/start", response_model=DaemonStatusModel)
@with_dbenv()
async def get_daemon_start() -> DaemonStatusModel:
    """Start the daemon."""
    client = get_daemon_client()

    if client.is_daemon_running:
        raise HTTPException(status_code=400, detail="The daemon is already running.")

    try:
        client.start_daemon()
    except DaemonException as exception:
        raise HTTPException(status_code=500, detail=str(exception)) from exception

    response = client.get_numprocesses()

    return DaemonStatusModel(running=True, num_workers=response["numprocesses"])


@router.post("/api/daemon/stop", response_model=DaemonStatusModel)
@with_dbenv()
async def get_daemon_stop() -> DaemonStatusModel:
    """Stop the daemon."""
    client = get_daemon_client()

    if not client.is_daemon_running:
        raise HTTPException(status_code=400, detail="The daemon is not running.")

    try:
        client.stop_daemon()
    except DaemonException as exception:
        raise HTTPException(status_code=500, detail=str(exception)) from exception

    return DaemonStatusModel(running=False, num_workers=None)


@router.post("/api/daemon/increase", response_model=DaemonStatusModel)
@with_dbenv()
async def increase_daemon_worker() -> DaemonStatusModel:
    """increase the daemon worker."""
    client = get_daemon_client()

    if not client.is_daemon_running:
        raise HTTPException(status_code=400, detail="The daemon is not running.")

    try:
        client.increase_workers(1)
    except DaemonException as exception:
        raise HTTPException(status_code=500, detail=str(exception)) from exception

    return DaemonStatusModel(running=False, num_workers=None)


@router.post("/api/daemon/decrease", response_model=DaemonStatusModel)
@with_dbenv()
async def decrease_daemon_worker() -> DaemonStatusModel:
    """decrease the daemon worker."""
    client = get_daemon_client()

    if not client.is_daemon_running:
        raise HTTPException(status_code=400, detail="The daemon is not running.")

    try:
        client.decrease_workers(1)
    except DaemonException as exception:
        raise HTTPException(status_code=500, detail=str(exception)) from exception

    return DaemonStatusModel(running=False, num_workers=None)
