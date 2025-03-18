from __future__ import annotations

from typing import Any, Dict, List, Optional
from uuid import UUID

from pydantic import BaseModel, Field, conint, constr

# TODO constr is deprecated

class Metadata(BaseModel):
    graph_type: Optional[str] = Field(
        None, description='Type of this graph, e.g. NORMAL'
    )
    group_properties: Optional[List[str]] = None
    group_inputs: Optional[List[str]] = None
    group_outputs: Optional[List[str]] = None

class Executor(BaseModel):
    mode: str
    module_path: Optional[str] = None
    callable_name: Optional[str]
    callable_kind: Optional[str] = None
    graph_data: Optional[str] = None
    pickled_callable: Optional[str] = Field(None, repr=False)
    source_code: Optional[str] = None
    metadata: Optional[str] = None


class Property(BaseModel):
    value: Any
    name: Optional[str] = None
    identifier: Optional[str] = None
    default: Optional[Any] = None
    metadata: Optional[Dict[str, Any]] = None
    arg_type: Optional[str] = None

class SocketLink(BaseModel):
    from_node: str
    from_socket: str

class Socket(BaseModel):
    name: str
    identifier: Optional[str] = None
    link_limit: Optional[float] = None
    links: Optional[List[SocketLink]] = None
    metadata: Optional[Dict[str, Any]] = None
    property: Optional[Property] = None
    sockets: Optional[Dict[str, Socket]] = Field(
        None, description='Nested sockets within this namespace socket'
    )

class Task(BaseModel):
    version: Optional[str] = Field(
        None, description='Version identifier of this node/task data format'
    )
    identifier: str = Field(
        ..., description='Underlying function or process this task references'
    )
    uuid: Optional[UUID] = None
    name: str
    state: Optional[str] = None
    action: Optional[str] = None
    error: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = Field(
        None, description='Metadata about the node/task'
    )
    properties: Optional[Dict[str, Any]] = None
    inputs: Socket = Field(..., description='Inputs socket that contains sub-sockets.')
    outputs: Socket = Field(..., description='Output socket that contains sub-sockets.')
    executor: Optional[Executor] = Field(
        None,
        description='Information about how to execute this task (calcfunction, process, etc.)',
        repr=False
    )
    position: Optional[List[float]] = Field(
        None,
        description='x,y position in a canvas or graph layout',
        max_items=2,
        min_items=2,
    )
    description: Optional[str] = None
    log: Optional[str] = None
    hash: Optional[str] = None
    context_mapping: Optional[Dict[str, Any]] = None
    wait: Optional[List[str]] = None
    children: Optional[List[str]] = None
    execution_count: Optional[int] = None
    parent_task: Optional[List[Optional[str]]] = None
    process: Optional[str] = None
    error_handlers: Optional[Dict[str, Any]] = None


class WorkGraphLink(BaseModel):
    from_socket: str
    from_node: str
    to_socket: str
    to_node: str

class WorkGraphData(BaseModel):
    platform_version: Optional[str] = Field(
        None,
        description='Version identifier for the aiida_workgraph. E.g. aiida_workgraph@0.0.1',
    )
    uuid: Optional[UUID] = Field(
        None, description='A unique identifier (UUID) for this top-level WorkGraph'
    )
    name: str = Field(..., description='Human-readable name for the WorkGraph')
    state: Optional[str] = Field(
        None,
        description='State of the WorkGraph, e.g. CREATED, RUNNING, FINISHED, etc.',
    )
    action: Optional[str] = Field(
        None, description='High-level command or request, e.g. NONE, STOP, etc.'
    )
    error: Optional[str] = Field(
        None, description='Error message if the WorkGraph encountered an issue'
    )
    metadata: Metadata = Field(..., description='Metadata about this WorkGraph')
    links: Optional[List[WorkGraphLink]] = Field(
        None,
        description='List of edges describing how node outputs link to other node inputs',
    )
    description: Optional[str] = None
    context: Optional[Dict[str, Any]] = Field(
        None, description='Any context or ephemeral state that might be needed'
    )
    restart_process: Optional[str] = Field(
        None, description='If relevant, references a process to restart'
    )
    max_iteration: Optional[conint(ge=0)] = None
    execution_count: Optional[conint(ge=0)] = None
    conditions: Optional[List[str]] = None
    max_number_jobs: Optional[conint(ge=0)] = None
    error_handlers: Optional[Dict[str, Any]] = Field(
        None, description='Holds any error handling instructions or metadata'
    )
    tasks: Dict[constr(pattern=r'^[A-Za-z0-9_-]+$'), Task] = Field(
        ..., description='Collection of tasks keyed by name/ID'
    )
