from __future__ import annotations
from dataclasses import replace
from typing import Any, Tuple
from node_graph.socket_spec import (
    SocketSpecMeta,
    SocketSpec,
    BaseSocketSpecAPI,
    BaseSpecInferAPI,
)
from aiida_workgraph.orm.mapping import type_mapping
from aiida.engine import Process
from aiida.engine.processes.process_spec import ProcessSpec
from plumpy.ports import Port, PortNamespace
from .socket import TaskSocketNamespace


class SocketSpecAPI(BaseSocketSpecAPI):
    TYPE_MAPPING = type_mapping
    SocketNamespace = TaskSocketNamespace


class SpecInferAPI(BaseSpecInferAPI):

    TYPE_MAPPING = type_mapping

    @classmethod
    def _identifier_from_valid_type(cls, valid_type: Any) -> str:
        """Map AiiDA Port.valid_type -> identifier with our mapping.
        - tuple of types: if len==1 map that; else -> any/default
        - None/empty: any/default
        - single type: mapped identifier
        """
        if isinstance(valid_type, tuple):
            if len(valid_type) == 1:
                return cls._map_identifier(valid_type[0])
            return cls.TYPE_MAPPING["default"]
        if valid_type in (None, Ellipsis):
            return cls.TYPE_MAPPING["default"]
        return cls._map_identifier(valid_type)

    @classmethod
    def _from_port(
        cls, port: Port | PortNamespace, *, parent_required: bool, role: str
    ) -> SocketSpec:
        """Recursively convert an AiiDA Port/PortNamespace to a SocketSpec.
        `role` is "input" or "output" (affects call_role metadata).
        """
        if isinstance(port, PortNamespace):
            required_here = bool(getattr(port, "required", True)) and bool(
                parent_required
            )

            # Build child fields by iterating explicit .ports mapping
            fields: dict[str, SocketSpec] = {}
            for name, child in port.ports.items():
                fields[name] = cls._from_port(
                    child, parent_required=required_here, role=role
                )

            ns = SocketSpec(
                identifier=cls.TYPE_MAPPING["namespace"],
                fields=fields,
                meta=SocketSpecMeta(
                    required=required_here,
                    is_metadata=getattr(port, "is_metadata", False),
                    call_role=("kwargs" if role == "input" else None),
                ),
            )

            # Dynamic namespace? (DynamicPortNamespace derives from PortNamespace)
            is_dyn = bool(getattr(port, "dynamic", False))
            if is_dyn:
                valid_type = getattr(port, "valid_type", None)
                item_ident = cls._identifier_from_valid_type(valid_type)
                ns = replace(ns, dynamic=True, item=SocketSpec(identifier=item_ident))
            return ns

        # Leaf Port (InputPort/OutputPort)
        required_here = bool(getattr(port, "required", True)) and bool(parent_required)
        valid_type = getattr(port, "valid_type", None)
        ident = cls._identifier_from_valid_type(valid_type)
        return SocketSpec(
            identifier=ident,
            meta=SocketSpecMeta(
                required=required_here,
                is_metadata=getattr(port, "is_metadata", False),
                call_role=("kwargs" if role == "input" else None),
            ),
        )

    @classmethod
    def from_aiida_process(
        cls, process_or_spec: type[Process] | Process | ProcessSpec
    ) -> Tuple[SocketSpec, SocketSpec]:
        """Return (inputs_spec, outputs_spec) for an AiiDA Process or its ProcessSpec.

        Accepts:
          - AiiDA Process subclass (e.g. CalcJob, WorkChain)
          - AiiDA Process instance
          - ProcessSpec object (as returned by `.spec()`)
        """
        # Normalize to a ProcessSpec
        if isinstance(process_or_spec, ProcessSpec):
            spec = process_or_spec
        elif isinstance(process_or_spec, type) and issubclass(process_or_spec, Process):
            spec = process_or_spec.spec()
        elif isinstance(process_or_spec, Process):
            spec = process_or_spec.spec()
        else:
            raise TypeError(
                "from_aiida_process expects an AiiDA Process class/instance or a ProcessSpec; "
                f"got {type(process_or_spec)!r}"
            )

        # Validate spec structure
        if not isinstance(spec, ProcessSpec):
            raise TypeError(f".spec() did not return a ProcessSpec; got {type(spec)!r}")

        inputs_ns = getattr(spec, "inputs", None)
        outputs_ns = getattr(spec, "outputs", None)
        if not isinstance(inputs_ns, PortNamespace) or not isinstance(
            outputs_ns, PortNamespace
        ):
            raise TypeError("Spec does not expose PortNamespace for inputs/outputs")

        in_spec = cls._from_port(inputs_ns, parent_required=True, role="input")
        out_spec = cls._from_port(outputs_ns, parent_required=True, role="output")
        # tag top-level outputs with 'return'
        out_spec = replace(out_spec, meta=replace(out_spec.meta, call_role="return"))
        return in_spec, out_spec


socket = SocketSpecAPI.socket
namespace = SocketSpecAPI.namespace
dynamic = SocketSpecAPI.dynamic
expose = SocketSpecAPI.expose
validate_socket_data = SocketSpecAPI.validate_socket_data
infer_specs_from_callable = SpecInferAPI.infer_specs_from_callable
from_aiida_process = SpecInferAPI.from_aiida_process
