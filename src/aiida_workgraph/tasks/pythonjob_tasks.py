from __future__ import annotations
from typing import Any, Dict, Optional, Callable
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida_pythonjob.data.serializer import general_serializer, all_serializers
from aiida_pythonjob.data.deserializer import deserialize_to_raw_python_data
from aiida_workgraph.utils import create_and_pause_process
from aiida.engine import run_get_node
from aiida_pythonjob import pyfunction, PythonJob, PyFunction
from aiida_workgraph.task import SpecTask
from node_graph.socket_spec import SocketSpec
from node_graph.node_spec import NodeSpec
from aiida_workgraph.socket_spec import namespace
from .function_task import build_callable_nodespec
from node_graph.error_handler import ErrorHandlerSpec


class BaseSerializablePythonTask(SpecTask):
    """
    A base Task that handles serialization and deserialization
    of Python data into AiiDA Data nodes, so that raw Python data
    can be stored/passed around by the WorkGraph engine.
    Subclasses must implement their own `execute` method.
    """

    def serialize_data(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Called during Task -> dict conversion. We walk over the input sockets
        and run our specialized Python serialization.
        """
        non_function_inputs = self.non_function_inputs
        use_pickle = self.inputs.metadata.use_pickle.value
        for name, socket in data["inputs"]["sockets"].items():
            # skip metadata etc
            if name in non_function_inputs:
                continue
            if socket["identifier"] == "workgraph.namespace":
                self._serialize_python_data(socket["sockets"], use_pickle=use_pickle)
            else:
                self._serialize_socket_data(socket, use_pickle=use_pickle)

    def update_from_dict(
        self, data: Dict[str, Any], **kwargs
    ) -> "BaseSerializablePythonTask":
        """
        Called when reloading from a dict. Note, we do not run `_deserialize_python_data` here.
        Thus, the value of the socket will be AiiDA data nodes.
        """
        super().update_from_dict(data, **kwargs)
        return self

    @classmethod
    def _serialize_python_data(
        cls, input_sockets: Dict[str, Any], use_pickle: bool | None = None
    ) -> None:
        """
        Recursively walk over the sockets and convert raw Python
        values to AiiDA Data nodes, if needed.
        """
        for socket in input_sockets.values():
            if not socket["metadata"].get("extras", {}).get("is_pythonjob", False):
                if socket["identifier"] == "workgraph.namespace":
                    cls._serialize_python_data(socket["sockets"], use_pickle=use_pickle)
                elif socket.get("property", {}).get("value") is not None:
                    cls._serialize_socket_data(socket, use_pickle=use_pickle)

    @classmethod
    def _deserialize_python_data(cls, input_sockets: Dict[str, Any]) -> None:
        """
        Recursively walk over the sockets and convert AiiDA Data nodes
        back into raw Python objects, if needed.
        """
        for socket in input_sockets.values():
            if not socket["metadata"].get("extras", {}).get("is_pythonjob", False):
                if socket["identifier"] == "workgraph.namespace":
                    cls._deserialize_python_data(socket["sockets"])
                else:
                    cls._deserialize_socket_data(socket)

    @classmethod
    def _serialize_socket_data(
        cls, socket: Dict[str, Any], use_pickle: bool | None = None
    ) -> Any:
        value = socket.get("property", {}).get("value")
        if value is None or isinstance(value, orm.Data):
            return  # Already stored or is None
        socket["property"]["value"] = general_serializer(
            value, serializers=all_serializers, use_pickle=use_pickle
        )

    @classmethod
    def _deserialize_socket_data(cls, socket: Dict[str, Any]) -> Any:
        value = socket.get("property", {}).get("value")
        if isinstance(value, orm.Data):
            socket["property"]["value"] = deserialize_to_raw_python_data(value)

    def execute(self, *args, **kwargs):
        """
        Subclasses must override.
        """
        raise NotImplementedError("Subclasses must implement `execute`.")

    def get_origin_specs(
        self, label: str = "function"
    ) -> tuple[Optional[SocketSpec], Optional[SocketSpec]]:
        block = self._spec.metadata.get("specs", {}).get(label, {})
        in_d = block.get("inputs")
        out_d = block.get("outputs")
        return (
            SocketSpec.from_dict(in_d) if in_d else None,
            SocketSpec.from_dict(out_d) if out_d else None,
        )

    @property
    def non_function_inputs(self):
        process_in, _ = self.get_origin_specs(label="process")
        additions_in, _ = self.get_origin_specs(label="additional")
        keys = []
        if process_in:
            keys.extend(process_in.fields.keys())
        if additions_in:
            keys.extend(additions_in.fields.keys())
        return keys


class PythonJobTask(BaseSerializablePythonTask):
    """PythonJob Task."""

    identifier = "workgraph.pythonjob"

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        """
        Here is the specialized 'execute' method for PythonJobTask,
        including the 'prepare_for_python_task' logic.
        """
        from aiida_pythonjob import prepare_pythonjob_inputs

        inputs_spec, outputs_spec = self.get_origin_specs()
        # Pull out code, computer, etc
        computer = kwargs.pop("computer", "localhost")
        if isinstance(computer, orm.Str):
            computer = computer.value
        command_info = kwargs.pop("command_info", {})
        register_pickle_by_value = kwargs.pop("register_pickle_by_value", False)
        upload_files = kwargs.pop("upload_files", {})
        metadata = kwargs.pop("metadata", {})
        metadata.update({"call_link_label": self.name})
        # Prepare inputs (similar to the old 'prepare_for_python_task' method):
        function_inputs = kwargs.pop("function_inputs", {}) or {}
        non_function_inputs = self.non_function_inputs
        for key in list(kwargs.keys()):
            if key not in non_function_inputs:
                function_inputs[key] = kwargs.pop(key)
        # Handle var_kwargs
        if self.get_args_data()["var_kwargs"] is not None:
            var_key = self.get_args_data()["var_kwargs"]
            function_inputs.pop(var_key, None)
            if var_kwargs:
                if isinstance(var_kwargs, (dict, AttributeDict)):
                    function_inputs.update(var_kwargs)
                elif isinstance(var_kwargs, orm.Data):
                    function_inputs.update(var_kwargs.value)
                else:
                    raise ValueError(f"Invalid var_kwargs type: {type(var_kwargs)}")
        # Resolve the actual function from the NodeExecutor
        func = self.get_executor().executor
        if hasattr(func, "_TaskCls") and hasattr(func, "_func"):
            func = func._func

        if hasattr(func, "is_process_function"):
            func = func.func

        # Prepare the final inputs for PythonJob
        inputs = prepare_pythonjob_inputs(
            function=func,
            function_inputs=function_inputs,
            inputs_spec=inputs_spec,
            outputs_spec=outputs_spec,
            code=kwargs.pop("code", None),
            command_info=command_info,
            computer=computer,
            metadata=metadata,
            upload_files=upload_files,
            process_label=f"PythonJob<{self.name}>",
            register_pickle_by_value=register_pickle_by_value,
            **kwargs,
        )

        # If we want to pause
        if self.action == "PAUSE":
            engine_process.report(f"Task {self.name} is created and paused.")
            process = create_and_pause_process(
                engine_process.runner,
                PythonJob,
                inputs,
                state_msg="Paused through WorkGraph",
            )
            state = "CREATED"
            process = process.node
        else:
            process = engine_process.submit(PythonJob, **inputs)
            state = "RUNNING"

        process.label = self.name
        return process, state


class PyFunctionTask(BaseSerializablePythonTask):
    """PyFunction Task."""

    identifier = "workgraph.pyfunction"

    def execute(self, args=None, kwargs=None, var_kwargs=None):
        """
        Specialized 'execute' for PyFunction. The logic is simpler than PythonJob.
        """
        from node_graph.node_spec import BaseHandle

        executor = self.get_executor().executor
        # If it's a wrapped function, unwrap
        if isinstance(executor, BaseHandle) and hasattr(executor, "_func"):
            executor = executor._func
        # Make sure it's process_function-decorated
        if not hasattr(executor, "is_process_function"):
            executor = pyfunction()(executor)

        kwargs = kwargs or {}
        kwargs.setdefault("metadata", {})
        kwargs["metadata"].update({"call_link_label": self.name})
        inputs_spec, outputs_spec = self.get_origin_specs()
        kwargs["inputs_spec"] = inputs_spec
        kwargs["outputs_spec"] = outputs_spec

        # If we have var_kwargs, pass them in
        if var_kwargs is None:
            _, process = run_get_node(executor, **kwargs)
        else:
            _, process = run_get_node(executor, **kwargs, **var_kwargs)

        process.label = self.name
        return process, "FINISHED"


def _build_pythonjob_nodespec(
    obj: Callable,
    identifier: Optional[str] = None,
    in_spec: Optional[SocketSpec] | list = None,
    out_spec: Optional[SocketSpec] | list = None,
    error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
) -> NodeSpec:
    # allow list specs just for PythonJob (keep existing behavior)
    from aiida_workgraph.socket_spec import validate_socket_data

    in_spec = validate_socket_data(in_spec)
    out_spec = validate_socket_data(out_spec)

    # additions specific to PythonJob
    add_in = namespace(
        computer=str,
        command_info=dict,
        register_pickle_by_value=bool,
    )

    return build_callable_nodespec(
        obj=obj,
        node_type="PYTHONJOB",
        base_class=PythonJobTask,
        identifier=identifier,
        process_cls=PythonJob,
        in_spec=in_spec,
        out_spec=out_spec,
        add_inputs=add_in,
        error_handlers=error_handlers,
    )


def _build_pyfunction_nodespec(
    obj: Callable,
    identifier: Optional[str] = None,
    in_spec: Optional[SocketSpec] = None,
    out_spec: Optional[SocketSpec] = None,
    error_handlers: Optional[Dict[str, ErrorHandlerSpec]] = None,
) -> NodeSpec:
    return build_callable_nodespec(
        obj=obj,
        node_type="PYFUNCTION",
        base_class=PyFunctionTask,
        identifier=identifier,
        process_cls=PyFunction,
        in_spec=in_spec,
        out_spec=out_spec,
        error_handlers=error_handlers,
    )
