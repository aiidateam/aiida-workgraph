from __future__ import annotations
from typing import Any, Dict, Optional, Callable, Annotated
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida_pythonjob.data.serializer import all_serializers
from aiida_pythonjob.data.deserializer import deserialize_to_raw_python_data
from aiida_workgraph.utils import create_and_pause_process
from aiida.engine import run_get_node
from aiida_pythonjob import pyfunction, PythonJob, PyFunction
from aiida_pythonjob.utils import serialize_ports
from aiida_workgraph.task import SpecTask
from node_graph.socket_spec import SocketSpec, SocketSpecSelect, SocketSpecMeta
from node_graph.node_spec import NodeSpec
from aiida_workgraph.socket_spec import namespace
from .function_task import build_callable_nodespec
from node_graph.error_handler import ErrorHandlerSpec
from node_graph.executor import RuntimeExecutor


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
        function_inputs = {
            key: data["inputs"][key]
            for key in data["inputs"]
            if key not in self.non_function_inputs
        }
        serialized_inputs = serialize_ports(
            python_data=function_inputs,
            port_schema=self._spec.inputs,
            serializers=all_serializers,
        )
        data["inputs"].update(serialized_inputs)

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
    def _deserialize_socket_data(cls, socket: Dict[str, Any]) -> Any:
        value = socket.get("property", {}).get("value")
        if isinstance(value, orm.Data):
            socket["property"]["value"] = deserialize_to_raw_python_data(value)

    def execute(self, *args, **kwargs):
        """
        Subclasses must override.
        """
        raise NotImplementedError("Subclasses must implement `execute`.")

    @property
    def non_function_inputs(self):
        return self._spec.metadata.get("non_function_inputs", [])

    @property
    def non_function_outputs(self):
        return self._spec.metadata.get("non_function_outputs", [])


class PythonJobTask(BaseSerializablePythonTask):
    """PythonJob Task."""

    identifier = "workgraph.pythonjob"

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        """
        Here is the specialized 'execute' method for PythonJobTask,
        including the 'prepare_for_python_task' logic.
        """
        from aiida_pythonjob import prepare_pythonjob_inputs

        inputs_spec = namespace(
            _=Annotated[
                Any,
                self._spec.inputs,
                SocketSpecSelect(exclude=self.non_function_inputs),
            ]
        ).fields["_"]
        outputs_spec = namespace(
            _=Annotated[
                Any,
                self._spec.outputs,
                SocketSpecSelect(exclude=self.non_function_outputs),
            ]
        ).fields["_"]

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
        for key in list(kwargs.keys()):
            if key not in self.non_function_inputs:
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
        executor = RuntimeExecutor(**self.get_executor().to_dict())
        func = executor.callable
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

        executor = RuntimeExecutor(**self.get_executor().to_dict()).callable
        # If it's a wrapped function, unwrap
        if isinstance(executor, BaseHandle) and hasattr(executor, "_func"):
            executor = executor._func
        # Make sure it's process_function-decorated
        if not hasattr(executor, "is_process_function"):
            executor = pyfunction()(executor)

        kwargs = kwargs or {}
        kwargs.setdefault("metadata", {})
        kwargs["metadata"].update({"call_link_label": self.name})
        inputs_spec = namespace(
            _=Annotated[
                Any,
                self._spec.inputs,
                SocketSpecSelect(exclude=self.non_function_inputs),
            ]
        ).fields["_"]
        outputs_spec = namespace(
            _=Annotated[
                Any,
                self._spec.outputs,
                SocketSpecSelect(exclude=self.non_function_outputs),
            ]
        ).fields["_"]

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
        computer=Annotated[str, SocketSpecMeta(required=False)],
        command_info=Annotated[dict, SocketSpecMeta(required=False)],
        register_pickle_by_value=Annotated[bool, SocketSpecMeta(required=False)],
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
