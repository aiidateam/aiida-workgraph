from typing import Any, Dict
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida_workgraph import Task
from aiida_pythonjob.data.serializer import general_serializer
from aiida_pythonjob.data.deserializer import deserialize_to_raw_python_data
from aiida_pythonjob import PythonJob
from aiida_workgraph.utils import create_and_pause_process
from node_graph.executor import NodeExecutor
from aiida.engine import run_get_node
from aiida_pythonjob import pyfunction


class BaseSerializablePythonTask(Task):
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
        self._serialize_python_data(data["inputs"]["sockets"])

    def update_from_dict(
        self, data: Dict[str, Any], **kwargs
    ) -> "BaseSerializablePythonTask":
        """
        Called when reloading from a dict. In principle, one could also do
        `_deserialize_python_data` here if you want to re-inflate the raw python,
        but many workflows prefer deferring that until run-time.
        """
        super().update_from_dict(data, **kwargs)
        return self

    @classmethod
    def _serialize_python_data(cls, input_sockets: Dict[str, Any]) -> None:
        """
        Recursively walk over the sockets and convert raw Python
        values to AiiDA Data nodes, if needed.
        """
        for socket in input_sockets.values():
            if not socket["metadata"].get("extras", {}).get("is_pythonjob", False):
                if socket["identifier"] == "workgraph.namespace":
                    cls._serialize_python_data(socket["sockets"])
                elif socket.get("property", {}).get("value") is not None:
                    cls._serialize_socket_data(socket)

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
    def _serialize_socket_data(cls, socket: Dict[str, Any]) -> Any:
        value = socket.get("property", {}).get("value")
        if value is None or isinstance(value, orm.Data):
            return  # Already stored or is None
        socket["property"]["value"] = general_serializer(value)

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

    def build_function_ports(self, socket):
        # Build an explicit list of function outputs
        port = {"name": socket._name, "identifier": socket._identifier}
        if hasattr(socket, "_sockets"):
            port["ports"] = []
            for name, sub_socket in socket._sockets.items():
                if not (
                    sub_socket._metadata.extras.get("is_pythonjob", False)
                    or sub_socket._metadata.builtin_socket
                ):
                    if sub_socket._identifier.upper() == "WORKGRAPH.NAMESPACE":
                        port["ports"].append(self.build_function_ports(sub_socket))
                    else:
                        port["ports"].append(
                            {"name": name, "identifier": sub_socket._identifier}
                        )
        return port


class PythonJobTask(BaseSerializablePythonTask):
    """PythonJob Task."""

    identifier = "workgraph.pythonjob"

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        """
        Here is the specialized 'execute' method for PythonJobTask,
        including the 'prepare_for_python_task' logic.
        """
        from aiida_pythonjob import prepare_pythonjob_inputs

        # Prepare inputs (similar to the old 'prepare_for_python_task' method):
        function_inputs = kwargs.pop("function_inputs", {}) or {}
        for socket_in in self.inputs:
            if not (
                socket_in._metadata.extras.get("is_pythonjob", False)
                or socket_in._metadata.builtin_socket
            ):
                if socket_in._name not in function_inputs:
                    function_inputs[socket_in._name] = kwargs.pop(socket_in._name, None)

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

        # Pull out code, computer, etc
        computer = kwargs.pop("computer", "localhost")
        command_info = kwargs.pop("command_info", {})
        register_pickle_by_value = kwargs.pop("register_pickle_by_value", False)
        upload_files = kwargs.pop("upload_files", {})

        metadata = kwargs.pop("metadata", {})
        metadata.update({"call_link_label": self.name})

        # Resolve the actual function from the NodeExecutor
        func = NodeExecutor(**self.get_executor()).executor
        if hasattr(func, "_TaskCls") and hasattr(func, "_func"):
            func = func._func

        if hasattr(func, "is_process_function"):
            func = func.func

        input_ports = self.build_function_ports(self.inputs)
        output_ports = self.build_function_ports(self.outputs)

        # Prepare the final inputs for PythonJob
        inputs = prepare_pythonjob_inputs(
            function=func,
            function_inputs=function_inputs,
            input_ports=input_ports["ports"],
            output_ports=output_ports["ports"],
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
        executor = NodeExecutor(**self.get_executor()).executor
        # If it's a wrapped function, unwrap
        if hasattr(executor, "_NodeCls") and hasattr(executor, "_func"):
            executor = executor._func
        # Make sure it's process_function-decorated
        if not hasattr(executor, "is_process_function"):
            executor = pyfunction()(executor)

        kwargs = kwargs or {}
        kwargs.setdefault("metadata", {})
        input_ports = self.build_function_ports(self.inputs)
        output_ports = self.build_function_ports(self.outputs)

        kwargs["metadata"].update({"call_link_label": self.name})
        kwargs["input_ports"] = input_ports["ports"]
        kwargs["output_ports"] = output_ports["ports"]

        # If we have var_kwargs, pass them in
        if var_kwargs is None:
            _, process = run_get_node(executor, **kwargs)
        else:
            _, process = run_get_node(executor, **kwargs, **var_kwargs)

        process.label = self.name
        return process, "FINISHED"
