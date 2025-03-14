from typing import Any, Callable, Dict, List, Optional, Tuple, Union
from .base import BaseTaskFactory
from aiida_pythonjob import PythonJob
from .function_task import DecoratedFunctionTaskFactory
from aiida_workgraph.tasks.factory.aiida_task import AiiDAComponentTaskFactory
from aiida_workgraph import Task
from aiida import orm
from aiida_pythonjob.data.serializer import general_serializer
from aiida_pythonjob.data.deserializer import deserialize_to_raw_python_data
from aiida.common.extendeddicts import AttributeDict
from node_graph.executor import NodeExecutor


additional_inputs = [
    {
        "identifier": "workgraph.string",
        "name": "computer",
        "metadata": {"is_pythonjob": True},
    },
    {
        "identifier": "workgraph.any",
        "name": "command_info",
        "metadata": {"is_pythonjob": True},
    },
    {
        "identifier": "workgraph.any",
        "name": "register_pickle_by_value",
        "metadata": {"is_pythonjob": True},
    },
]
additional_outputs = [
    {
        "identifier": "workgraph.any",
        "name": "exit_code",
        "metadata": {"is_pythonjob": True},
    }
]


class PythonJobTask(Task):
    """PythonJob Task."""

    identifier = "workgraph.pythonjob"

    def to_dict(self, short: bool = False) -> Dict[str, Any]:
        """Overwrite the to_dict method to handle the PythonJob data.
        Because the data will be passed as input of the WorkGraphEngine,
        all raw data need to be serialized."""
        data = super().to_dict(short=short)
        self.serialize_pythonjob_data(data["inputs"]["sockets"])

        return data

    def update_from_dict(self, data: Dict[str, Any], **kwargs) -> "PythonJobTask":
        """Overwrite the update_from_dict method to handle the PythonJob data."""
        self.deserialize_pythonjob_data(data["inputs"]["sockets"])
        super().update_from_dict(data)

    @classmethod
    def serialize_pythonjob_data(cls, input_data: Dict[str, Any]) -> None:
        """Serialize the properties for PythonJob."""

        for input in input_data.values():
            if not input["metadata"].get("is_pythonjob", False):
                if input["identifier"] == "workgraph.namespace":
                    cls.serialize_pythonjob_data(input["sockets"])
                elif input.get("property", {}).get("value") is not None:
                    cls.serialize_socket_data(input)

    @classmethod
    def deserialize_pythonjob_data(cls, input_data: Dict[str, Any]) -> None:
        """
        Process the task data dictionary for a PythonJob.
        It load the orignal Python data from the AiiDA Data node for the
        args and kwargs of the function.

        Args:
            tdata (Dict[str, Any]): The input data dictionary.

        Returns:
            Dict[str, Any]: The processed data dictionary.
        """

        for input in input_data.values():
            if not input["metadata"].get("is_pythonjob", False):
                if input["identifier"] == "workgraph.namespace":
                    # print("deserialize namespace: ", input["name"])
                    cls.deserialize_pythonjob_data(input["sockets"])
                else:
                    # print("deserialize socket: ", input["name"])
                    cls.deserialize_socket_data(input)

    @classmethod
    def serialize_socket_data(cls, data: Dict[str, Any]) -> Any:
        value = data.get("property", {}).get("value")
        if value is None or isinstance(value, orm.Data):
            return
        data["property"]["value"] = general_serializer(value)

    @classmethod
    def deserialize_socket_data(cls, data: Dict[str, Any]) -> Any:
        value = data.get("property", {}).get("value")
        if isinstance(value, orm.Data):
            data["property"]["value"] = deserialize_to_raw_python_data(value)

    def prepare_for_python_task(self, kwargs: dict, var_kwargs: dict) -> dict:
        """Prepare the inputs for PythonJob"""
        from aiida_pythonjob import prepare_pythonjob_inputs

        function_inputs = kwargs.pop("function_inputs", {})
        for input in self.inputs:
            if not (
                input._metadata.get("is_pythonjob", False)
                or input._metadata.get("is_builtin", False)
            ):
                # if the input is not in the function_inputs, we need try to retrieve it from kwargs
                if input._name not in function_inputs:
                    function_inputs[input._name] = kwargs.pop(input._name, None)
        # if the var_kwargs is not None, we need to pop the var_kwargs from the kwargs
        # then update the function_inputs if var_kwargs is not None
        task_var_kwargs = self.get_args_data()["var_kwargs"]
        if task_var_kwargs is not None:
            function_inputs.pop(task_var_kwargs, None)
            if var_kwargs:
                # var_kwargs can be AttributeDict if it get data from the previous task output
                if isinstance(var_kwargs, (dict, AttributeDict)):
                    function_inputs.update(var_kwargs)
                # otherwise, it should be a Data node
                elif isinstance(var_kwargs, orm.Data):
                    function_inputs.update(var_kwargs.value)
                else:
                    raise ValueError(f"Invalid var_kwargs type: {type(var_kwargs)}")
        # setup code
        code = kwargs.pop("code", None)
        computer = kwargs.pop("computer", "localhost")
        command_info = kwargs.pop("command_info", {})
        register_pickle_by_value = kwargs.pop("register_pickle_by_value", False)
        upload_files = kwargs.pop("upload_files", {})

        metadata = kwargs.pop("metadata", {})
        metadata.update({"call_link_label": self.name})
        # get the function from executor
        func = NodeExecutor(**self.get_executor()).executor
        function_outputs = []
        for output_name in self.outputs._get_all_keys():
            output = self.outputs[output_name]
            if not (
                output._metadata.get("is_pythonjob", False)
                or output._metadata.get("is_builtin", False)
            ):
                # if the output is WORKGRAPH.NAMESPACE, we need to change it to NAMESPACE
                if output._identifier.upper() == "WORKGRAPH.NAMESPACE":
                    function_outputs.append(
                        {"name": output_name, "identifier": "NAMESPACE"}
                    )
                else:
                    function_outputs.append(
                        {"name": output_name, "identifier": output._identifier}
                    )
        # delete workgraph related attributes of the func if exist
        for attr in ["TaskCls", "_ndata", "NodeCls"]:
            if hasattr(func, attr):
                delattr(func, attr)
        inputs = prepare_pythonjob_inputs(
            function=func,
            function_inputs=function_inputs,
            function_outputs=function_outputs,
            code=code,
            command_info=command_info,
            computer=computer,
            metadata=metadata,
            upload_files=upload_files,
            process_label=f"PythonJob<{self.name}>",
            register_pickle_by_value=register_pickle_by_value,
            **kwargs,
        )

        return inputs

    def execute(self, engine_process, args=None, kwargs=None, var_kwargs=None):
        from aiida_workgraph.utils import create_and_pause_process
        from aiida_pythonjob import PythonJob

        inputs = self.prepare_for_python_task(kwargs, var_kwargs)
        inputs["metadata"].update({"call_link_label": self.name})

        kwargs.setdefault("metadata", {})
        kwargs["metadata"].update({"call_link_label": self.name})
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


class PythonJobTaskFactory(BaseTaskFactory):
    """A factory to create PythonJobTask from functions."""

    @classmethod
    def create_class(cls, inputs: dict) -> Union[Task, None]:
        return PythonJobTaskFactory.from_function(
            func=inputs.pop("function", None),
        )

    @classmethod
    def from_function(
        cls,
        func: Callable,
        identifier: Optional[str] = None,
        properties: Optional[List[Tuple[str, str]]] = None,
        inputs: Optional[List[Union[str, dict]]] = None,
        outputs: Optional[List[Union[str, dict]]] = None,
        error_handlers: Optional[List[Dict[str, Any]]] = None,
        catalog: str = "Others",
        additional_data: Optional[Dict[str, Any]] = None,
    ):
        """PythonJobTask is a combination of function task and AiiDA component task."""

        if not hasattr(func, "TaskCls"):
            outputs = outputs or [{"identifier": "workgraph.any", "name": "result"}]
            TaskCls0 = DecoratedFunctionTaskFactory.from_function(
                func,
                identifier=identifier,
                inputs=inputs,
                outputs=outputs,
                properties=properties,
                error_handlers=error_handlers,
            )
        else:
            TaskCls0 = func.TaskCls
        if TaskCls0.node_type.upper() == "GRAPH_BUILDER":
            raise ValueError(
                "GraphBuilder task cannot be run remotely. Please remove 'PythonJob'."
            )
        TaskCls = AiiDAComponentTaskFactory.from_aiida_component(PythonJob)
        tdata = TaskCls0._ndata
        # merge the inputs and outputs from the PythonJob task to the function task
        # skip the already existed inputs and outputs
        for input in additional_inputs:
            tdata["inputs"]["sockets"][input["name"]] = input.copy()
        for input in TaskCls._ndata["inputs"]["sockets"].values():
            if input["name"] not in tdata["inputs"]["sockets"]:
                input["metadata"]["is_pythonjob"] = True
                tdata["inputs"]["sockets"][input["name"]] = input
        for output in TaskCls._ndata["outputs"]["sockets"].values():
            if output["name"] not in tdata["outputs"]["sockets"]:
                output["metadata"]["is_pythonjob"] = True
                tdata["outputs"]["sockets"][output["name"]] = output
        for output in additional_outputs:
            tdata["outputs"]["sockets"][output["name"]] = output.copy()
        # change "copy_files" link_limit to 1e6
        tdata["inputs"]["sockets"]["copy_files"]["link_limit"] = 1e6
        tdata["metadata"]["node_type"] = "PYTHONJOB"
        tdata["metadata"]["catalog"] = catalog
        tdata["default_name"] = func.__name__
        tdata["identifier"] = "workgraph.pythonjob"
        tdata["metadata"]["node_class"] = PythonJobTask
        tdata["metadata"]["class_name"] = "PythonJobTask"
        additional_data = additional_data or {}
        tdata.update(additional_data)

        TaskCls = cls(tdata)
        func.identifier = TaskCls0.identifier
        return TaskCls
