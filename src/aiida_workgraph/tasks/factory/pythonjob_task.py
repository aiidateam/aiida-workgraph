from typing import Any, Callable, Dict, List, Optional, Tuple, Union
from .base import BaseTaskFactory
from aiida_pythonjob import PythonJob
from aiida_workgraph.tasks.pythonjob import PythonJobTask
from .function_task import DecoratedFunctionTaskFactory
from .aiida_task import AiiDAComponentTaskFactory
from aiida_workgraph import Task


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
        tdata["inputs"].extend(
            [
                {"identifier": "workgraph.string", "name": "computer"},
                {"identifier": "workgraph.any", "name": "command_info"},
            ]
        )
        for input in TaskCls._ndata["inputs"]:
            if input["name"] not in tdata["inputs"]:
                tdata["inputs"].append(input)
        for output in TaskCls._ndata["outputs"]:
            if output["name"] not in tdata["outputs"]:
                tdata["outputs"].append(output)
        tdata["outputs"].append({"identifier": "workgraph.any", "name": "exit_code"})
        # change "copy_files" link_limit to 1e6
        for input in tdata["inputs"]:
            if input["name"] == "copy_files":
                input["link_limit"] = 1e6
        tdata["metadata"]["node_type"] = "PYTHONJOB"
        tdata["metadata"]["catalog"] = catalog
        tdata["identifier"] = "workgraph.pythonjob"
        tdata["metadata"]["node_class"] = PythonJobTask
        tdata["metadata"]["class_name"] = "PythonJobTask"
        additional_data = additional_data or {}
        tdata.update(additional_data)

        TaskCls = cls(tdata)
        func.identifier = TaskCls0.identifier
        return TaskCls
