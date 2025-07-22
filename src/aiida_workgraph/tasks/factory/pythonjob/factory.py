from typing import Any, Callable, Dict, List, Optional, Tuple, Union
from aiida_workgraph.tasks.factory.aiida_task import AiiDAComponentTaskFactory
from aiida_workgraph.tasks.factory.function_task import DecoratedFunctionTaskFactory
from aiida_pythonjob import PythonJob
from aiida_pythonjob.calculations.pyfunction import PyFunction
from aiida_workgraph.tasks.factory.pythonjob import PythonJobTask, PyFunctionTask
from aiida_workgraph.tasks.factory.base import BaseTaskFactory
from aiida_workgraph import Task


class BasePythonTaskFactory(BaseTaskFactory):
    """A factory to create PythonJobTask and PyFunctionTask."""

    process_class = None
    task_class = None
    task_type = None
    identifier = None
    additional_inputs = {}
    additional_outputs = {}

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
        """The task is a combination of function task and AiiDA component task."""

        if not hasattr(func, "_TaskCls"):
            outputs = outputs or {
                "result": {"identifier": "workgraph.any", "name": "result"}
            }
            TaskCls0 = DecoratedFunctionTaskFactory.from_function(
                func,
                identifier=identifier,
                inputs=inputs,
                outputs=outputs,
                properties=properties,
                error_handlers=error_handlers,
            )
        else:
            TaskCls0 = func._TaskCls
        if TaskCls0.node_type.upper() == "GRAPH_TASK":
            raise ValueError(
                "Graph task cannot be run remotely. Please remove 'PythonJob'."
            )
        TaskCls = AiiDAComponentTaskFactory.from_aiida_component(cls.process_class)
        tdata = TaskCls0._ndata
        # merge the inputs and outputs from the process_class task to the function task
        # skip the already existed inputs and outputs
        for name, input_data in cls.additional_inputs.items():
            tdata["inputs"]["sockets"][name] = input_data.copy()
        for name, input_data in TaskCls._ndata["inputs"]["sockets"].items():
            if name not in tdata["inputs"]["sockets"]:
                input_data["metadata"].setdefault("extras", {})
                input_data["metadata"]["extras"]["is_pythonjob"] = True
                tdata["inputs"]["sockets"][name] = input_data
        for name, output in TaskCls._ndata["outputs"]["sockets"].items():
            if name not in tdata["outputs"]["sockets"]:
                output["metadata"].setdefault("extras", {})
                output["metadata"]["extras"]["is_pythonjob"] = True
                tdata["outputs"]["sockets"][name] = output
        for name, output in cls.additional_outputs.items():
            tdata["outputs"]["sockets"][name] = output.copy()
        tdata["metadata"]["node_type"] = cls.task_type
        tdata["metadata"]["catalog"] = catalog
        tdata["default_name"] = func.__name__
        tdata["identifier"] = cls.identifier
        tdata["metadata"]["node_class"] = cls.task_class
        additional_data = additional_data or {}
        tdata.update(additional_data)

        TaskCls = cls(tdata)
        func.identifier = TaskCls0.identifier
        return TaskCls


class PythonJobTaskFactory(BasePythonTaskFactory):
    """A factory to create PythonJobTask from functions."""

    process_class = PythonJob
    task_class = PythonJobTask
    task_type = "PythonJob"
    identifier = "workgraph.pythonjob"

    additional_inputs = {
        "computer": {
            "identifier": "workgraph.string",
            "metadata": {"extras": {"is_pythonjob": True}},
        },
        "command_info": {
            "identifier": "workgraph.any",
            "metadata": {"extras": {"is_pythonjob": True}},
        },
        "register_pickle_by_value": {
            "identifier": "workgraph.any",
            "metadata": {"extras": {"is_pythonjob": True}},
        },
    }
    additional_outputs = {
        "exit_code": {
            "identifier": "workgraph.any",
            "metadata": {"extras": {"is_pythonjob": True}},
        }
    }


class PyFunctionTaskFactory(BasePythonTaskFactory):

    process_class = PyFunction
    task_class = PyFunctionTask
    task_type = "PyFunction"
    identifier = "workgraph.pyfunction"

    additional_inputs = {}
    additional_outputs = {
        "exit_code": {
            "identifier": "workgraph.any",
            "metadata": {"is_pythonjob": True},
        }
    }
