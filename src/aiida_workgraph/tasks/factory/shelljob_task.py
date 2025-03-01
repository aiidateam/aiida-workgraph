from typing import List, Optional, Union
from .base import BaseTaskFactory
from aiida_shell import ShellJob
from .aiida_task import AiiDAComponentTaskFactory
from aiida_workgraph.utils import validate_task_inout


class ShellJobTaskFactory(BaseTaskFactory):
    """A factory to create ShellJobTask."""

    @classmethod
    def create_task(
        cls,
        outputs: Optional[List[Union[str, dict]]] = None,
        parser_outputs: Optional[List[Union[str, dict]]] = None,
    ):
        from aiida_shell.parsers.shell import ShellParser

        outputs = outputs or []
        parser_outputs = parser_outputs or []
        parser_outputs = validate_task_inout(parser_outputs, "parser_outputs")
        TaskCls = AiiDAComponentTaskFactory.from_aiida_component(ShellJob)
        tdata = TaskCls._ndata
        tdata["outputs"].extend(
            [
                {"identifier": "workgraph.any", "name": "stdout"},
                {"identifier": "workgraph.any", "name": "stderr"},
            ]
        )
        outputs = [
            {
                "identifier": "workgraph.any",
                "name": ShellParser.format_link_label(output),
            }
            for output in outputs
        ]
        outputs.extend(parser_outputs)
        for output in outputs:
            if output["name"] not in tdata["outputs"]:
                tdata["outputs"].append(output)
        #
        tdata["identifier"] = "ShellJob"
        tdata["inputs"].extend(
            [
                {"identifier": "workgraph.any", "name": "command"},
                {"identifier": "workgraph.any", "name": "resolve_command"},
            ]
        )
        tdata["metadata"]["node_type"] = "SHELLJOB"

        TaskCls = cls(tdata)
        return TaskCls
