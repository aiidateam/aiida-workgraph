from typing import Dict
from aiida_workgraph.node import Node


class AiiDAGather(Node):
    """AiiDAGather"""

    identifier = "AiiDAGather"
    name = "AiiDAGather"
    node_type = "workchain"
    catalog = "AiiDA"
    kwargs = ["datas"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        inp = self.inputs.new("General", "datas")
        inp.link_limit = 100000
        self.outputs.new("General", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida_workgraph.executors.builtin",
            "name": "GatherWorkChain",
        }


class AiiDAToCtx(Node):
    """AiiDAToCtx"""

    identifier = "ToCtx"
    name = "ToCtx"
    node_type = "Control"
    catalog = "AiiDA"
    args = ["key", "value"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("General", "key")
        self.inputs.new("General", "value")
        self.outputs.new("General", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "builtins",
            "name": "setattr",
        }


class AiiDAFromCtx(Node):
    """AiiDAFromCtx"""

    identifier = "FromCtx"
    name = "FromCtx"
    node_type = "Control"
    catalog = "AiiDA"
    args = ["key"]

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        self.inputs.new("General", "key")
        self.outputs.new("General", "result")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "builtins",
            "name": "getattr",
        }


shelljob_inputs = [
    "metadata",
    "metadata.store_provenance",
    "metadata.description",
    "metadata.label",
    "metadata.call_link_label",
    "metadata.dry_run",
    "metadata.computer",
    "metadata.options",
    "metadata.options.input_filename",
    "metadata.options.output_filename",
    "metadata.options.submit_script_filename",
    "metadata.options.scheduler_stdout",
    "metadata.options.scheduler_stderr",
    "metadata.options.resources",
    "metadata.options.max_wallclock_seconds",
    "metadata.options.custom_scheduler_commands",
    "metadata.options.queue_name",
    "metadata.options.rerunnable",
    "metadata.options.account",
    "metadata.options.qos",
    "metadata.options.withmpi",
    "metadata.options.mpirun_extra_params",
    "metadata.options.import_sys_environment",
    "metadata.options.environment_variables",
    "metadata.options.environment_variables_double_quotes",
    "metadata.options.priority",
    "metadata.options.max_memory_kb",
    "metadata.options.prepend_text",
    "metadata.options.append_text",
    "metadata.options.parser_name",
    "metadata.options.additional_retrieve_list",
    "metadata.options.stash",
    "metadata.options.stash.target_base",
    "metadata.options.stash.source_list",
    "metadata.options.stash.stash_mode",
    "metadata.options.redirect_stderr",
    "metadata.options.filename_stdin",
    "metadata.options.additional_retrieve",
    "metadata.options.use_symlinks",
    "code",
    "monitors",
    "remote_folder",
    "nodes",
    "filenames",
    "arguments",
    "outputs",
    "parser",
]

shelljob_outputs = [
    "remote_folder",
    "remote_stash",
    "retrieved",
    "stdout",
    "stderr",
]


class AiiDAShell(Node):
    """AiiDAShell"""

    identifier = "AiiDAShell"
    name = "AiiDAShell"
    node_type = "calcjob"
    catalog = "AiiDA"
    kwargs = shelljob_inputs

    def create_sockets(self) -> None:
        self.inputs.clear()
        self.outputs.clear()
        for inp in shelljob_inputs:
            self.inputs.new("General", inp)
        for out in shelljob_outputs:
            self.outputs.new("General", out)

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida_shell.calculations.shell",
            "name": "ShellJob",
        }
