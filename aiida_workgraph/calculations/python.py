"""Calcjob to run a Python function on a remote computer."""
from __future__ import annotations

import pathlib
import typing as t

from aiida.common.datastructures import CalcInfo, CodeInfo
from aiida.common.folders import Folder
from aiida.engine import CalcJob, CalcJobProcessSpec
from aiida.orm import Data, SinglefileData

from .general_data import GeneralData

__all__ = ("PythonCalculation",)


class PythonCalculation(CalcJob):
    """Calcjob to run a Python function on a remote computer."""

    _DEFAULT_INPUT_FILE = "script.py"
    _DEFAULT_OUTPUT_FILE = "aiida.out"

    @classmethod
    def define(cls, spec: CalcJobProcessSpec) -> None:  # type: ignore[override]
        """Define the process specification, including its inputs, outputs and known exit codes.

        :param spec: the calculation job process spec to define.
        """
        super().define(spec)
        spec.input("function", required=True)
        spec.input_namespace("kwargs", valid_type=Data, required=False)
        spec.input(
            "code",
            required=True,
            help="Python executable to run the script with.",
        )
        spec.output("results", required=True)
        # set default options (optional)
        spec.inputs["metadata"]["options"]["parser_name"].default = "workgraph.python"
        spec.inputs["metadata"]["options"]["input_filename"].default = "script.py"
        spec.inputs["metadata"]["options"]["output_filename"].default = "aiida.out"
        spec.inputs["metadata"]["options"]["resources"].default = {
            "num_machines": 1,
            "num_mpiprocs_per_machine": 1,
        }

        # start exit codes - marker for docs
        spec.exit_code(
            310,
            "ERROR_READING_OUTPUT_FILE",
            invalidates_cache=True,
            message="The output file could not be read.",
        )
        spec.exit_code(
            320,
            "ERROR_INVALID_OUTPUT",
            invalidates_cache=True,
            message="The output file contains invalid output.",
        )

    def prepare_for_submission(self, folder: Folder) -> CalcInfo:
        """Prepare the calculation for submission.

        1) Write the python script to the folder.
        2) Write the inputs to a pickle file and save it to the folder.

        :param folder: A temporary folder on the local file system.
        :returns: A :class:`aiida.common.datastructures.CalcInfo` instance.
        """
        import inspect
        import cloudpickle as pickle
        import textwrap

        dirpath = pathlib.Path(folder._abspath)
        inputs: dict[str, t.Any]

        if self.inputs.kwargs:
            inputs = dict(self.inputs.kwargs)
        else:
            inputs = {}
        # get the value of pickled function
        function = self.inputs.function.value
        # get the source code of the function
        source_code = inspect.getsource(function)
        source_code_lines = source_code.split("\n")
        function_source_code = "\n".join(source_code_lines[1:])
        function_source_code = textwrap.dedent(function_source_code)
        # create python script to run the function
        script = f"""
import pickle

# define the function
{function_source_code}

# load the inputs from the pickle file
with open('inputs.pickle', 'rb') as handle:
    inputs = pickle.load(handle)

# run the function
result = {function.__name__}(**inputs)
# save the result as a pickle file
with open('result.pickle', 'wb') as handle:
    pickle.dump(result, handle)
"""
        # write the script to the folder
        with folder.open(self.options.input_filename, "w", encoding="utf8") as handle:
            handle.write(script)

        local_copy_list = []
        # create pickle file for the inputs
        input_values = {}
        for key, value in inputs.items():
            if isinstance(value, GeneralData):
                # get the value of the pickled data
                input_values[key] = value.value
            else:
                raise ValueError(f"Unsupported data type: {type(value)}")
            # save the value as a pickle file, the path is absolute
        filename = "inputs.pickle"
        with folder.open(filename, "wb") as handle:
            pickle.dump(input_values, handle)
            # create a singlefiledata object for the pickled data
            file_data = SinglefileData(file=f"{dirpath}/{filename}")
            local_copy_list.append((file_data.uuid, file_data.filename, filename))

        codeinfo = CodeInfo()
        codeinfo.stdin_name = self.options.input_filename
        codeinfo.stdout_name = self.options.output_filename

        if "code" in self.inputs:
            codeinfo.code_uuid = self.inputs.code.uuid

        calcinfo = CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = local_copy_list
        calcinfo.retrieve_list = ["result.pickle", self.options.output_filename]

        return calcinfo
