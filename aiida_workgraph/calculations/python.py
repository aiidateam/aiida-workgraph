"""Calcjob to run a Python function on a remote computer."""
from __future__ import annotations

import pathlib
import typing as t

from aiida.common.datastructures import CalcInfo, CodeInfo
from aiida.common.folders import Folder
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import CalcJob, CalcJobProcessSpec
from aiida.orm import (
    Data,
    SinglefileData,
    Str,
    List,
    FolderData,
    RemoteData,
    to_aiida_type,
)


__all__ = ("PythonJob",)


class PythonJob(CalcJob):
    """Calcjob to run a Python function on a remote computer."""

    _internal_retrieve_list = []
    _retrieve_singlefile_list = []
    _retrieve_temporary_list = []

    _DEFAULT_INPUT_FILE = "script.py"
    _DEFAULT_OUTPUT_FILE = "aiida.out"

    _default_parser = "workgraph.python"

    @classmethod
    def define(cls, spec: CalcJobProcessSpec) -> None:  # type: ignore[override]
        """Define the process specification, including its inputs, outputs and known exit codes.

        :param spec: the calculation job process spec to define.
        """
        super().define(spec)
        spec.input(
            "function_source_code",
            valid_type=Str,
            serializer=to_aiida_type,
            required=False,
        )
        spec.input(
            "function_name", valid_type=Str, serializer=to_aiida_type, required=False
        )
        spec.input_namespace(
            "function_kwargs", valid_type=Data, required=False
        )  # , serializer=serialize_to_aiida_nodes)
        spec.input(
            "output_info",
            valid_type=List,
            required=False,
            serializer=to_aiida_type,
            help="The information of the output ports",
        )
        spec.input(
            "parent_folder",
            valid_type=(RemoteData, FolderData, SinglefileData),
            required=False,
            help="Use a local or remote folder as parent folder (for restarts and similar)",
        )
        spec.input(
            "parent_folder_name",
            valid_type=Str,
            default=lambda: Str("./parent_folder/"),
            required=False,
            serializer=to_aiida_type,
            help="""Default name of the subfolder that you want to create in the working directory,
            in which you want to place the files taken from parent_folder""",
        )
        spec.input(
            "parent_output_folder",
            valid_type=Str,
            default=None,
            required=False,
            serializer=to_aiida_type,
            help="Name of the subfolder inside 'parent_folder' from which you want to copy the files",
        )
        spec.input_namespace(
            "upload_files",
            valid_type=(FolderData, SinglefileData),
            required=False,
            help="The folder/files to upload",
        )
        spec.input_namespace(
            "copy_files",
            valid_type=(RemoteData,),
            required=False,
            help="The folder/files to copy from the remote computer",
        )
        spec.input(
            "additional_retrieve_list",
            valid_type=List,
            default=None,
            required=False,
            serializer=to_aiida_type,
            help="The names of the files to retrieve",
        )
        spec.outputs.dynamic = True
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
        spec.exit_code(
            321,
            "ERROR_RESULT_OUTPUT_MISMATCH",
            invalidates_cache=True,
            message="The number of results does not match the number of outputs.",
        )

    def _build_process_label(self) -> str:
        """Use the function name as the process label.

        :returns: The process label to use for ``ProcessNode`` instances.
        """
        return f"PythonJob<{self.inputs.function_name.value}>"

    def prepare_for_submission(self, folder: Folder) -> CalcInfo:
        """Prepare the calculation for submission.

        1) Write the python script to the folder.
        2) Write the inputs to a pickle file and save it to the folder.

        :param folder: A temporary folder on the local file system.
        :returns: A :class:`aiida.common.datastructures.CalcInfo` instance.
        """
        import cloudpickle as pickle

        dirpath = pathlib.Path(folder._abspath)
        inputs: dict[str, t.Any]

        if self.inputs.function_kwargs:
            inputs = dict(self.inputs.function_kwargs)
        else:
            inputs = {}
        # get the value of pickled function
        function_source_code = self.inputs.function_source_code.value
        # create python script to run the function
        script = f"""
import pickle

# define the function
{function_source_code}

# load the inputs from the pickle file
with open('inputs.pickle', 'rb') as handle:
    inputs = pickle.load(handle)

# run the function
result = {self.inputs.function_name.value}(**inputs)
# save the result as a pickle file
with open('results.pickle', 'wb') as handle:
    pickle.dump(result, handle)
"""
        # write the script to the folder
        with folder.open(self.options.input_filename, "w", encoding="utf8") as handle:
            handle.write(script)
        # symlink = settings.pop('PARENT_FOLDER_SYMLINK', False)
        symlink = True

        remote_copy_list = []
        local_copy_list = []
        remote_symlink_list = []
        remote_list = remote_symlink_list if symlink else remote_copy_list

        source = self.inputs.get("parent_folder", None)

        if source is not None:
            if isinstance(source, RemoteData):
                dirpath = pathlib.Path(source.get_remote_path())
                if self.inputs.parent_output_folder is not None:
                    dirpath = (
                        pathlib.Path(source.get_remote_path())
                        / self.inputs.parent_output_folder.value
                    )
                remote_list.append(
                    (
                        source.computer.uuid,
                        str(dirpath),
                        self.inputs.parent_folder_name.value,
                    )
                )
            elif isinstance(source, FolderData):
                dirname = (
                    self.inputs.parent_output_folder.value
                    if self.inputs.parent_output_folder is not None
                    else ""
                )
                local_copy_list.append(
                    (source.uuid, dirname, self.inputs.parent_folder_name.value)
                )
            elif isinstance(source, SinglefileData):
                local_copy_list.append((source.uuid, source.filename, source.filename))
        if self.inputs.upload_files:
            upload_files = self.inputs.upload_files
            for key, source in upload_files.items():
                # replace "_dot_" with "." in the key
                key = key.replace("_dot_", ".")
                if isinstance(source, FolderData):
                    local_copy_list.append((source.uuid, "", key))
                elif isinstance(source, SinglefileData):
                    local_copy_list.append(
                        (source.uuid, source.filename, source.filename)
                    )
                else:
                    raise ValueError(
                        f"""Input folder/file: {source} is not supported.
Only AiiDA SinglefileData and FolderData are allowed."""
                    )
        if "copy_files" in self.inputs:
            copy_files = self.inputs.copy_files
            for key, source in copy_files.items():
                # replace "_dot_" with "." in the key
                key = key.replace("_dot_", ".")
                dirpath = pathlib.Path(source.get_remote_path())
                remote_list.append((source.computer.uuid, str(dirpath), key))
        # create pickle file for the inputs
        input_values = {}
        for key, value in inputs.items():
            if isinstance(value, Data) and hasattr(value, "value"):
                # get the value of the pickled data
                input_values[key] = value.value
            # TODO: should check this recursively
            elif isinstance(value, (AttributeDict, dict)):
                # if the value is an AttributeDict, use recursively
                input_values[key] = {k: v.value for k, v in value.items()}
            else:
                raise ValueError(
                    f"Input data {value} is not supported. Only AiiDA data Node with a value attribute is allowed. "
                )
        # save the value as a pickle file, the path is absolute
        filename = "inputs.pickle"
        dirpath = pathlib.Path(folder._abspath)
        with folder.open(filename, "wb") as handle:
            pickle.dump(input_values, handle)
            # create a singlefiledata object for the pickled data
            file_data = SinglefileData(file=f"{dirpath}/{filename}")
            local_copy_list.append((file_data.uuid, file_data.filename, filename))

        codeinfo = CodeInfo()
        codeinfo.stdin_name = self.options.input_filename
        codeinfo.stdout_name = self.options.output_filename
        codeinfo.code_uuid = self.inputs.code.uuid

        calcinfo = CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.remote_symlink_list = remote_symlink_list
        calcinfo.retrieve_list = ["results.pickle", self.options.output_filename]
        if self.inputs.additional_retrieve_list is not None:
            calcinfo.retrieve_list += self.inputs.additional_retrieve_list.get_list()
        calcinfo.retrieve_list += self._internal_retrieve_list

        calcinfo.retrieve_temporary_list = self._retrieve_temporary_list
        calcinfo.retrieve_singlefile_list = self._retrieve_singlefile_list

        return calcinfo
