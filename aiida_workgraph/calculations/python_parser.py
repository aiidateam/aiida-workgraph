"""Parser for an `PythonJob` job."""
from aiida.parsers.parser import Parser
from aiida_workgraph.orm import general_serializer
from aiida.engine import ExitCode


class PythonParser(Parser):
    """Parser for an `PythonJob` job."""

    def parse(self, **kwargs):
        """Parse the contents of the output files stored in the `retrieved` output node.

        The function_outputs could be a namespce, e.g.,
        function_outputs=[
            {"identifier": "workgraph.namespace", "name": "add_multiply"},
            {"name": "add_multiply.add"},
            {"name": "add_multiply.multiply"},
            {"name": "minus"},
        ]
        """
        import pickle

        function_outputs = self.node.inputs.function_outputs.get_list()
        # function_outputs exclude ['_wait', '_outputs', 'remote_folder', 'remote_stash', 'retrieved']
        self.output_list = [
            data
            for data in function_outputs
            if data["name"]
            not in [
                "_wait",
                "_outputs",
                "remote_folder",
                "remote_stash",
                "retrieved",
                "exit_code",
            ]
        ]
        # first we remove nested outputs, e.g., "add_multiply.add"
        top_level_output_list = [
            output for output in self.output_list if "." not in output["name"]
        ]
        exit_code = 0
        try:
            with self.retrieved.base.repository.open("results.pickle", "rb") as handle:
                results = pickle.load(handle)
                if isinstance(results, tuple):
                    if len(top_level_output_list) != len(results):
                        self.exit_codes.ERROR_RESULT_OUTPUT_MISMATCH
                    for i in range(len(top_level_output_list)):
                        top_level_output_list[i]["value"] = self.serialize_output(
                            results[i], top_level_output_list[i]
                        )
                elif isinstance(results, dict) and len(top_level_output_list) > 1:
                    # pop the exit code if it exists
                    exit_code = results.pop("exit_code", 0)
                    for output in top_level_output_list:
                        if output.get("required", False):
                            if output["name"] not in results:
                                self.exit_codes.ERROR_MISSING_OUTPUT
                        output["value"] = self.serialize_output(
                            results.pop(output["name"]), output
                        )
                    # if there are any remaining results, raise an warning
                    if results:
                        self.logger.warning(
                            f"Found extra results that are not included in the output: {results.keys()}"
                        )
                elif isinstance(results, dict) and len(top_level_output_list) == 1:
                    exit_code = results.pop("exit_code", 0)
                    # if output name in results, use it
                    if top_level_output_list[0]["name"] in results:
                        top_level_output_list[0]["value"] = self.serialize_output(
                            results[top_level_output_list[0]["name"]],
                            top_level_output_list[0],
                        )
                    # otherwise, we assume the results is the output
                    else:
                        top_level_output_list[0]["value"] = self.serialize_output(
                            results, top_level_output_list[0]
                        )
                elif len(top_level_output_list) == 1:
                    # otherwise, we assume the results is the output
                    top_level_output_list[0]["value"] = self.serialize_output(
                        results, top_level_output_list[0]
                    )
                else:
                    raise ValueError(
                        "The number of results does not match the number of outputs."
                    )
                for output in top_level_output_list:
                    self.out(output["name"], output["value"])
                if exit_code:
                    if isinstance(exit_code, dict):
                        exit_code = ExitCode(exit_code["status"], exit_code["message"])
                    elif isinstance(exit_code, int):
                        exit_code = ExitCode(exit_code)
                    return exit_code
        except OSError:
            return self.exit_codes.ERROR_READING_OUTPUT_FILE
        except ValueError as exception:
            self.logger.error(exception)
            return self.exit_codes.ERROR_INVALID_OUTPUT

    def find_output(self, name):
        """Find the output with the given name."""
        for output in self.output_list:
            if output["name"] == name:
                return output
        return None

    def serialize_output(self, result, output):
        """Serialize outputs."""

        name = output["name"]
        if output.get("identifier", "Any").upper() == "WORKGRAPH.NAMESPACE":
            if isinstance(result, dict):
                serialized_result = {}
                for key, value in result.items():
                    full_name = f"{name}.{key}"
                    full_name_output = self.find_output(full_name)
                    if (
                        full_name_output
                        and full_name_output.get("identifier", "Any").upper()
                        == "WORKGRAPH.NAMESPACE"
                    ):
                        serialized_result[key] = self.serialize_output(
                            value, full_name_output
                        )
                    else:
                        serialized_result[key] = general_serializer(value)
                return serialized_result
            else:
                self.exit_codes.ERROR_INVALID_OUTPUT
        else:
            return general_serializer(result)
