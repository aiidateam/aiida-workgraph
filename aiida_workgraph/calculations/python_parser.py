"""Parser for an `PythonJob` job."""
from aiida.parsers.parser import Parser
from aiida_workgraph.orm import general_serializer


class PythonParser(Parser):
    """Parser for an `PythonJob` job."""

    def parse(self, **kwargs):
        """Parse the contents of the output files stored in the `retrieved` output node."""
        import pickle

        try:
            with self.retrieved.base.repository.open("results.pickle", "rb") as handle:
                results = pickle.load(handle)
                output_info = self.node.inputs.output_info.get_list()
                # output_info exclude ['_wait', '_outputs', 'remote_folder', 'remote_stash', 'retrieved']
                outputs = [
                    data
                    for data in output_info
                    if data["name"]
                    not in [
                        "_wait",
                        "_outputs",
                        "remote_folder",
                        "remote_stash",
                        "retrieved",
                    ]
                ]
                if isinstance(results, tuple):
                    if len(outputs) != len(results):
                        raise ValueError(
                            "The number of results does not match the number of outputs."
                        )
                    for i in range(len(outputs)):
                        outputs[i]["value"] = self.serialize_output(
                            results[i], outputs[i]["identifier"]
                        )
                elif isinstance(results, dict) and len(outputs) > 1:
                    for output in outputs:
                        if output.get("required", False):
                            if output["name"] not in results:
                                self.exit_codes.ERROR_MISSING_OUTPUT
                        output["value"] = self.serialize_output(
                            results.pop(output["name"]), output["identifier"]
                        )
                    # if there are any remaining results, raise an warning
                    if results:
                        self.logger.warning(
                            f"Found extra results that are not included in the output: {results.keys()}"
                        )
                elif isinstance(results, dict) and len(outputs) == 1:
                    # if output name in results, use it
                    if outputs[0]["name"] in results:
                        outputs[0]["value"] = self.serialize_output(
                            results[outputs[0]["name"]], outputs[0]["identifier"]
                        )
                    # otherwise, we assume the results is the output
                    else:
                        outputs[0]["value"] = self.serialize_output(
                            results, outputs[0]["identifier"]
                        )
                elif len(outputs) == 1:
                    # otherwise, we assume the results is the output
                    outputs[0]["value"] = self.serialize_output(
                        results, outputs[0]["identifier"]
                    )
                else:
                    raise ValueError(
                        "The number of results does not match the number of outputs."
                    )
                for output in outputs:
                    self.out(output["name"], output["value"])
        except OSError:
            return self.exit_codes.ERROR_READING_OUTPUT_FILE
        except ValueError:
            return self.exit_codes.ERROR_INVALID_OUTPUT

    def serialize_output(self, result, identifier):
        """Serialize outputs."""
        if identifier.upper() == "NAMESPACE":
            if isinstance(result, dict):
                output = {}
                for key, value in result.items():
                    output[key] = general_serializer(value)
                return output
            else:
                self.exit_codes.ERROR_INVALID_OUTPUT
        else:
            return general_serializer(result)
