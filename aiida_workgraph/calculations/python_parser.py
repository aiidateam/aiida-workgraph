"""Parser for an `PythonJob` job."""
from aiida.parsers.parser import Parser
from aiida_workgraph.orm import serialize_to_aiida_nodes


class PythonParser(Parser):
    """Parser for an `PythonJob` job."""

    def parse(self, **kwargs):
        """Parse the contents of the output files stored in the `retrieved` output node."""
        import pickle

        try:
            with self.retrieved.base.repository.open("results.pickle", "rb") as handle:
                results = pickle.load(handle)
                output_name_list = self.node.inputs.output_name_list.get_list()
                # output_name_list exclude ['_wait', '_outputs', 'remote_folder', 'remote_stash', 'retrieved']
                output_name_list = [
                    name
                    for name in output_name_list
                    if name
                    not in [
                        "_wait",
                        "_outputs",
                        "remote_folder",
                        "remote_stash",
                        "retrieved",
                    ]
                ]
                outputs = {}
                if isinstance(results, tuple):
                    if len(output_name_list) != len(results):
                        raise ValueError(
                            "The number of results does not match the number of output_name_list."
                        )
                    for i in range(len(output_name_list)):
                        outputs[output_name_list[i].name] = results[i]
                        outputs = serialize_to_aiida_nodes(outputs)
                elif isinstance(results, dict) and len(results) == len(
                    output_name_list
                ):
                    outputs = serialize_to_aiida_nodes(results)
                else:
                    outputs = serialize_to_aiida_nodes({"result": results})
                for key, value in outputs.items():
                    self.out(key, value)
        except OSError:
            return self.exit_codes.ERROR_READING_OUTPUT_FILE
        except ValueError:
            return self.exit_codes.ERROR_INVALID_OUTPUT
