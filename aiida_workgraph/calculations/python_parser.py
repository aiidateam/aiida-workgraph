"""Parser for an `PythonJob` job."""
from aiida.parsers.parser import Parser
from .general_data import GeneralData


class PythonParser(Parser):
    """Parser for an `PythonJob` job."""

    def parse(self, **kwargs):
        """Parse the contents of the output files stored in the `retrieved` output node."""
        import pickle

        try:
            with self.retrieved.base.repository.open("results.pickle", "rb") as handle:
                results = pickle.load(handle)
                output_name_list = self.node.inputs.output_name_list.get_list()
                if isinstance(results, tuple):
                    if len(output_name_list) != len(results):
                        raise ValueError(
                            "The number of results does not match the number of output_name_list."
                        )
                    for i in range(len(output_name_list)):
                        self.out(output_name_list[i].name, GeneralData(results[i]))
                elif isinstance(results, dict):
                    for key, value in results.items():
                        self.out(key, GeneralData(value))
                else:
                    self.out("result", GeneralData(results))
        except OSError:
            return self.exit_codes.ERROR_READING_OUTPUT_FILE
        except ValueError:
            return self.exit_codes.ERROR_INVALID_OUTPUT
