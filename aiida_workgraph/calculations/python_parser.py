"""Parser for an `PythonCalculation` job."""
from aiida.parsers.parser import Parser


class PythonParser(Parser):
    """Parser for an `PythonCalculation` job."""

    def parse(self, **kwargs):
        """Parse the contents of the output files stored in the `retrieved` output node."""
        import pickle
        from aiida_ase.data.general import GeneralData

        try:
            with self.retrieved.base.repository.open("results.pickle", "rb") as handle:
                results = pickle.load(handle)
                outputs = self.node.inputs.outputs.value
                if isinstance(results, tuple):
                    if len(outputs) != len(results):
                        raise ValueError(
                            "The number of results does not match the number of outputs."
                        )
                    for i in range(len(outputs)):
                        self.out(outputs[i].name, GeneralData(results[i]))
                elif isinstance(results, dict):
                    for key, value in results.items():
                        self.out(key, GeneralData(value))
                else:
                    self.out("result", GeneralData(results))
        except OSError:
            return self.exit_codes.ERROR_READING_OUTPUT_FILE
        except ValueError:
            return self.exit_codes.ERROR_INVALID_OUTPUT
