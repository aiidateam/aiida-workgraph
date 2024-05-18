"""Parser for an `PythonCalculation` job."""
from aiida.parsers.parser import Parser


class PythonParser(Parser):
    """Parser for an `PythonCalculation` job."""

    def parse(self, **kwargs):
        """Parse the contents of the output files stored in the `retrieved` output node."""
        import pickle
        from aiida_ase.data.general import GeneralData

        try:
            with self.retrieved.base.repository.open("result.pickle", "rb") as handle:
                result = pickle.load(handle)
        except OSError:
            return self.exit_codes.ERROR_READING_OUTPUT_FILE
        except ValueError:
            return self.exit_codes.ERROR_INVALID_OUTPUT
        self.out("results", GeneralData(result))
