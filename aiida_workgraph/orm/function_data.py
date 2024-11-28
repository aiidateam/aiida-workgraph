import inspect
import textwrap
from typing import Callable, Dict, Any, get_type_hints, _SpecialForm
from .pickled_data import PickledData


class PickledFunction(PickledData):
    """Data class to represent a pickled Python function."""

    def __init__(self, value=None, **kwargs):
        """Initialize a PickledFunction node instance.

        :param value: a Python function
        """
        super().__init__(**kwargs)
        if not callable(value):
            raise ValueError("value must be a callable Python function")
        self.set_value(value)
        self.set_attribute(value)

    def __str__(self):
        return (
            f"PickledFunction<{self.base.attributes.get('function_name')}> pk={self.pk}"
        )

    @property
    def metadata(self):
        """Return a dictionary of metadata."""
        return {
            "name": self.base.attributes.get("name"),
            "import_statements": self.base.attributes.get("import_statements"),
            "source_code": self.base.attributes.get("source_code"),
            "source_code_without_decorator": self.base.attributes.get(
                "source_code_without_decorator"
            ),
            "type": "function",
            "is_pickle": True,
        }

    @classmethod
    def build_callable(cls, func):
        """Return the executor for this node."""
        import cloudpickle as pickle

        executor = {
            "executor": pickle.dumps(func),
            "type": "function",
            "is_pickle": True,
        }
        executor.update(cls.inspect_function(func))
        return executor

    def set_attribute(self, value):
        """Set the contents of this node by pickling the provided function.

        :param value: The Python function to pickle and store.
        """
        # Serialize the function and extract metadata
        serialized_data = self.inspect_function(value)

        # Store relevant metadata
        self.base.attributes.set("name", serialized_data["name"])
        self.base.attributes.set(
            "import_statements", serialized_data["import_statements"]
        )
        self.base.attributes.set("source_code", serialized_data["source_code"])
        self.base.attributes.set(
            "source_code_without_decorator",
            serialized_data["source_code_without_decorator"],
        )

    @classmethod
    def inspect_function(cls, func: Callable) -> Dict[str, Any]:
        """Serialize a function for storage or transmission."""
        try:
            # we need save the source code explicitly, because in the case of jupyter notebook,
            # the source code is not saved in the pickle file
            source_code = inspect.getsource(func)
            # Split the source into lines for processing
            source_code_lines = source_code.split("\n")
            function_source_code = "\n".join(source_code_lines)
            # Find the first line of the actual function definition
            for i, line in enumerate(source_code_lines):
                if line.strip().startswith("def "):
                    break
            function_source_code_without_decorator = "\n".join(source_code_lines[i:])
            function_source_code_without_decorator = textwrap.dedent(
                function_source_code_without_decorator
            )
            # we also need to include the necessary imports for the types used in the type hints.
            try:
                required_imports = cls.get_required_imports(func)
            except Exception as e:
                required_imports = {}
                print(
                    f"Failed to get required imports for function {func.__name__}: {e}"
                )
            # Generate import statements
            import_statements = "\n".join(
                f"from {module} import {', '.join(types)}"
                for module, types in required_imports.items()
            )
        except Exception as e:
            print(f"Failed to inspect function {func.__name__}: {e}")
            function_source_code = ""
            function_source_code_without_decorator = ""
            import_statements = ""
        return {
            "name": func.__name__,
            "source_code": function_source_code,
            "source_code_without_decorator": function_source_code_without_decorator,
            "import_statements": import_statements,
        }

    @classmethod
    def get_required_imports(cls, func: Callable) -> Dict[str, set]:
        """Retrieve type hints and the corresponding modules."""
        type_hints = get_type_hints(func)
        imports = {}

        def add_imports(type_hint):
            if isinstance(
                type_hint, _SpecialForm
            ):  # Handle special forms like Any, Union, Optional
                module_name = "typing"
                type_name = type_hint._name or str(type_hint)
            elif hasattr(
                type_hint, "__origin__"
            ):  # This checks for higher-order types like List, Dict
                module_name = type_hint.__module__
                type_name = getattr(type_hint, "_name", None) or getattr(
                    type_hint.__origin__, "__name__", None
                )
                for arg in getattr(type_hint, "__args__", []):
                    if arg is type(None):  # noqa: E721
                        continue
                    add_imports(arg)  # Recursively add imports for each argument
            elif hasattr(type_hint, "__module__"):
                module_name = type_hint.__module__
                type_name = type_hint.__name__
            else:
                return  # If no module or origin, we can't import it, e.g., for literals

            if type_name is not None:
                if module_name not in imports:
                    imports[module_name] = set()
                imports[module_name].add(type_name)

        for _, type_hint in type_hints.items():
            add_imports(type_hint)

        return imports


def to_pickled_function(value):
    """Convert a Python function to a `PickledFunction` instance."""
    return PickledFunction(value)


class PickledLocalFunction(PickledFunction):
    """PickledFunction subclass for local functions."""
