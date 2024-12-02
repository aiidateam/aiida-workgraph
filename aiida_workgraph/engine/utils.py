from aiida import orm
from aiida.common.extendeddicts import AttributeDict


def prepare_for_workgraph_task(task: dict, kwargs: dict) -> tuple:
    """Prepare the inputs for WorkGraph task"""
    from aiida_workgraph.utils import organize_nested_inputs, serialize_properties
    from aiida.orm.utils.serialize import deserialize_unsafe

    wgdata = deserialize_unsafe(task["executor"]["wgdata"])
    wgdata["name"] = task["name"]
    wgdata["metadata"]["group_outputs"] = task["metadata"]["group_outputs"]
    # update the workgraph data by kwargs
    for task_name, data in kwargs.items():
        # because kwargs is updated using update_nested_dict_with_special_keys
        # which means the data is grouped by the task name
        for socket_name, value in data.items():
            wgdata["tasks"][task_name]["inputs"][socket_name]["property"][
                "value"
            ] = value
    # merge the properties
    organize_nested_inputs(wgdata)
    serialize_properties(wgdata)
    metadata = {"call_link_label": task["name"]}
    inputs = {"wg": wgdata, "metadata": metadata}
    return inputs, wgdata


def sort_socket_data(socket_data: dict) -> dict:
    """Sort the socket data by the list_index"""
    data = [
        {"name": data["name"], "identifier": data["identifier"]}
        for data, _ in sorted(
            ((data, data["list_index"]) for data in socket_data.values()),
            key=lambda x: x[1],
        )
    ]
    return data


def prepare_for_python_task(task: dict, kwargs: dict, var_kwargs: dict) -> dict:
    """Prepare the inputs for PythonJob"""
    from aiida_pythonjob import prepare_pythonjob_inputs

    # get the names kwargs for the PythonJob, which are the inputs before _wait
    function_inputs = kwargs.pop("function_inputs", {})
    sorted_inputs = sort_socket_data(task["inputs"])
    # TODO better way to find the function_inputs
    for input in sorted_inputs:
        if input["name"] == "_wait":
            break
        function_inputs[input["name"]] = kwargs.pop(input["name"], None)

    # if the var_kwargs is not None, we need to pop the var_kwargs from the kwargs
    # then update the function_inputs if var_kwargs is not None
    if task["var_kwargs"] is not None:
        function_inputs.pop(task["var_kwargs"], None)
        if var_kwargs:
            # var_kwargs can be AttributeDict if it get data from the previous task output
            if isinstance(var_kwargs, (dict, AttributeDict)):
                function_inputs.update(var_kwargs)
            # otherwise, it should be a Data node
            elif isinstance(var_kwargs, orm.Data):
                function_inputs.update(var_kwargs.value)
            else:
                raise ValueError(f"Invalid var_kwargs type: {type(var_kwargs)}")
    # setup code
    code = kwargs.pop("code", None)
    computer = kwargs.pop("computer", "localhost")
    command_info = kwargs.pop("command_info", {})
    upload_files = kwargs.pop("upload_files", {})

    metadata = kwargs.pop("metadata", {})
    metadata.update({"call_link_label": task["name"]})
    # get the source code of the function
    executor = task["executor"]

    # outputs
    outputs = [
        {"name": output["name"], "identifier": output["identifier"]}
        for output, _ in sorted(
            ((output, output["list_index"]) for output in task["outputs"].values()),
            key=lambda x: x[1],
        )
    ]
    # only the output before _wait is the function_outputs
    function_outputs = []
    for output in outputs:
        if output["name"] == "_wait":
            break
        # if the output is WORKGRAPH.NAMESPACE, we need to change it to NAMESPACE
        if output["identifier"].upper() == "WORKGRAPH.NAMESPACE":
            function_outputs.append({"name": output["name"], "identifier": "NAMESPACE"})
        else:
            function_outputs.append(output)

    inputs = prepare_pythonjob_inputs(
        pickled_function=executor,
        function_inputs=function_inputs,
        function_outputs=function_outputs,
        code=code,
        command_info=command_info,
        computer=computer,
        metadata=metadata,
        upload_files=upload_files,
        process_label=f"PythonJob<{task['name']}>",
        **kwargs,
    )

    return inputs


def prepare_for_shell_task(task: dict, inputs: dict) -> dict:
    """Prepare the inputs for ShellJob"""
    from aiida_shell.launch import prepare_shell_job_inputs
    import inspect

    # Retrieve the signature of `prepare_shell_job_inputs` to determine expected input parameters.
    signature = inspect.signature(prepare_shell_job_inputs)
    aiida_shell_input_keys = signature.parameters.keys()

    # Iterate over all WorkGraph `inputs`, and extract the ones which are expected by `prepare_shell_job_inputs`
    inputs_aiida_shell_subset = {
        key: inputs[key] for key in inputs.keys() if key in aiida_shell_input_keys
    }

    try:
        aiida_shell_inputs = prepare_shell_job_inputs(**inputs_aiida_shell_subset)
    except ValueError:
        raise

    # We need to remove the original input-keys, as they might be offending for the call to `launch_shell_job`
    # E.g., `inputs` originally can contain `command`, which gets, however, transformed to #
    # `code` by `prepare_shell_job_inputs`
    for key in inputs_aiida_shell_subset.keys():
        inputs.pop(key)

    # Finally, we update the original `inputs` with the modified ones from the call to `prepare_shell_job_inputs`
    inputs = {**inputs, **aiida_shell_inputs}

    inputs.setdefault("metadata", {})
    inputs["metadata"].update({"call_link_label": task["name"]})
    return inputs
