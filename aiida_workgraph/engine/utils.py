from aiida_workgraph.orm.serializer import serialize_to_aiida_nodes
from aiida import orm


def prepare_for_workgraph_node(node, kwargs):
    """Prepare the inputs for WorkGraph node"""
    from aiida_workgraph.utils import merge_properties

    print("node  type: workgraph.")
    wgdata = node["executor"]["wgdata"]
    wgdata["name"] = node["name"]
    wgdata["metadata"]["group_outputs"] = node["metadata"]["group_outputs"]
    # update the workgraph data by kwargs
    for node_name, data in kwargs.items():
        # because kwargs is updated using update_nested_dict_with_special_keys
        # which means the data is grouped by the node name
        for socket_name, value in data.items():
            wgdata["nodes"][node_name]["properties"][socket_name]["value"] = value
    # merge the properties
    merge_properties(wgdata)
    metadata = {"call_link_label": node["name"]}
    inputs = {"wg": wgdata, "metadata": metadata}
    return inputs, wgdata


def prepare_for_pythonjob(node, executor, args, kwargs, var_args, var_kwargs):
    """Prepare the inputs for PythonJob"""
    from aiida_workgraph.utils import get_or_create_code
    import os

    print("node  type: Python.")
    # get the names kwargs for the PythonJob, which are the inputs before _wait
    input_kwargs = {}
    for input in node["inputs"]:
        if input["name"] == "_wait":
            break
        input_kwargs[input["name"]] = kwargs.pop(input["name"], None)
    # setup code
    code = kwargs.pop("code", None)
    computer = kwargs.pop("computer", None)
    code_label = kwargs.pop("code_label", None)
    code_path = kwargs.pop("code_path", None)
    prepend_text = kwargs.pop("prepend_text", None)
    upload_files = kwargs.pop("upload_files", {})
    new_upload_files = {}
    # change the string in the upload files to SingleFileData, or FolderData
    for key, source in upload_files.items():
        # only alphanumeric and underscores are allowed in the key
        # replace all "." with "_dot_"
        new_key = key.replace(".", "_dot_")
        if isinstance(source, str):
            if os.path.isfile(source):
                new_upload_files[new_key] = orm.SinglefileData(file=source)
            elif os.path.isdir(source):
                new_upload_files[new_key] = orm.FolderData(tree=source)
        elif isinstance(source, (orm.SinglefileData, orm.FolderData)):
            new_upload_files[new_key] = source
        else:
            raise ValueError(f"Invalid upload file type: {type(source)}, {source}")
    #
    if code is None:
        code = get_or_create_code(
            computer=computer if computer else "localhost",
            code_label=code_label if code_label else "python3",
            code_path=code_path if code_path else None,
            prepend_text=prepend_text if prepend_text else None,
        )
    parent_folder = kwargs.pop("parent_folder", None)
    metadata = kwargs.pop("metadata", {})
    metadata.update({"call_link_label": node["name"]})
    # get the source code of the function
    function_name = executor.__name__
    function_source_code = (
        node["executor"]["import_statements"]
        + "\n"
        + node["executor"]["function_source_code"]
    )
    # outputs
    output_name_list = [output["name"] for output in node["outputs"]]
    # serialize the kwargs into AiiDA Data
    input_kwargs = serialize_to_aiida_nodes(input_kwargs)
    # transfer the args to kwargs
    inputs = {
        "function_source_code": orm.Str(function_source_code),
        "function_name": orm.Str(function_name),
        "code": code,
        "kwargs": input_kwargs,
        "upload_files": new_upload_files,
        "output_name_list": orm.List(output_name_list),
        "parent_folder": parent_folder,
        "metadata": metadata,
        **kwargs,
    }
    return inputs
