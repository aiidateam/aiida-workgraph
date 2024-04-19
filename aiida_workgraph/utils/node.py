from aiida_workgraph.utils import get_executor


def serialize_item(data):
    # print(f"serialize_item: {data}")
    Executor, _ = get_executor(data["serialize"])
    data["value"] = Executor(data["value"])
    return data


def deserialize_item(data):
    # print("deserialize_item: ", data)
    Executor, _ = get_executor(data["deserialize"])
    data["value"] = Executor(data["value"])
    return data
