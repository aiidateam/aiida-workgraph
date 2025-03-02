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
