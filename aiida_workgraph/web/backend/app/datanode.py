from typing import List, Dict, Any
from fastapi import APIRouter, HTTPException, Query
from aiida import orm

router = APIRouter()


@router.get("/api/datanode-data")
async def read_datanode_data(
    typeSearch: str = Query(None),
    labelSearch: str = Query(None),
) -> List[Dict[str, Any]]:
    from aiida.orm import QueryBuilder, Data
    from aiida_workgraph.web.backend.app.utils import time_ago

    try:
        builder = QueryBuilder()
        filters = {}

        if typeSearch:
            filters["node_type"] = {"like": f"%{typeSearch}%"}

        if labelSearch:
            filters["label"] = {"like": f"%{labelSearch}%"}

        builder.append(
            Data,
            filters=filters,
            project=["id", "uuid", "ctime", "node_type", "label"],
            tag="data",
        )
        builder.order_by({"data": {"ctime": "desc"}})
        records = builder.all()
        data = [
            {
                "pk": pk,
                "uuid": uuid,
                "ctime": time_ago(ctime),
                "node_type": node_type,
                "label": label,
            }
            for pk, uuid, ctime, node_type, label in records
        ]
        return data
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Data node {id} not found")


@router.get("/api/datanode/{id}")
async def read_data_node_item(id: int) -> Dict[str, Any]:

    try:
        node = orm.load_node(id)
        content = node.backend_entity.attributes
        content["node_type"] = node.node_type
        return content
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Data node {id} not found")


# Route for deleting a datanode item
@router.delete("/api/datanode/delete/{id}")
async def delete_data_node(
    id: int,
) -> Dict[str, str]:
    from aiida.tools import delete_nodes

    try:
        # Perform the delete action here
        _, was_deleted = delete_nodes([id], dry_run=False)

        if was_deleted:
            return {"message": f"Deleted datanode {id}"}
        else:
            return {"message": f"Failed to delete datanode {id}"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
