from fastapi import APIRouter, HTTPException, Query
from aiida import orm

router = APIRouter()


@router.get("/api/datanode-data")
async def read_datanode_data(search: str = Query(None)):
    from aiida.orm import QueryBuilder, Data
    from aiida_worktree.web.backend.app.utils import time_ago

    try:
        builder = QueryBuilder()
        builder.append(
            Data,
            filters={"node_type": {"like": f"%{search}%"}} if search else None,
            project=["id", "uuid", "ctime", "node_type", "label"],
            tag="data",
        )

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
async def read_data_node_item(id: int):

    try:
        node = orm.load_node(id)
        content = {}
        if isinstance(node, orm.Int):
            content["node_type"] = "Int"
            content["value"] = node.value
        elif isinstance(node, orm.Float):
            content["node_type"] = "Float"
            content["value"] = node.value
        elif isinstance(node, orm.Str):
            content["node_type"] = "Str"
            content["value"] = node.value
        elif isinstance(node, orm.Bool):
            content["node_type"] = "Bool"
            content["value"] = node.value
        elif isinstance(node, orm.StructureData):
            content["node_type"] = "StructureData"
            atoms = node.get_ase()
            content["cell"] = atoms.cell.flatten().tolist()
            content["pbc"] = atoms.pbc.tolist()
            content["positions"] = atoms.positions.tolist()
            numbers = atoms.get_atomic_numbers()
            symbols = atoms.get_chemical_symbols()
            content["species"] = {}
            content["speciesArray"] = symbols
            for i in range(len(numbers)):
                if symbols[i] not in content["species"]:
                    content["species"][symbols[i]] = {
                        "atomicNumber": int(numbers[i]),
                        "symbol": symbols[i],
                    }
        return content
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Data node {id} not found")


# Route for deleting a datanode item
@router.delete("/api/datanode/delete/{id}")
async def delete_data_node(
    id: int,
):
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
