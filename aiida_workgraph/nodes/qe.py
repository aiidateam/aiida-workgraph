from typing import Dict
from aiida_workgraph.node import Node


class AiiDAKpoint(Node):
    identifier = "AiiDAKpoint"
    name = "AiiDAKpoint"
    node_type = "data"
    catalog = "Test"

    kwargs = ["mesh", "offset"]

    def create_properties(self) -> None:
        self.properties.new("AiiDAIntVector", "mesh", default=[1, 1, 1], size=3)
        self.properties.new("AiiDAIntVector", "offset", default=[0, 0, 0], size=3)

    def create_sockets(self) -> None:
        self.outputs.new("General", "Kpoint")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "KpointsData",
        }


class AiiDAStructure(Node):
    identifier = "AiiDAStructure"
    name = "AiiDAStructure"
    node_type = "data"
    catalog = "Test"

    kwargs = ["cell", "kinds", "pbc1", "pbc2", "pbc3", "sites"]

    def create_properties(self) -> None:
        self.properties.new("BaseList", "cell", default=[])
        self.properties.new("BaseList", "kinds", default=[])
        self.properties.new("Bool", "pbc1", default=True)
        self.properties.new("Bool", "pbc2", default=True)
        self.properties.new("Bool", "pbc3", default=True)
        self.properties.new("BaseList", "sites", default=[])

    def create_sockets(self) -> None:
        self.outputs.new("General", "Structure")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida.orm",
            "name": "StructureData",
        }


class AiiDAPWPseudo(Node):
    identifier = "AiiDAPWPseudo"
    name = "AiiDAPWPseudo"
    node_type = "Normal"
    catalog = "Test"

    args = ["psuedo_familay", "structure"]

    def create_properties(self) -> None:
        self.properties.new(
            "AiiDAString", "psuedo_familay", default="SSSP/1.3/PBEsol/efficiency"
        )

    def create_sockets(self) -> None:
        self.inputs.new("General", "structure")
        self.outputs.new("General", "Pseudo")

    def get_executor(self) -> Dict[str, str]:
        return {
            "path": "aiida_workgraph.executors.qe",
            "name": "get_pseudo_from_structure",
            "type": "function",
        }
