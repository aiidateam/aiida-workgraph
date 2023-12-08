from aiida_worktree.node import Node


class AiiDAKpoint(Node):
    identifier = "AiiDAKpoint"
    name = "AiiDAKpoint"
    node_type = "data"
    catalog = "Test"

    kwargs = ["mesh", "offset"]

    def create_properties(self):
        self.properties.new("AiiDAIntVector", "mesh", default=[1, 1, 1], size=3)
        self.properties.new("AiiDAIntVector", "offset", default=[0, 0, 0], size=3)

    def create_sockets(self):
        self.outputs.new("General", "Kpoint")

    def get_executor(self):
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

    def create_properties(self):
        self.properties.new("BaseList", "cell", default=[])
        self.properties.new("BaseList", "kinds", default=[])
        self.properties.new("Bool", "pbc1", default=True)
        self.properties.new("Bool", "pbc2", default=True)
        self.properties.new("Bool", "pbc3", default=True)
        self.properties.new("BaseList", "sites", default=[])

    def create_sockets(self):
        self.outputs.new("General", "Structure")

    def get_executor(self):
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

    def create_properties(self):
        self.properties.new(
            "AiiDAString", "psuedo_familay", default="SSSP/1.2/PBEsol/efficiency"
        )

    def create_sockets(self):
        self.inputs.new("General", "structure")
        self.outputs.new("General", "Pseudo")

    def get_executor(self):
        return {
            "path": "aiida_worktree.executors.qe",
            "name": "get_pseudo_from_structure",
            "type": "function",
        }


class AiiDAPW(Node):
    """DFT calculation using PW code."""

    identifier = "AiiDAPW"
    name = "PW"
    node_type = "calcfunction"
    catalog = "QE"
    args = []
    kwargs = ["kpoints", "parameters", "pseudos", "structure", "code", "metadata"]

    def create_properties(self):
        self.properties.new("BaseDict", "metadata", default={})

    def create_sockets(self):
        self.inputs.new("General", "parameters")
        self.inputs.new("General", "pseudos")
        self.inputs.new("General", "code")
        self.inputs.new("General", "structure")
        self.inputs.new("General", "kpoints")
        self.outputs.new("General", "output_parameters")
        self.outputs.new("General", "remote_folder")
        self.outputs.new("General", "output_structure")
        self.outputs.new("General", "output_trajectory")
        self.outputs.new("General", "output_band")

    def get_executor(self):
        return {
            "name": "quantumespresso.pw",
            "type": "CalculationFactory",
        }


class AiiDADos(Node):
    """DFT calculation using dos code."""

    identifier = "AiiDADos"
    name = "Dos"
    node_type = "calcfunction"
    catalog = "QE"
    args = []
    kwargs = ["parent_folder", "code", "parameters", "metadata"]

    def create_properties(self):
        self.properties.new("BaseDict", "metadata", default={})

    def create_sockets(self):
        self.inputs.new("General", "parent_folder")
        self.inputs.new("General", "code")
        inp = self.inputs.new("General", "parameters")
        inp.add_property("AiiDADict", "parameters", default={})
        self.outputs.new("General", "output_dos")
        self.outputs.new("General", "output_parameters")
        self.outputs.new("General", "remote_folder")

    def get_executor(self):
        return {
            "name": "quantumespresso.dos",
            "type": "CalculationFactory",
        }


class AiiDAProjwfc(Node):
    """DFT calculation using Projwfc code."""

    identifier = "AiiDAProjwfc"
    name = "Projwfc"
    node_type = "calcfunction"
    catalog = "QE"
    args = []
    kwargs = ["parent_folder", "code", "parameters", "metadata"]

    def create_properties(self):
        self.properties.new("BaseDict", "metadata", default={})

    def create_sockets(self):
        self.inputs.new("General", "parent_folder")
        self.inputs.new("General", "code")
        inp = self.inputs.new("General", "parameters")
        inp.add_property("AiiDADict", "parameters", default={})
        self.outputs.new("General", "Dos")
        self.outputs.new("General", "output_parameters")
        self.outputs.new("General", "remote_folder")

    def get_executor(self):
        return {
            "name": "quantumespresso.projwfc",
            "type": "CalculationFactory",
        }
