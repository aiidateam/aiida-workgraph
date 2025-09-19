from typing import Any, Dict
from aiida_workgraph.task import Task, ChildTaskSet, SpecTask
from node_graph.nodes.builtins import _GraphIOSharedMixin
from node_graph import RuntimeExecutor


class GraphLevelTask(_GraphIOSharedMixin, SpecTask):
    """Graph level task variant with shared IO."""

    catalog = 'Builtins'
    is_dynamic: bool = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._unify_io()


class Zone(Task):
    """
    Extend the Task class to include a 'children' attribute.
    """

    identifier = 'workgraph.zone'
    name = 'Zone'
    node_type = 'ZONE'
    catalog = 'Control'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.children = ChildTaskSet(parent=self)

    def add_task(self, *args, **kwargs):
        """Syntactic sugar to add a task to the zone."""
        task = self.graph.add_task(*args, **kwargs)
        self.children.add(task)
        task.parent = self
        return task

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_output('workgraph.any', '_wait')

    def to_dict(self, **kwargs) -> Dict[str, Any]:
        tdata = super().to_dict(**kwargs)
        tdata['children'] = [task.name for task in self.children]
        return tdata


class While(Zone):
    """While"""

    identifier = 'workgraph.while_zone'
    name = 'While'
    node_type = 'WHILE'
    catalog = 'Control'

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_input('workgraph.int', 'max_iterations', property={'default': 10000})
        self.add_input('workgraph.any', 'conditions', link_limit=100000)
        self.add_output('workgraph.any', '_wait')


class If(Zone):
    """If task"""

    identifier = 'workgraph.if_zone'
    name = 'If'
    node_type = 'IF'
    catalog = 'Control'

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_input('workgraph.any', 'conditions')
        self.add_input('workgraph.bool', 'invert_condition', property={'default': False})
        self.add_output('workgraph.any', '_wait')


class Map(Zone):
    """Map"""

    identifier = 'workgraph.map_zone'
    name = 'Map'
    node_type = 'MAP'
    catalog = 'Control'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def item(self):
        for child in self.children:
            if child.identifier == 'workgraph.map_item':
                return child.outputs.item
        # create a child map_item_task if it does not exist
        map_item_task = self.add_task('workgraph.map_item')
        return map_item_task.outputs.item

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.any', 'source', link_limit=100000)
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_output('workgraph.any', '_wait')


class MapItem(Task):
    """MapItem"""

    identifier = 'workgraph.map_item'
    name = 'MapItem'
    node_type = 'Normal'
    catalog = 'Control'

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.any', 'source', link_limit=100000)
        self.add_input('workgraph.any', 'key')
        self.add_output('workgraph.any', 'item')
        self.add_output('workgraph.any', '_wait')

    def get_executor(self):
        from aiida_workgraph.executors.builtins import get_item

        return RuntimeExecutor.from_callable(get_item)


class SetContext(Task):
    """SetContext"""

    identifier = 'workgraph.set_context'
    name = 'SetContext'
    node_type = 'Normal'
    catalog = 'Control'

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.any', 'context')
        self.add_input('workgraph.any', 'key')
        self.add_input('workgraph.any', 'value')
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_output('workgraph.any', '_wait')

    def get_executor(self):
        from aiida_workgraph.executors.builtins import update_ctx

        return RuntimeExecutor.from_callable(update_ctx)


class GetContext(Task):
    """GetContext"""

    identifier = 'workgraph.get_context'
    name = 'GetContext'
    node_type = 'Normal'
    catalog = 'Control'

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.any', 'context')
        self.add_input('workgraph.any', 'key')
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_output('workgraph.any', 'result')
        self.add_output('workgraph.any', '_wait')

    def get_executor(self):
        from aiida_workgraph.executors.builtins import get_context

        return RuntimeExecutor.from_callable(get_context)


class AiiDAInt(Task):
    identifier = 'workgraph.aiida_int'
    name = 'AiiDAInt'
    node_type = 'Normal'
    catalog = 'Test'

    def update_sockets(self) -> None:
        self.add_input('workgraph.any', 'value', property={'default': 0.0})
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_output('workgraph.int', 'result')
        self.add_output('workgraph.any', '_wait')

    def get_executor(self):
        from aiida import orm

        return RuntimeExecutor.from_callable(orm.Int)


class AiiDAFloat(Task):
    identifier = 'workgraph.aiida_float'
    name = 'AiiDAFloat'
    node_type = 'Normal'
    catalog = 'Test'

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.float', 'value', property={'default': 0.0})
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_output('workgraph.float', 'result')
        self.add_output('workgraph.any', '_wait')

    def get_executor(self):
        from aiida import orm

        return RuntimeExecutor.from_callable(orm.Float)


class AiiDAString(Task):
    identifier = 'workgraph.aiida_string'
    name = 'AiiDAString'
    node_type = 'Normal'
    catalog = 'Test'

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.string', 'value', property={'default': ''})
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_output('workgraph.string', 'result')
        self.add_output('workgraph.any', '_wait')

    def get_executor(self):
        from aiida import orm

        return RuntimeExecutor.from_callable(orm.Str)


class AiiDAList(Task):
    identifier = 'workgraph.aiida_list'
    name = 'AiiDAList'
    node_type = 'Normal'
    catalog = 'Test'

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.any', 'value', property={'default': []})
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_output('workgraph.list', 'result')
        self.add_output('workgraph.any', '_wait')

    def get_executor(self):
        from aiida import orm

        return RuntimeExecutor.from_callable(orm.List)


class AiiDADict(Task):
    identifier = 'workgraph.aiida_dict'
    name = 'AiiDADict'
    node_type = 'Normal'
    catalog = 'Test'

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.any', 'value', property={'default': {}})
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_output('workgraph.dict', 'result')
        self.add_output('workgraph.any', '_wait')

    def get_executor(self):
        from aiida import orm

        return RuntimeExecutor.from_callable(orm.Dict)


class AiiDANode(Task):
    """AiiDANode"""

    identifier = 'workgraph.load_node'
    name = 'AiiDANode'
    node_type = 'Normal'
    catalog = 'Test'

    def create_properties(self) -> None:
        pass

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.any', 'identifier')
        self.add_input('workgraph.any', 'pk')
        self.add_input('workgraph.any', 'uuid')
        self.add_input('workgraph.any', 'label')
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_output('workgraph.any', 'node')
        self.add_output('workgraph.any', '_wait')

    def get_executor(self):
        from aiida import orm

        return RuntimeExecutor.from_callable(orm.load_node)


class AiiDACode(Task):
    """AiiDACode"""

    identifier = 'workgraph.load_code'
    name = 'AiiDACode'
    node_type = 'Normal'
    catalog = 'Test'

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.any', 'identifier')
        self.add_input('workgraph.any', 'pk')
        self.add_input('workgraph.any', 'uuid')
        self.add_input('workgraph.any', 'label')
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_output('workgraph.any', 'Code')
        self.add_output('workgraph.any', '_wait')

    def get_executor(self):
        from aiida import orm

        return RuntimeExecutor.from_callable(orm.load_code)


class Select(Task):
    """Select"""

    identifier = 'workgraph.select'
    name = 'Select'
    node_type = 'Normal'
    catalog = 'Control'

    def update_sockets(self) -> None:
        self.inputs._clear()
        self.outputs._clear()
        self.add_input('workgraph.any', 'condition')
        self.add_input('workgraph.any', 'true')
        self.add_input('workgraph.any', 'false')
        self.add_input('workgraph.any', '_wait', link_limit=100000, metadata={'arg_type': 'none'})
        self.add_output('workgraph.any', 'result')
        self.add_output('workgraph.any', '_wait')

    def get_executor(self):
        from aiida_workgraph.executors.builtins import select

        return RuntimeExecutor.from_callable(select)
