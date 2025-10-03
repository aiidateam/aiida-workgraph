from typing import Any, Dict
from aiida_workgraph.task import ChildTaskSet, Task
from aiida_workgraph import task, namespace, meta
from node_graph.nodes.builtins import _GraphIOSharedMixin
from node_graph import RuntimeExecutor
from aiida import orm
from node_graph.node_spec import NodeSpec
from node_graph.socket_spec import SocketSpec, SocketSpecMeta
from typing import Annotated
from aiida_workgraph.executors.builtins import get_item, update_ctx, get_context, select


class GraphLevelTask(_GraphIOSharedMixin, Task):
    """Graph level task variant with shared IO."""

    _default_spec = NodeSpec(
        identifier='workgraph.graph_level_task',
        catalog='Builtins',
        base_class_path='aiida_workgraph.tasks.builtins.GraphLevelTask',
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._unify_io()


class Zone(Task):
    """
    Extend the Task class to include a 'children' attribute.
    """

    _default_spec = NodeSpec(
        identifier='workgraph.zone',
        node_type='ZONE',
        catalog='Control',
        base_class_path='aiida_workgraph.tasks.builtins.Zone',
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.children = ChildTaskSet(parent=self)

    def add_task(self, *args, **kwargs):
        """Syntactic sugar to add a task to the zone."""
        task = self.graph.add_task(*args, **kwargs)
        self.children.add(task)
        task.parent = self
        return task

    def to_dict(self, **kwargs) -> Dict[str, Any]:
        tdata = super().to_dict(**kwargs)
        tdata['children'] = [task.name for task in self.children]
        return tdata


class While(Zone):
    """While"""

    _default_spec = NodeSpec(
        identifier='workgraph.while_zone',
        node_type='WHILE',
        catalog='Control',
        inputs=namespace(
            max_iterations=Annotated[int, SocketSpec('workgraph.any', default=10000)],
            conditions=Annotated[Any, SocketSpec('workgraph.any', link_limit=100000)],
        ),
        base_class_path='aiida_workgraph.tasks.builtins.While',
    )


class If(Zone):
    """If task"""

    _default_spec = NodeSpec(
        identifier='workgraph.if_zone',
        node_type='IF',
        catalog='Control',
        inputs=namespace(
            invert_condition=Annotated[bool, SocketSpec('workgraph.bool', default=False)],
            conditions=Annotated[Any, SocketSpec('workgraph.any', link_limit=100000)],
        ),
        base_class_path='aiida_workgraph.tasks.builtins.If',
    )


class Map(Zone):
    """Map"""

    _default_spec = NodeSpec(
        identifier='workgraph.map_zone',
        node_type='MAP',
        catalog='Control',
        inputs=namespace(
            source=Annotated[Any, SocketSpec('workgraph.any', link_limit=100000)],
        ),
        base_class_path='aiida_workgraph.tasks.builtins.Map',
    )

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


class MapItem(Task):
    """MapItem"""

    _default_spec = NodeSpec(
        identifier='workgraph.map_item',
        node_type='Normal',
        catalog='Control',
        inputs=namespace(
            source=SocketSpec('workgraph.any', link_limit=100000, meta=SocketSpecMeta(required=False)),
            key=SocketSpec('workgraph.string', meta=SocketSpecMeta(required=False)),
        ),
        outputs=namespace(item=SocketSpec('workgraph.any')),
        executor=RuntimeExecutor.from_callable(get_item),
        base_class_path='aiida_workgraph.tasks.builtins.MapItem',
    )


class SetContext(Task):
    """SetContext"""

    _default_spec = NodeSpec(
        identifier='workgraph.set_context',
        node_type='Normal',
        catalog='Control',
        inputs=namespace(
            context=SocketSpec('workgraph.any', meta=SocketSpecMeta(required=False)),
            key=any,
            value=any,
        ),
        executor=RuntimeExecutor.from_callable(update_ctx),
        base_class_path='aiida_workgraph.tasks.builtins.SetContext',
    )


class GetContext(Task):
    """GetContext"""

    _default_spec = NodeSpec(
        identifier='workgraph.get_context',
        node_type='Normal',
        catalog='Control',
        inputs=namespace(context=SocketSpec('workgraph.any', meta=SocketSpecMeta(required=False)), key=any),
        outputs=namespace(result=any),
        executor=RuntimeExecutor.from_callable(get_context),
        base_class_path='aiida_workgraph.tasks.builtins.GetContext',
    )


class Select(Task):
    """Select"""

    _default_spec = NodeSpec(
        identifier='workgraph.select',
        node_type='Normal',
        catalog='Control',
        inputs=namespace(
            condition=any,
            true=any,
            false=any,
        ),
        outputs=namespace(result=any),
        executor=RuntimeExecutor.from_callable(select),
        base_class_path='aiida_workgraph.tasks.builtins.Select',
    )


@task(identifier='workgraph.aiida_int')
def aiida_int(value: int) -> orm.Int:
    return orm.Int(value)


@task(identifier='workgraph.aiida_float')
def aiida_float(value: float) -> orm.Float:
    return orm.Float(value)


@task(identifier='workgraph.aiida_string')
def aiida_string(value: str) -> orm.Str:
    return orm.Str(value)


@task(identifier='workgraph.aiida_list')
def aiida_list(value: list) -> orm.List:
    return orm.List(value)


@task(identifier='workgraph.aiida_dict')
def aiida_dict(value: dict) -> orm.Dict:
    return orm.Dict(value)


class AiiDANode(Task):
    """AiiDANode"""

    identifier = 'workgraph.load_node'
    name = 'AiiDANode'
    catalog = 'Test'

    _default_spec = NodeSpec(
        identifier=identifier,
        node_type='Normal',
        inputs=namespace(
            pk=Annotated[int, meta(required=False)],
            uuid=Annotated[str, meta(required=False)],
        ),
        outputs=namespace(code=orm.Code),
        executor=RuntimeExecutor.from_callable(orm.load_node),
        base_class_path='aiida_workgraph.task.Task',
    )


class AiiDACode(Task):
    """AiiDACode"""

    identifier = 'workgraph.load_code'
    name = 'AiiDACode'
    catalog = 'Test'

    _default_spec = NodeSpec(
        identifier=identifier,
        node_type='Normal',
        inputs=namespace(
            pk=Annotated[int, meta(required=False)],
            uuid=Annotated[str, meta(required=False)],
            label=Annotated[str, meta(required=False)],
        ),
        outputs=namespace(code=orm.Code),
        executor=RuntimeExecutor.from_callable(orm.load_code),
        base_class_path='aiida_workgraph.task.Task',
    )
