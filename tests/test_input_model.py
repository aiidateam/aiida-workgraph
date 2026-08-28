"""A Pydantic model as a task's wire contract, under the AiiDA engine.

The contract itself -- which sockets a model declares, and when its rules are
held to -- is node-graph's, and its tests live there. What is asserted here is
what only this package can answer: that the model decides the *stored* form of
a value, that the body gets the rich object back out of storage, and that a
rule broken inside a submitted process fails that process with a message
naming the task.
"""

from __future__ import annotations

import enum
from decimal import Decimal
from typing import Any

import pytest
from aiida import orm
from node_graph.input_model import BODY_RECEIVES, ModelContractError, TaskInputValidationError
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    field_serializer,
    field_validator,
    model_validator,
)

from aiida_workgraph import WorkGraph, task
from aiida_workgraph.socket_spec import spec_from_model

#: aiida-pythonjob's exit status for a body that raised.
FUNCTION_FAILED = 323


class Color(enum.Enum):
    RED = 'red'
    BLUE = 'blue'


# --------------------------------------------------------------------------
# 1. The model builds this package's sockets
# --------------------------------------------------------------------------


class AddInputs(BaseModel):
    """Two summands, the second optional."""

    x: int
    y: int = 7


@task(input_model=AddInputs)
def add(x, y):
    return x + y


def test_the_sockets_carry_this_packages_identifiers():
    """The spec is built through aiida-workgraph's socket vocabulary, not node-graph's."""
    fields = add._spec.inputs.fields
    assert fields['x'].identifier == 'workgraph.int'
    assert fields['x'].meta.required is True
    assert fields['y'].identifier == 'workgraph.int'
    assert fields['y'].meta.required is False
    assert fields['y'].default == 7


def test_an_omitted_input_runs_on_the_models_default():
    wg = WorkGraph('add_default')
    node = wg.add_task(add, name='add', x=2)
    wg.run()
    assert node.outputs.result.value.value == 9


class NudgeInputs(BaseModel):
    """``by`` defaults to ``None``, the value the engine declines to store."""

    x: int
    by: int | None = None


@task(input_model=NudgeInputs)
def nudge(x, by):
    return x if by is None else x + by


def test_a_field_defaulting_to_none_is_not_a_missing_required_input():
    assert nudge._spec.inputs.fields['by'].meta.required is False
    wg = WorkGraph('nudge')
    node = wg.add_task(nudge, name='nudge', x=4)
    wg.run()
    assert node.outputs.result.value.value == 4


def test_a_model_on_a_process_task_is_refused():
    from aiida.calculations.arithmetic.add import ArithmeticAddCalculation

    with pytest.raises(ModelContractError, match='plain Python function task'):
        task(input_model=AddInputs)(ArithmeticAddCalculation)


# --------------------------------------------------------------------------
# 2. An Enum: stored as its bare value, delivered as the member
# --------------------------------------------------------------------------


class PaintInputs(BaseModel):
    color: Color


@task(input_model=PaintInputs, outputs=['is_member', 'seen'])
def paint(color):
    # ``is`` and not ``==``: the body gets the member itself, not a look-alike.
    return {'is_member': color is Color.RED, 'seen': type(color).__name__}


@task(inputs=spec_from_model(PaintInputs), outputs=['is_member', 'seen'])
def paint_without_model(color):
    """Same sockets, no contract: the control for :func:`paint`."""
    return {'is_member': color is Color.RED, 'seen': type(color).__name__}


def test_an_enum_is_stored_as_its_bare_value():
    wg = WorkGraph('paint_store')
    node = wg.add_task(paint, name='paint', color=Color.RED)
    wg.run()
    stored = dict(node.process.inputs.function_inputs)['color']
    assert stored.value == 'red'


def test_the_body_receives_the_enum_member_and_only_with_the_model():
    wg = WorkGraph('paint_body')
    with_model = wg.add_task(paint, name='with_model', color=Color.RED)
    without_model = wg.add_task(paint_without_model, name='without_model', color=Color.RED)
    wg.run()

    assert with_model.outputs.is_member.value.value is True
    assert with_model.outputs.seen.value.value == 'Color'
    # The control proves the model is what rebuilt the member: the same socket
    # hands the bare value to a body no model stands in front of.
    assert without_model.outputs.is_member.value.value is False
    assert without_model.outputs.seen.value.value == 'str'


def test_an_enum_survives_a_workgraph_round_trip():
    wg = WorkGraph('paint_round_trip')
    wg.add_task(paint, name='paint', color=Color.BLUE)
    rebuilt = WorkGraph.from_dict(wg.to_dict())
    assert rebuilt.tasks.paint.inputs.color.value == Color.BLUE
    rebuilt.run()
    assert rebuilt.tasks.paint.outputs.is_member.value.value is False
    assert rebuilt.tasks.paint.outputs.seen.value.value == 'Color'


def test_no_class_path_is_written_into_the_task():
    """The spec stored with the task names no class: the model is reached through the code."""
    assert 'structured_type' not in paint._spec.inputs.fields['color'].meta.extras


# --------------------------------------------------------------------------
# 3. A type JSON cannot hold, carried by the model's own serializer
# --------------------------------------------------------------------------


class MoneyInputs(BaseModel):
    """A ``Decimal`` amount, stored as the string the model renders."""

    amount: Decimal

    @field_serializer('amount')
    def _dump_amount(self, value: Decimal) -> str:
        return str(value)

    @field_validator('amount', mode='before')
    @classmethod
    def _load_amount(cls, value):
        return value if isinstance(value, Decimal) else Decimal(str(value))


@task(input_model=MoneyInputs, outputs=['kind', 'doubled'])
def double_money(amount):
    return {'kind': type(amount).__name__, 'doubled': str(amount * 2)}


@task(inputs=spec_from_model(MoneyInputs), outputs=['kind', 'doubled'])
def double_money_without_model(amount):
    """Same sockets, no contract: the control for :func:`double_money`."""
    return {'kind': type(amount).__name__, 'doubled': str(amount * 2)}


def test_a_field_serializer_decides_the_stored_form():
    wg = WorkGraph('money_store')
    node = wg.add_task(double_money, name='money', amount=Decimal('0.10'))
    wg.run()
    stored = dict(node.process.inputs.function_inputs)['amount']
    assert stored.value == '0.10'


def test_the_body_receives_the_decimal_and_only_with_the_model():
    wg = WorkGraph('money_body')
    node = wg.add_task(double_money, name='money', amount=Decimal('0.10'))
    wg.run()
    assert node.outputs.kind.value.value == 'Decimal'
    # Exact, because a Decimal round-tripped as a string never became a float.
    assert node.outputs.doubled.value.value == '0.20'


def test_without_the_model_the_same_value_cannot_even_be_stored():
    """The control: nothing else in the stack knows how to write a ``Decimal``."""
    wg = WorkGraph('money_control')
    wg.add_task(double_money_without_model, name='money', amount=Decimal('0.10'))
    with pytest.raises(ValueError, match='decimal.Decimal'):
        wg.run()


# --------------------------------------------------------------------------
# 4. Rules the socket layer cannot see fail the process that broke them
# --------------------------------------------------------------------------


class RangeInputs(BaseModel):
    """``low`` and ``high`` are ints the socket accepts; their order is the model's rule."""

    low: int
    high: int = Field(le=100)

    @model_validator(mode='after')
    def _ordered(self):
        if self.low >= self.high:
            raise ValueError('low must be below high')
        return self


@task(input_model=RangeInputs)
def span(low, high):
    return high - low


def test_a_cross_field_rule_fails_the_process_naming_the_task():
    wg = WorkGraph('span_bad_order')
    node = wg.add_task(span, name='span', low=9, high=3)
    wg.run()
    assert node.process.exit_status == FUNCTION_FAILED
    assert "Task 'span' got inputs RangeInputs rejects" in node.process.exit_message
    assert 'low must be below high' in node.process.exit_message


def test_a_field_constraint_fails_where_it_is_written():
    """``le=100`` is a model rule the socket layer cannot see, and ``add_task`` still refuses 500."""
    wg = WorkGraph('span_too_high')
    with pytest.raises(TaskInputValidationError, match='less than or equal to 100'):
        wg.add_task(span, name='span', low=1, high=500)


def test_without_checkpoint_a_that_constraint_reaches_the_run_edge(monkeypatch):
    """The control: with the write unchecked, ``le=100`` is first seen where the body runs."""
    from node_graph import input_model

    monkeypatch.setattr(input_model, 'validate_task_inputs', lambda task, inputs: None)
    wg = WorkGraph('span_too_high_unchecked')
    node = wg.add_task(span, name='span', low=1, high=500)
    wg.run()
    assert node.process.exit_status == FUNCTION_FAILED
    assert 'less than or equal to 100' in node.process.exit_message


# --------------------------------------------------------------------------
# 5. A graph task's contract, held where the engine expands it
# --------------------------------------------------------------------------


class WindowInputs(BaseModel):
    """The window a graph opens; ``lower`` below ``upper`` is the graph's own rule."""

    lower: int
    upper: int

    @model_validator(mode='after')
    def _ordered(self):
        if self.lower >= self.upper:
            raise ValueError('lower must be below upper')
        return self


@task.graph(input_model=WindowInputs)
def window(lower, upper):
    return add(x=lower, y=upper)


@task()
def shrink(value):
    return value - 10


@task.graph()
def window_of_a_computed_bound(lower, upper):
    """The subgraph's bound is decided by a task, so it is a value only at run time."""
    shrunk = shrink(value=upper)
    return window(lower=lower, upper=shrunk.result)


def test_a_graph_contract_holds_when_the_graph_is_submitted():
    wg = WorkGraph('window_ok')
    node = wg.add_task(window, name='window', lower=1, upper=3)
    wg.run()
    assert node.process.exit_status == 0


def test_a_graph_contract_fails_the_submitted_graph():
    """The graph is refused as it is expanded, so it never becomes a process at all."""
    wg = WorkGraph('window_bad')
    node = wg.add_task(window, name='window', lower=9, upper=3)
    wg.run()
    assert node.process is None
    assert node.state == 'FAILED'
    assert wg.process.exit_status != 0


class Named(BaseModel):
    """A `str` field, which is what the engine's own wrappers get in the way of."""

    label: str
    count: int


@task()
def echo(label, count):
    return f'{label}-{count}'


@task.graph(input_model=Named)
def named(label, count):
    return echo(label=label, count=count)


def test_a_graph_contract_reads_through_the_engines_wrappers():
    """A graph body is handed storage nodes; the contract is checked against what they hold."""
    wg = WorkGraph('named_ok')
    node = wg.add_task(named, name='named', label='silicon', count=2)
    wg.run()
    assert node.process.exit_status == 0


class Priced(BaseModel):
    """A `Decimal` beside a `str`: one kind no socket identifier can carry."""

    label: str
    amount: Decimal


@task()
def show(label, amount):
    return f'{label}-{amount}'


@task.graph(input_model=Priced)
def priced(label, amount):
    return show(label=label, amount=str(amount))


def test_a_graph_contract_reads_a_field_no_identifier_carries():
    """`Decimal` is stored as the string the model rendered, and read back through the model."""
    wg = WorkGraph('priced_ok')
    node = wg.add_task(priced, name='priced', label='silicon', amount=Decimal('1.50'))
    wg.run()
    assert node.process.exit_status == 0


def test_without_the_unwrap_that_field_would_be_refused(monkeypatch):
    """The control: leave the node on and `Decimal` refuses the `orm.Str` holding its rendering.

    The field is deliberately one a socket identifier cannot answer for -- it is
    `workgraph.annotated`, not `workgraph.string` -- so the control still
    discriminates under an adapter that unwraps by identifier.
    """
    from node_graph.serializer import SerializationAdapter

    from aiida_workgraph.serialization import AiidaSerializationAdapter

    monkeypatch.setattr(AiidaSerializationAdapter, 'deserialize', SerializationAdapter.deserialize)
    wg = WorkGraph('priced_unread')
    node = wg.add_task(priced, name='priced', label='silicon', amount=Decimal('1.50'))
    wg.run()
    assert node.process is None
    assert node.state == 'FAILED'


ARRIVED: dict = {}


def _arrival(value):
    """Return the type of what a body was handed, seeing through the tag it wears."""
    return type(getattr(value, '__wrapped__', value)).__name__


class EveryKind(BaseModel):
    """One field per kind the read edge has to tell apart."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    text: str
    number: int
    fraction: float
    flag: bool
    mapping: dict
    items: list
    anything: Any
    node: orm.Int


@task.graph(input_model=EveryKind)
def records_arrivals(text, number, fraction, flag, mapping, items, anything, node):
    ARRIVED.clear()
    ARRIVED.update(
        text=_arrival(text),
        number=_arrival(number),
        fraction=_arrival(fraction),
        flag=_arrival(flag),
        mapping=_arrival(mapping),
        items=_arrival(items),
        anything=_arrival(anything),
        node=_arrival(node),
    )


def test_a_body_receives_what_its_field_declares():
    """Declared Python arrives as Python; `Any` and an AiiDA type arrive as nodes."""
    wg = WorkGraph('arrivals')
    wg.add_task(
        records_arrivals,
        name='r',
        text='silicon',
        number=3,
        fraction=1.5,
        flag=True,
        mapping={'k': 1},
        items=[1, 2],
        anything='whatever',
        node=orm.Int(7),
    )
    wg.run()
    assert ARRIVED == {
        'text': 'str',
        'number': 'int',
        'fraction': 'float',
        'flag': 'bool',
        'mapping': 'dict',
        'items': 'list',
        'anything': 'Str',
        'node': 'Int',
    }


def test_the_socket_identifier_alone_cannot_tell_int_from_orm_int():
    """Both are ``workgraph.int``; only the model says which of the two was declared."""
    from aiida_workgraph.socket_spec import spec_from_model

    spec = spec_from_model(EveryKind)
    assert spec.fields['number'].identifier == spec.fields['node'].identifier
    assert spec.fields['number'].meta.extras[BODY_RECEIVES] == 'python'
    assert spec.fields['node'].meta.extras[BODY_RECEIVES] == 'node'


def _naive_unwrap(value):
    """Unwrap the way an identifier-keyed read edge does, handing the tag back to no one."""
    from node_graph.socket import TaggedValue

    plain = value.__wrapped__ if isinstance(value, TaggedValue) else value
    if isinstance(plain, orm.BaseType):
        return plain.value
    if isinstance(plain, orm.Dict):
        return plain.get_dict()
    if isinstance(plain, orm.List):
        return plain.get_list()
    return plain


def _label_nodes(node):
    """Return the node written into the graph's ``label`` and the one its subtask got."""
    outer = inner = None
    for link in node.process.base.links.get_incoming().all():
        if link.link_label.endswith('label'):
            outer = link.node
    for child in node.process.called:
        for link in child.base.links.get_incoming().all():
            if link.link_label.endswith('label'):
                inner = link.node
    return outer, inner


def test_an_unwrapped_value_keeps_the_tag_that_draws_its_link():
    """The body forwards a `str` field, and the child must be linked to it, not given a copy.

    Provenance is what says which: the node the subtask reads has to be the
    very node the graph was given, not a second one holding the same string.
    """
    wg = WorkGraph('linked')
    node = wg.add_task(named, name='named', label='silicon', count=2)
    wg.run()
    assert node.process.exit_status == 0
    outer, inner = _label_nodes(node)
    assert outer is not None and inner is not None
    assert outer.uuid == inner.uuid


def test_without_the_tag_the_subtask_reads_a_node_nobody_produced(monkeypatch):
    """The control: unwrap without retagging and the body holds a copy, not a reference."""
    import aiida_workgraph.serialization as serialization

    monkeypatch.setattr(serialization, '_to_declared_python', _naive_unwrap)
    wg = WorkGraph('linked_naive')
    node = wg.add_task(named, name='named', label='silicon', count=2)
    wg.run()
    outer, inner = _label_nodes(node)
    assert outer is not None and inner is not None
    assert outer.uuid != inner.uuid
    assert len(inner.base.links.get_incoming().all()) == 0


def test_a_value_that_needed_no_unwrapping_keeps_its_tag_too():
    """A plain `str` has no node to take off, and dropping the tag there drops the link."""
    from node_graph.socket import TaggedValue

    from aiida_workgraph.serialization import _to_declared_python

    tagged = TaggedValue('silicon', socket=object())
    unwrapped = _to_declared_python(tagged)
    assert isinstance(unwrapped, TaggedValue)
    assert unwrapped._socket is tagged._socket
    assert unwrapped._uuid == tagged._uuid


def test_a_value_that_was_unwrapped_keeps_the_uuid_it_arrived_with():
    """One value, one uuid: a new one would make the body's value a second value."""
    from node_graph.socket import TaggedValue

    from aiida_workgraph.serialization import _to_declared_python

    tagged = TaggedValue(orm.Int(7), socket=object())
    unwrapped = _to_declared_python(tagged)
    assert unwrapped == 7
    assert unwrapped._uuid == tagged._uuid


def test_a_runtime_value_is_checked_at_the_graph_it_reaches():
    """Nothing knows ``upper`` is 5 until ``shrink`` has run, so this is the first chance."""
    wg = WorkGraph('window_computed')
    node = wg.add_task(window_of_a_computed_bound, name='outer', lower=9, upper=15)
    wg.run()
    assert node.process.exit_status != 0


class Spin(enum.Enum):
    """The four spin treatments a workflow accepts."""

    NONE = 'none'
    COLLINEAR = 'collinear'
    NON_COLLINEAR = 'non_collinear'
    SPIN_ORBIT = 'spin_orbit'


class DielectricInputs(BaseModel):
    """Two of the four members are ph.x's own limit, written where ph.x is declared."""

    spin: Spin = Spin.NONE
    structure: str

    @field_validator('spin')
    @classmethod
    def _supported(cls, value):
        if value in (Spin.NON_COLLINEAR, Spin.SPIN_ORBIT):
            raise ValueError('ph.x has no electric-field perturbation for noncollinear magnetism')
        return value


@task(input_model=DielectricInputs, outputs=['seen'])
def dielectric(spin, structure):
    return {'seen': type(spin).__name__}


class EverySpin(BaseModel):
    spin: Spin = Spin.NONE
    structure: str


@task.graph(input_model=EverySpin)
def eps(spin, structure):
    return dielectric(spin=spin, structure=structure).seen


def test_a_rule_on_an_inner_task_fails_the_expansion():
    """The graph takes every spin, so the value meets ph.x's rule as the body wires it."""
    wg = WorkGraph('eps_noncollinear')
    node = wg.add_task(eps, name='eps', spin=Spin.NON_COLLINEAR, structure='si')
    wg.run()
    assert node.state == 'FAILED'
    assert node.process is None
    # Neither the subgraph nor the task it would have held became a process.
    assert [called.process_label for called in wg.process.called_descendants] == []


def test_the_spin_the_rule_admits_reaches_the_body_as_the_member():
    """The control: the same wiring expands, runs, and the enum survives storage."""
    wg = WorkGraph('eps_collinear')
    node = wg.add_task(eps, name='eps', spin=Spin.COLLINEAR, structure='si')
    wg.run()
    assert node.process.exit_status == 0
    assert sorted(called.process_label for called in wg.process.called_descendants) == [
        'WorkGraph<eps>',
        'dielectric',
    ]
    assert node.outputs.result.value.value == 'Spin'


# --------------------------------------------------------------------------
# 6. Output models
# --------------------------------------------------------------------------


class SumAndProduct(BaseModel):
    sum: int
    product: int


@task(output_model=SumAndProduct)
def combine(x, y):
    return {'sum': x + y, 'product': x * y}


@task(output_model=SumAndProduct)
def combine_forgetting_product(x, y):
    return {'sum': x + y}


@task(output_model=SumAndProduct)
def combine_with_a_bad_type(x, y):
    return {'sum': 'not a number', 'product': x * y}


def test_the_output_sockets_come_from_the_output_model():
    assert set(combine._spec.outputs.fields) == {'sum', 'product'}


def test_a_return_the_model_accepts_lands_on_the_sockets():
    wg = WorkGraph('combine_ok')
    node = wg.add_task(combine, name='combine', x=2, y=3)
    wg.run()
    assert node.outputs.sum.value.value == 5
    assert node.outputs.product.value.value == 6


def test_a_missing_output_fails_at_the_source_task():
    wg = WorkGraph('combine_missing')
    node = wg.add_task(combine_forgetting_product, name='combine', x=2, y=3)
    wg.run()
    assert node.process.exit_status == FUNCTION_FAILED
    assert "Task 'combine_forgetting_product' returned outputs SumAndProduct rejects" in node.process.exit_message
    assert 'product' in node.process.exit_message


def test_a_mistyped_output_fails_at_the_source_task():
    wg = WorkGraph('combine_bad_type')
    node = wg.add_task(combine_with_a_bad_type, name='combine', x=2, y=3)
    wg.run()
    assert node.process.exit_status == FUNCTION_FAILED
    assert 'sum' in node.process.exit_message


# --------------------------------------------------------------------------
# 7. Mappings whose size is only known at runtime
# --------------------------------------------------------------------------


class Block(BaseModel):
    """One member of a mapping the task fills in at runtime."""

    width: int
    label: str


class Blocks(BaseModel):
    blocks: dict[str, Block]


class Recipe(BaseModel):
    """Two mappings on the way in: one of plain strings, one needing the model."""

    pseudos: dict[str, str]
    amounts: dict[str, Decimal] = {}

    @field_serializer('amounts')
    def _dump_amounts(self, value: dict[str, Decimal]) -> dict[str, str]:
        return {key: str(amount) for key, amount in value.items()}

    @field_validator('amounts', mode='before')
    @classmethod
    def _load_amounts(cls, value):
        return {key: amount if isinstance(amount, Decimal) else Decimal(str(amount)) for key, amount in value.items()}


@task(output_model=Blocks)
def cut_blocks(n):
    """Decide how many members the mapping has, inside the body."""
    return {'blocks': {f'b{i}': {'width': i + 1, 'label': f'block {i}'} for i in range(n)}}


@task(output_model=Blocks)
def cut_one_bad_block(n):
    return {'blocks': {'b0': {'width': 1, 'label': 'fine'}, 'b1': {'width': 'wide', 'label': 'bad'}}}


@task(input_model=Blocks, outputs=['widest'])
def widest_block(blocks):
    # Each member arrived through ``Block`` as the keys that were written.
    return {'widest': max(item['width'] for item in blocks.values())}


@task(input_model=Recipe, outputs=['names', 'kinds'])
def describe_recipe(pseudos, amounts):
    return {
        'names': ','.join(sorted(pseudos)),
        'kinds': ','.join(sorted({type(amount).__name__ for amount in amounts.values()})),
    }


def test_a_typed_mapping_field_becomes_a_dynamic_namespace():
    blocks = cut_blocks._spec.outputs.fields['blocks']
    assert blocks.is_namespace()
    assert blocks.meta.dynamic is True
    # Every member is typed, so a key's sockets are known before any key is.
    assert set(blocks.item.fields) == {'width', 'label'}
    assert blocks.item.fields['width'].identifier == 'workgraph.int'


def test_a_mapping_becomes_one_socket_per_key_after_the_run():
    wg = WorkGraph('blocks_out')
    node = wg.add_task(cut_blocks, name='cut', n=3)
    wg.run()
    assert set(node.outputs.blocks._sockets) == {'b0', 'b1', 'b2'}
    assert node.outputs.blocks.b1.width.value.value == 2
    assert node.outputs.blocks.b1.label.value.value == 'block 1'


def test_a_bad_member_fails_at_the_source_task_naming_its_key():
    wg = WorkGraph('blocks_bad')
    node = wg.add_task(cut_one_bad_block, name='cut', n=2)
    wg.run()
    assert node.process.exit_status == FUNCTION_FAILED
    # The key is in the path pydantic reports, so the reader knows which member.
    assert 'blocks.b1.width' in node.process.exit_message


def test_a_downstream_task_consumes_the_mapping():
    wg = WorkGraph('blocks_chain')
    source = wg.add_task(cut_blocks, name='cut', n=3)
    consumer = wg.add_task(widest_block, name='widest', blocks=source.outputs.blocks)
    wg.run()
    assert consumer.outputs.widest.value.value == 3


def test_addressing_one_future_member_at_build_time_is_not_available_yet():
    """The per-key sockets appear when the task runs, so a build-time link cannot name one.

    Ordered and by-name member access on a socket that has not produced its
    members yet is scinode/node-graph#160.
    """
    wg = WorkGraph('blocks_future')
    source = wg.add_task(cut_blocks, name='cut', n=3)
    with pytest.raises(AttributeError, match="has no sub-socket 'b1'"):
        source.outputs.blocks.b1


def test_a_mapping_input_becomes_one_socket_per_key():
    wg = WorkGraph('recipe_in')
    node = wg.add_task(describe_recipe, name='recipe', pseudos={'Si': 'si.upf', 'O': 'o.upf'})
    assert set(node.inputs.pseudos._sockets) == {'Si', 'O'}
    wg.run()
    assert node.outputs.names.value.value == 'O,Si'


def test_a_member_of_a_mapping_is_stored_through_the_model():
    wg = WorkGraph('recipe_decimal')
    node = wg.add_task(
        describe_recipe,
        name='recipe',
        pseudos={'Si': 'si.upf'},
        amounts={'Si': Decimal('1.50'), 'O': Decimal('0.25')},
    )
    wg.run()
    stored = dict(node.process.inputs.function_inputs)['amounts']
    assert {key: value.value for key, value in dict(stored).items()} == {'Si': '1.50', 'O': '0.25'}
    assert node.outputs.kinds.value.value == 'Decimal'


class Amount(BaseModel):
    """A model used as a field, carrying a type JSON cannot hold."""

    value: Decimal

    @field_serializer('value')
    def _dump_value(self, value: Decimal) -> str:
        return str(value)

    @field_validator('value', mode='before')
    @classmethod
    def _load_value(cls, value):
        return value if isinstance(value, Decimal) else Decimal(str(value))


class Nested(BaseModel):
    cfg: Amount


class Mapped(BaseModel):
    cfgs: dict[str, Amount]


class Layer(BaseModel):
    inner: Amount


class DeepMapped(BaseModel):
    items: dict[str, Layer]


@task(input_model=Nested, outputs=['kind', 'stored'])
def read_nested(cfg):
    # A nested model reaches the body as the members that were written.
    return {'kind': type(cfg['value']).__name__, 'stored': str(cfg['value'])}


@task(input_model=Mapped, outputs=['kind', 'stored'])
def read_mapped(cfgs):
    # Each item of a mapping reaches the body as the members written into it.
    return {'kind': type(cfgs['a']['value']).__name__, 'stored': str(cfgs['a']['value'])}


@task(input_model=DeepMapped, outputs=['kind', 'stored'])
def read_deep(items):
    value = items['a']['inner']['value']
    return {'kind': type(value).__name__, 'stored': str(value)}


@pytest.mark.parametrize(
    'entry_point, payload, path',
    [
        (read_nested, {'cfg': {'value': Decimal('1.50')}}, ('cfg', 'value')),
        (read_mapped, {'cfgs': {'a': {'value': Decimal('1.50')}}}, ('cfgs', 'a', 'value')),
        (read_deep, {'items': {'a': {'inner': {'value': Decimal('1.50')}}}}, ('items', 'a', 'inner', 'value')),
    ],
    ids=['nested-model', 'mapping-of-models', 'mapping-of-nested-models'],
)
def test_a_model_renders_its_own_leaf_however_deep_it_sits(entry_point, payload, path):
    """Whatever models and mappings the path crosses, the model declaring the leaf renders it."""
    wg = WorkGraph(f'depth_{path[0]}')
    node = wg.add_task(entry_point, name='read', **payload)
    wg.run()

    assert node.process.exit_status == 0
    stored = node.process.inputs.function_inputs
    for name in path:
        stored = stored[name]
    assert stored.value == '1.50'
    # And the body gets the Decimal back, not the string that was stored.
    assert node.outputs.kind.value.value == 'Decimal'
    assert node.outputs.stored.value.value == '1.50'


# --------------------------------------------------------------------------
# 8. What a task without a model does
# --------------------------------------------------------------------------


@task
def plain_add(x: int, y: int = 7) -> int:
    return x + y


def test_a_task_without_a_model_is_untouched():
    fields = plain_add._spec.inputs.fields
    assert fields['x'].identifier == 'workgraph.int'
    assert fields['y'].default == 7
    wg = WorkGraph('plain')
    node = wg.add_task(plain_add, name='add', x=2)
    wg.run()
    assert node.outputs.result.value.value == 9


# --------------------------------------------------------------------------
# 9. The node a body was promised, and the edge that takes it off
# --------------------------------------------------------------------------


class StructureInputs(BaseModel):
    """One socket the model can only be satisfied by the node itself."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    structure: orm.StructureData


class PayloadInputs(BaseModel):
    """``Any`` declares nothing to rebuild, so its socket arrives as stored."""

    payload: Any = None


@task.calcfunction(input_model=StructureInputs)
def reads_a_structure(structure):
    return orm.Str(type(structure).__name__)


@task.calcfunction(input_model=PayloadInputs)
def reads_a_payload(payload):
    return orm.Str(type(payload).__name__)


@task(outputs=['seen'])
def reads_a_structure_by_annotation(structure: orm.StructureData):
    return {'seen': type(structure).__name__}


def a_silicon_structure() -> orm.StructureData:
    """Return a two-atom silicon cell, stored."""
    structure = orm.StructureData(cell=[[0.0, 2.7, 2.7], [2.7, 0.0, 2.7], [2.7, 2.7, 0.0]])
    structure.append_atom(position=(0.0, 0.0, 0.0), symbols='Si')
    structure.append_atom(position=(1.35, 1.35, 1.35), symbols='Si')
    return structure.store()


def test_the_spec_says_a_node_typed_socket_reaches_the_body_as_the_node():
    """What the contract promises for these two sockets, before anything runs."""
    structure_spec = spec_from_model(StructureInputs)
    payload_spec = spec_from_model(PayloadInputs)
    assert structure_spec.fields['structure'].meta.extras[BODY_RECEIVES] == 'node'
    assert payload_spec.fields['payload'].meta.extras[BODY_RECEIVES] == 'node'


def test_a_node_typed_field_is_refused_on_a_pyfunction():
    """A PyFunction reads its inputs out of their nodes, so this field cannot be met."""
    with pytest.raises(ModelContractError, match='@task.calcfunction'):

        @task(input_model=StructureInputs, outputs=['seen'])
        def reads_it(structure):
            return {'seen': type(structure).__name__}


def test_a_field_under_any_is_refused_on_a_pyfunction_too():
    """``Any`` declares nothing to rebuild, so its socket is node-typed as well."""
    with pytest.raises(ModelContractError, match='@task.calcfunction'):

        @task(input_model=PayloadInputs, outputs=['seen'])
        def reads_it(payload):
            return {'seen': type(payload).__name__}


def test_a_calcfunction_body_is_handed_the_node_the_model_declared():
    """Where the rule sends such a field: the body reads the node itself."""
    wg = WorkGraph('structure_leaf')
    node = wg.add_task(reads_a_structure, name='leaf', structure=a_silicon_structure())
    wg.run()
    assert node.process.exit_status == 0
    assert node.outputs.result.value.value == 'StructureData'


def test_a_calcfunction_is_handed_the_node_under_any_as_well():
    """The field kind no per-type override could cover reaches the body whole."""
    wg = WorkGraph('payload_leaf')
    node = wg.add_task(reads_a_payload, name='leaf', payload=orm.Dict(dict={'a': 1}).store())
    wg.run()
    assert node.process.exit_status == 0
    assert node.outputs.result.value.value == 'Dict'


def test_the_same_socket_declared_without_a_model_loses_the_node_in_silence():
    """The control: the edge behaves the same way with no model to notice it."""
    wg = WorkGraph('structure_annotated')
    node = wg.add_task(reads_a_structure_by_annotation, name='leaf', structure=a_silicon_structure())
    wg.run()
    assert node.process.exit_status == 0
    assert node.outputs.seen.value.value == 'Atoms'


# --------------------------------------------------------------------------
# 10. A graph body is handed the type its field declares
# --------------------------------------------------------------------------


class SpinOnly(BaseModel):
    spin: Spin = Spin.NONE


@task(outputs=['seen'])
def report(text):
    return {'seen': text}


def describe(value) -> str:
    """Return what a body can tell about the value it was handed."""
    from node_graph.socket import TaggedValue

    # By kind, not by class name: which class tags a value is the tagging
    # layer's business, and it has more than one.
    kind = 'tagged' if isinstance(value, TaggedValue) else type(value).__name__
    return f'{kind}|{value.__class__.__name__}|{value == Spin.COLLINEAR}|{value in (Spin.NONE, Spin.COLLINEAR)}'


@task.graph(input_model=SpinOnly)
def modelled_spin(spin):
    return report(text=describe(spin)).seen


@task.graph()
def annotated_spin(spin: Spin = Spin.NONE):
    return report(text=describe(spin)).seen


def test_a_modelled_graph_body_is_handed_the_member_under_its_tag():
    """The member, so ``==`` and ``in`` answer; the tag, so the link is still drawn."""
    wg = WorkGraph('spin_modelled')
    node = wg.add_task(modelled_spin, name='g', spin=Spin.COLLINEAR)
    wg.run()
    assert node.outputs.result.value.value == 'tagged|Spin|True|True'


def test_a_graph_body_declaring_the_same_field_by_annotation_agrees():
    """The control: the two ways of declaring one field hand the body one value."""
    wg = WorkGraph('spin_annotated')
    node = wg.add_task(annotated_spin, name='g', spin=Spin.COLLINEAR)
    wg.run()
    assert node.outputs.result.value.value == 'tagged|Spin|True|True'


# --------------------------------------------------------------------------
# 11. A handle that does not take the name of the function it decorates
# --------------------------------------------------------------------------


class BandWindow(BaseModel):
    """A field rule for the write, and a cross-field rule only the run edge holds."""

    lower: int = 0
    upper: int = 1

    @field_validator('lower')
    @classmethod
    def _counted(cls, value):
        if value < 0:
            raise ValueError('lower counts bands, so it cannot be negative')
        return value

    @model_validator(mode='after')
    def _ordered(self):
        if self.upper <= self.lower:
            raise ValueError('upper must be above lower')
        return self


def _bands_body(lower, upper):
    return upper - lower


#: The decorated name is never rebound, so the module still binds the function.
named_apart_bands = task(input_model=BandWindow)(_bands_body)


@task(input_model=BandWindow)
def rebound_bands(lower, upper):
    """The ordinary spelling: the handle replaces the module global."""
    return upper - lower


SPELLINGS = pytest.mark.parametrize(
    'handle',
    [named_apart_bands, rebound_bands],
    ids=['handle named apart', 'handle rebound to the name'],
)


@SPELLINGS
def test_the_write_is_refused_whatever_the_handle_is_called(handle):
    """The write reads the model off the executor, and both spellings store one."""
    wg = WorkGraph('bands_write')
    with pytest.raises(TaskInputValidationError, match='cannot be negative'):
        wg.add_task(handle, name='w', lower=-5, upper=9)


@SPELLINGS
def test_the_submitted_task_holds_the_rule_the_write_could_not(handle):
    """A cross-field rule waits for the run edge, which the engine reaches the same way."""
    wg = WorkGraph('bands_run')
    node = wg.add_task(handle, name='w', lower=9, upper=3)
    wg.run()
    assert node.process.exit_status == FUNCTION_FAILED
    assert 'upper must be above lower' in node.process.exit_message


@SPELLINGS
def test_the_bounds_the_model_admits_still_run(handle):
    """The control: what the model accepts reaches the body and comes back."""
    wg = WorkGraph('bands_ok')
    node = wg.add_task(handle, name='w', lower=1, upper=9)
    wg.run()
    assert node.process.exit_status == 0
    assert node.outputs.result.value.value == 8


def test_the_process_is_still_labelled_with_the_name_of_its_function():
    """What the executor is stored under is its own business, not the label's."""
    wg = WorkGraph('bands_label')
    node = wg.add_task(named_apart_bands, name='w', lower=1, upper=9)
    wg.run()
    assert node.process.process_label == '_bands_body'


# --------------------------------------------------------------------------
# 12. What a member nobody wrote is worth at this engine's run edge
# --------------------------------------------------------------------------


class EngineSystem(BaseModel):
    """Five members, each with a default the model answers for."""

    nbnd: int = 1
    nosym: bool = False
    ecutwfc: float = 60.0
    occupations: str = 'fixed'
    degauss: float = 0.0


class EngineRoute(BaseModel):
    spin: Spin = Spin.NONE
    system: EngineSystem = EngineSystem()


@task(input_model=EngineRoute, outputs=['report'])
def reports_its_system(spin, system):
    return {'report': f'{spin.value}|{sorted(system)}'}


def test_the_leaf_body_is_handed_the_members_that_were_written():
    """One member written, and the body reads that member alone out of storage.

    The top-level field is the other half of the rule: nobody wrote ``spin``,
    so the model answers for it and the body is handed its default.
    """
    wg = WorkGraph('written_members')
    node = wg.add_task(reports_its_system, name='leaf', system={'nbnd': 20})
    wg.run()
    assert node.process.exit_status == 0
    assert node.outputs.report.value.value == "none|['nbnd']"


class EngineBlock(BaseModel):
    """One member of a mapping, with defaults of its own."""

    num_iter: int = 100
    dis_froz_max: float = 0.0
    num_wann: int = 4


class EngineBlocks(BaseModel):
    blocks: dict[str, EngineBlock] = {}


@task(input_model=EngineBlocks, outputs=['report'])
def reports_its_blocks(blocks):
    return {'report': '|'.join(f'{key}:{sorted(item)}' for key, item in sorted(blocks.items()))}


def test_each_item_of_a_mapping_reaches_the_body_as_what_was_written():
    """A mapping is one namespace per key, so a key carries only its written members."""
    wg = WorkGraph('written_items')
    node = wg.add_task(reports_its_blocks, name='leaf', blocks={'occ_1': {'num_iter': 42}})
    wg.run()
    assert node.process.exit_status == 0
    assert node.outputs.report.value.value == "occ_1:['num_iter']"


# --------------------------------------------------------------------------
# 13. What a calcfunction's parameters can be, and what a broken rule leaves
# --------------------------------------------------------------------------


class NamespaceInputs(BaseModel):
    """A parameter a process function has no port for."""

    system: EngineSystem = EngineSystem()


class WindowInSteps(BaseModel):
    """Two nodes and a rule that cannot answer until both are in hand."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    lower: orm.Int
    upper: orm.Int

    @model_validator(mode='after')
    def _ordered(self):
        if self.upper.value <= self.lower.value:
            raise ValueError('upper must be above lower')
        return self


@task.calcfunction(input_model=WindowInSteps)
def width(lower, upper):
    return orm.Int(upper.value - lower.value)


def test_a_namespace_field_is_refused_on_a_calcfunction():
    """A calcfunction parameter is one port carrying one node, never a namespace."""
    with pytest.raises(ModelContractError, match='one port carrying one node'):

        @task.calcfunction(input_model=NamespaceInputs)
        def reads_it(system):
            return orm.Str(repr(sorted(system)))


def test_the_same_namespace_is_what_a_pyfunction_body_reads():
    """The control, and where the rule sends such a field."""
    wg = WorkGraph('namespace_pyfunction')
    node = wg.add_task(reports_its_system, name='leaf', system={'nbnd': 20})
    wg.run()
    assert node.process.exit_status == 0


def test_a_calcfunction_runs_what_its_rule_admits():
    """The scalar-node shape a calcfunction does carry, rule and all."""
    wg = WorkGraph('width_ok')
    node = wg.add_task(width, name='t', lower=orm.Int(1).store(), upper=orm.Int(9).store())
    wg.run()
    assert node.process.exit_status == 0
    assert node.outputs.result.value.value == 8


def test_a_rule_broken_inside_a_calcfunction_excepts_it():
    """What a calcfunction can carry when a rule fails: an exception, not an exit status.

    A process function has no controlled-failure channel for an exception in
    its body, so the model's report reaches the excepted ``CalcFunctionNode``
    rather than an exit message.
    """
    wg = WorkGraph('width_broken')
    wg.add_task(width, name='t', lower=orm.Int(9).store(), upper=orm.Int(3).store())
    wg.run()
    excepted = (
        orm.QueryBuilder()
        .append(
            orm.CalcFunctionNode,
            tag='n',
            filters={'label': 'width', 'attributes.process_state': 'excepted'},
        )
        .order_by({'n': {'id': 'desc'}})
        .first()[0]
    )
    assert excepted.exit_status is None
    assert 'upper must be above lower' in (excepted.exception or '')
