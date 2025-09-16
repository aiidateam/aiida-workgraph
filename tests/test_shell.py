import pytest
from aiida_workgraph import WorkGraph, task, shelljob
from aiida_shell.launch import prepare_code
from aiida.orm import SinglefileData, load_computer, Int


@pytest.mark.skip(reason='need rewrite')
def test_nonexistent_command():
    """Check that the `ValueError` raised by `aiida-shell` for a non-extistent executable is captured by WorkGraph."""

    with pytest.raises(ValueError, match='failed to determine the absolute path'):
        wg = WorkGraph(name='test_shell_command')
        wg.add_task(shelljob, command='abc42')
        wg.run()


def test_shell_command(fixture_localhost):
    """Test the ShellJob with command as a string."""
    wg = WorkGraph(name='test_shell_command')
    job1 = wg.add_task(
        shelljob,
        command='cat',
        resolve_command=True,
        arguments=['{file_a}', '{file_b}'],
        nodes={
            'file_a': SinglefileData.from_string('string a'),
            'file_b': SinglefileData.from_string('string b'),
        },
    )
    # also check if we can set the computer explicitly
    job1.set_inputs({'metadata.computer': load_computer('localhost')})
    wg.run()
    assert job1.outputs.stdout.value.get_content() == 'string astring b'


def test_shell_code():
    """Test the ShellJob with code."""
    cat_code = prepare_code('cat')
    with WorkGraph(name='test_shell_code') as wg:
        # use the code object directly
        outputs = shelljob(
            command=cat_code,
            arguments=['{file_a}', '{file_b}'],
            nodes={
                'file_a': SinglefileData.from_string('string a'),
                'file_b': SinglefileData.from_string('string b'),
            },
        )
        wg.run()
        assert outputs.stdout.value.get_content() == 'string astring b'


def test_dynamic_port():
    """Set the nodes during/after the creation of the task."""
    wg = WorkGraph(name='test_dynamic_port')
    echo_task = wg.add_task(
        shelljob,
        name='echo',
        command='cp',
        arguments=['{file}', 'copied_file'],
        nodes={'file': SinglefileData.from_string('1 5 1')},
        outputs=['copied_file'],
    )

    cat_task = wg.add_task(
        shelljob,
        name='cat',
        command='cat',
        arguments=['{input}'],
        nodes={'input1': None, 'input2': Int(2), 'input3': echo_task.outputs['_wait']},
    )
    wg.add_link(echo_task.outputs['copied_file'], cat_task.inputs['nodes.input1'])
    # task will create input for each item in the dynamic port (nodes)
    assert 'nodes.input1' in cat_task.inputs
    assert 'nodes.input2' in cat_task.inputs
    # if the value of the item is a Socket, then it will create a link, and pop the item
    assert 'nodes.input3' in cat_task.inputs
    assert cat_task.inputs['nodes']._value == {'input2': Int(2)}


def test_outputs_with_dot():
    """Test the outputs with dot in the name."""
    wg = WorkGraph(name='test_outputs_with_dot')
    job1 = wg.add_task(shelljob, command='cat', resolve_command=False, outputs=['file.txt'])

    assert 'file_txt' in job1.outputs
    assert 'file.txt' not in job1.outputs


def test_shell_graph_task():
    """Test the ShellJob inside a graph task.
    And the parser is also defined in the graph task."""
    from aiida.orm import Int

    @task.graph()
    def add_multiply(x, y, command1: str, command2: str) -> Int:
        """Add two numbers and multiply the result by 2."""

        # define the parser function
        def parser(dirpath):
            from aiida.orm import Int

            return {'result': Int((dirpath / 'stdout').read_text().strip())}

        # echo x + y expression
        out1 = shelljob(
            command=command1,
            arguments=['{x}', '+', '{y}'],
            nodes={
                'x': x,
                'y': y,
            },
        )
        # bc command to calculate the expression
        out2 = shelljob(
            command=command2,
            arguments=['{expression}'],
            nodes={'expression': out1.stdout},
            parser=parser,
            parser_outputs=['result'],  # add a "result" output socket from the parser
        )
        return out2.result

    wg = add_multiply.build(x=Int(2), y=Int(3), command1='echo', command2='bc')
    wg.run()
    # wg.submit(wait=True, timeout=60)
    assert wg.outputs.result.value.value == 5
