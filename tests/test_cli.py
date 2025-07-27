from click.testing import CliRunner
from aiida_workgraph.cli.cmd_workgraph import workgraph


def test_workgraph():
    """Test ``verdi group path ls``"""
    cli_runner = CliRunner()
    result = cli_runner.invoke(workgraph, ["--help"])
    assert result.exit_code == 0, result.exception
    print(result.output)


def test_graph():
    """Test ``verdi group path ls``"""
    cli_runner = CliRunner()
    result = cli_runner.invoke(workgraph, ["graph", "--help"])
    assert result.exit_code == 0, result.exception
    print(result.output)


def test_task():
    """Test ``verdi group path ls``"""
    cli_runner = CliRunner()
    result = cli_runner.invoke(workgraph, ["task", "--help"])
    assert result.exit_code == 0, result.exception
    print(result.output)
