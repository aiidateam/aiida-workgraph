import pathlib
from setuptools import setup, find_packages


def test_suite():
    import unittest

    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover("tests", pattern="test_*.py")
    return test_suite


# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name="aiida-worktree",
    version="0.0.3",
    description="Design flexible node-based workflow for AiiDA calculation.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/superstart54/aiida-worktree",
    author="Xing Wang",
    author_email="xingwang1991@gmail.com",
    license="MIT License",
    classifiers=[],
    packages=find_packages(),
    install_requires=[
        "numpy",
        "aiida-core",
        "scinode",
        "cloudpickle",
        "aiida-pseudo",
        "aiida-quantumespresso",
        "pytest",
        "pre-commit",
    ],
    entry_points={
        "aiida.node": [
            "process.workflow.worktree = aiida_worktree.orm.worktree:WorkTreeNode",
        ],
        "scinode_node": [
            "aiida = aiida_worktree.nodes:node_list",
        ],
        "scinode_property": [
            "aiida = aiida_worktree.properties.built_in:property_list",
        ],
    },
    package_data={},
    python_requires=">=3.6",
    test_suite="setup.test_suite",
)
